import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from itertools import combinations

def main():
    parser = argparse.ArgumentParser(description='Process pseudobulk DESeq2 analysis for a specific tissue.')
    parser.add_argument('--tissue', '-t', type=str, help='Name of the tissue for analysis', required=True)
    
    args = parser.parse_args()
    
    tissue = args.tissue

    file = f'../IGVF_analysis/annotated_tissues/{tissue}_annotated.h5ad'
    outdirectory = f'degs/{tissue}/'

    group = 'Genotype'
    sexes = ['Male','Female']
    sample_key = 'Mouse_Tissue_ID'
    condition_key = 'Genotype'
    groupby = 'celltype'

    adata = sc.read_h5ad(file)
    
    adata = adata[~adata.obs['celltype'].isin(["igvf_003", 
                                               "igvf_004", 
                                               "igvf_005", 
                                               "igvf_008",
                                               "igvf_008b",
                                               "igvf_009",
                                               "igvf_010",
                                               "igvf_011",
                                               "igvf_012",
                                               "Bleedthrough",
                                               "Doublet"])]

    if tissue == "Heart":
        adata = adata[~adata.obs['celltype'].isin(["Schwann","Epithelial"])]


    # sample genotypes for each 'Multiplexed_sampleN'
    ms1 = ['B6J','AJ','WSBJ','129S1J']
    ms2 = ['NODJ','PWKJ','NZOJ','CASTJ']

    # Define a function to update 'Mouse_Tissue_ID' based on conditions
    def update_mouse_tissue_id(row):
        if row['plate'] == 'igvf_010' and row['Column'] in [9.0, 10.0, 11.0, 12.0]:
            if row['Genotype'] in ms1:
                return row['Multiplexed_sample1']
            elif row['Genotype'] in ms2:
                return row['Multiplexed_sample2']
        return row['Mouse_Tissue_ID']


    # Apply the function to update the 'Mouse_Tissue_ID' column
    adata.obs['Mouse_Tissue_ID'] = adata.obs.apply(update_mouse_tissue_id, axis=1)

    ######################################################################################################################
    # Differential Gene analysis loops
    # performs DGA for all celltypes for between all genotype pairs within males and females 

    for sex in sexes:
        filtering = {'Sex': ['Male', 'Female'],
                'Genotype': ['CASTJ', 'B6J', 'AJ', 'PWKJ', 'WSBJ', '129S1J', 'NODJ', 'NZOJ']}

        filtering['Sex'] = [sex]
        genotype_combinations = combinations(filtering['Genotype'], 2)

        for genotype1, genotype2 in genotype_combinations:
            print(tissue)
            print(genotype1, ' ', genotype2, ' ' , sex)
            print('###################')

            filtering[group] = [genotype1, genotype2]

            adatac = adata[adata.obs['Sex'].isin(filtering['Sex'])]
            adatac = adatac[adatac.obs['Genotype'].isin(filtering['Genotype'])]

            adatac.obs[sample_key] = adatac.obs[sample_key].astype(str)

            #####
            # Create pseudobulk 
            pdata = dc.get_pseudobulk(
                adatac,
                sample_col=sample_key,
                groups_col=groupby,
                obs=adatac.obs,
                layer='raw_counts',
                mode='sum',
                min_cells=10,
                min_counts=10000
            )

            print(pdata.shape)

            #for Bug fix below
            unique_df = adatac.obs[['Mouse_Tissue_ID','Genotype']].drop_duplicates(subset=['Mouse_Tissue_ID'])
            unique_dict = unique_df.set_index('Mouse_Tissue_ID')['Genotype'].to_dict()

            dea_results = {}
            for cell_group in pdata.obs[groupby].unique():
                # Select cell type profiles
                print(cell_group)
                ctdata = pdata[pdata.obs[groupby] == cell_group].copy()

                ctdata.obs['Genotype'] = ctdata.obs['Mouse_Tissue_ID'].map(unique_dict)

    #             # Obtain genes that pass the edgeR-like thresholds
    #             genes = dc.filter_by_expr(ctdata,
    #                                       group=condition_key,
    #                                       min_count=5, # a minimum number of counts in a number of samples
    #                                       min_total_count=10 # a minimum total number of reads across samples
    #                                       )
    #             print(f"Filtering out {genes.shape[0]} lowly-expressed genes.")

    #             # Filter by these genes
    #             ctdata = ctdata[:, genes].copy()

                # Build DESeq2 object
                dds = DeseqDataSet(
                    adata=ctdata,
                    design_factors=condition_key,
                    ref_level=[condition_key, genotype2], 
                    refit_cooks=True, # got KeyError on heart schwann... 
                    quiet = True
                )

                # Compute LFCs
                dds.deseq2()

                stat_res = DeseqStats(dds, contrast=[condition_key, genotype1, genotype2])
                # Compute Wald test
                stat_res.summary()

                # Shrink LFCs
                coeff_name = condition_key+'_' + genotype1 + '_vs_' + genotype2
                stat_res.lfc_shrink(coeff= coeff_name) 
                dea_results[cell_group] = stat_res.results_df
                dea_results[cell_group]['comparison'] = coeff_name

                # concat results across cell types
            dea_df = pd.concat(dea_results)

            dea_df = dea_df.reset_index().rename(columns={'level_0': groupby}).set_index('celltype')
            print(dea_df.head())

            dea_df['genotype1'] = genotype1
            dea_df['genotype2'] = genotype2
            dea_df['sex'] = sex

            dea_df.to_csv(outdirectory + genotype1 + '_' + genotype2 + '_' + sex + '_DEG_results.csv')

if __name__ == "__main__":
    main()

