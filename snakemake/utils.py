import pandas as pd
import anndata
import scanpy as sc
import scrublet as scr

from bc_utils import *

def get_bc1(text):
    return text[-8:]

def get_bc2(text):
    return text[8:16]

def get_bc3(text):
    return text[:8]

# TODO hopefully won't need this at some pt
def get_genotype_path_dict():
    d = {"129S1J_fasta/ncbi_dataset/data/GCA_029255695.1/GCA_029255695.1_ASM2925569v1_genomic.fna": "129S1J",
        "CASTJ_fasta/ncbi_dataset/data/GCA_029237265.1/GCA_029237265.1_ASM2923726v1_genomic.fna": "CASTJ",
        "AJ_fasta/ncbi_dataset/data/GCA_029255665.1/GCA_029255665.1_ASM2925566v1_genomic.fna": "AJ",
        "NODJ_fasta/ncbi_dataset/data/GCA_029234005.1/GCA_029234005.1_ASM2923400v1_genomic.fna": "NODJ",
        "NZOJ_fasta/ncbi_dataset/data/GCA_029233705.1/GCA_029233705.1_ASM2923370v1_genomic.fna": "NZOJ",
        "PWKJ_fasta/ncbi_dataset/data/GCA_029233695.1/GCA_029233695.1_ASM2923369v1_genomic.fna": "PWKJ",
        "WSBJ_fasta/ncbi_dataset/data/GCA_029233295.1/GCA_029233295.1_ASM2923329v1_genomic.fna": "WSBJ",
        "GRCm39.primary_assembly.genome.fa.gz": "B6J"}
    return d

def get_genotype_counts(files, ofile):
    """
    Given a list of klue anndata objects, merge
    the genotype counts together for each cell and output
    as a tsv
    """
    for i, f in enumerate(files):
        if "WSBJ_CASTJ" not in f:
            adata = sc.read(f)
            if i == 0:
                df = adata.to_df()
            else:
                df = df.merge(adata.to_df(),
                              how='outer',
                              left_index=True,
                              right_index=True)
    df.fillna(0, inplace=True)
    df.to_csv(ofile, index=True, sep='\t')

def rename_klue_genotype_cols(adata):
    """
    """
    d = get_genotype_path_dict()
    adata.var['genotype'] = adata.var.index.map(d)
    adata.var.set_index('genotype', inplace=True)
    adata.var.drop(adata.var.columns, axis=1, inplace=True)
    return adata

def add_meta_filter(mtx,
                    cgb,
                    cggn,
                    cgg,
                    wc,
                    bc_df,
                    kit,
                    chemistry,
                    sample_df,
                    min_counts,
                    ofile):
    
    """
    Add gene names to the kallisto output.
    Very initial filtering on min_counts.
    """
    
    data = sc.read_mtx(mtx)
    obs = pd.read_csv(cgb, header = None, sep="\t")  
    obs.columns = ["bc"]
    var = pd.read_csv(cgg, header = None, sep="\t")
    var.columns = ["gene_name"]
    genes = pd.read_csv(cggn, header = None, sep="\t")
    genes.columns = ["gene_id"]
    var['gene_id'] = genes['gene_id']
    
    X = data.X
    adata = anndata.AnnData(X=X, obs=obs, var=var)  
    print(adata.var)
    adata.obs.index = obs["bc"]

    adata.obs['bc1_sequence'] = adata.obs['bc'].apply(get_bc1)
    adata.obs['bc2_sequence'] = adata.obs['bc'].apply(get_bc2)
    adata.obs['bc3_sequence'] = adata.obs['bc'].apply(get_bc3)

    adata.obs['subpool'] = wc.subpool

    # merge in bc metadata
    temp = bc_df.copy(deep=True)
    temp = temp[['bc1_dt', 'well']].rename({'bc1_dt': 'bc1_sequence',
                                             'well': 'bc1_well'}, axis=1)
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_sequence')

    # merge in other bc information
    for bc in [2,3]:
        bc_name = f'bc{bc}'
        seq_col = f'bc{bc}_sequence'
        well_col = f'bc{bc}_well'
        bc_df = get_bcs(bc, kit, chemistry)
        bc_df.rename({'well': well_col,
                      bc_name: seq_col}, axis=1, inplace=True)
        adata.obs = adata.obs.merge(bc_df, how='left', on=seq_col)

    # merge in w/ sample-level metadata
    temp = sample_df.copy(deep=True)
    temp = temp.loc[(sample_df.plate==wc.plate)]
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_well')
    adata.var.set_index('gene_id', inplace=True)
    print(adata.var)

    # make all object columns string columns
    for c in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[c].dtype):
            adata.obs[c] = adata.obs[c].fillna('NA')

    # create new index for each cell
    adata.obs['cellID'] = adata.obs['bc1_well']+'_'+\
        adata.obs['bc2_well']+'_'+\
        adata.obs['bc3_well']+'_'+\
        adata.obs['subpool']+'_'+\
        adata.obs['plate']
    adata.obs.reset_index(drop=True)

    # make sure these are unique + set as index
    assert len(adata.obs.index) == len(adata.obs.cellID.unique().tolist())
    adata.obs.set_index('cellID', inplace=True)

<<<<<<< HEAD
    # remove non-multiplexed cells if from klue
    if klue:
        inds = adata.obs.loc[adata.obs.well_type=='Multiplexed'].index
        adata = adata[inds, :].copy()
        adata = rename_klue_genotype_cols(adata)

    if not klue:
        # filter based on min_counts in Snakefile
        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
        adata_filt = adata[adata.obs.n_counts >= min_counts,:]
        #adata.X = adata.layers['unspliced']
=======
    # filter based on min_counts in Snakefile         
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata = adata[adata.obs.n_counts >= min_counts,:]
    
    # filter out sample swaps with wrong multiplexed genotype 
    inds = adata.obs.loc[adata.obs['Genotype'] != 'WSBJ/CASTJ'].index
    adata = adata[inds, :].copy()
>>>>>>> b1a9c80c5c38be23112184abd83b6cb5a1b67e93

    adata.write(ofile)
    
def add_meta_klue(adata,
                         wc,
                         bc_df,
                         kit,
                         chemistry,
                         sample_df,
                         ofile):
    
    """
    Merge metadata in with klue output
    """
    
    adata = sc.read(adata)
    print(adata.obs.head())
    adata.obs.reset_index(inplace=True)
    adata.obs.columns = ['bc']

    adata.obs['bc1_sequence'] = adata.obs['bc'].apply(get_bc1)
    adata.obs['bc2_sequence'] = adata.obs['bc'].apply(get_bc2)
    adata.obs['bc3_sequence'] = adata.obs['bc'].apply(get_bc3)
    adata.obs['subpool'] = wc.subpool

    # merge in bc metadata
    temp = bc_df.copy(deep=True)
    temp = temp[['bc1_dt', 'well']].rename({'bc1_dt': 'bc1_sequence',
                                             'well': 'bc1_well'}, axis=1)
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_sequence')

    # merge in other bc information
    for bc in [2,3]:
        bc_name = f'bc{bc}'
        seq_col = f'bc{bc}_sequence'
        well_col = f'bc{bc}_well'
        bc_df = get_bcs(bc, kit, chemistry)
        bc_df.rename({'well': well_col,
                      bc_name: seq_col}, axis=1, inplace=True)
        adata.obs = adata.obs.merge(bc_df, how='left', on=seq_col)

    # merge in w/ sample-level metadata
    temp = sample_df.copy(deep=True)
    temp = temp.loc[(sample_df.plate==wc.plate)]
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_well')

    # make all object columns string columns
    for c in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[c].dtype):
            adata.obs[c] = adata.obs[c].fillna('NA')

    # create new index for each cell
    adata.obs['cellID'] = adata.obs['bc1_well']+'_'+\
        adata.obs['bc2_well']+'_'+\
        adata.obs['bc3_well']+'_'+\
        adata.obs['subpool']+'_'+\
        adata.obs['plate']
    adata.obs.reset_index(drop=True)

    # make sure these are unique + set as index
    assert len(adata.obs.index) == len(adata.obs.cellID.unique().tolist())
    adata.obs.set_index('cellID', inplace=True)
    
    # filter out sample swaps with wrong multiplexed genotype     
    inds = adata.obs.loc[(adata.obs.well_type=='Multiplexed') & (adata.obs['Genotype'] != 'WSBJ/CASTJ')].index
    adata = adata[inds, :].copy()
    adata = rename_klue_genotype_cols(adata)
        
    adata.write(ofile)


def make_subpool_sample_adata(infile, wc, ofile):
    adata = sc.read(infile)
    inds = []

    # filter on sample
    # TODO sample will be combo of {tissue}_{genotype}_{sex}_{rep}
    inds += adata.obs.loc[adata.obs['Mouse_Tissue_ID']==wc.sample].index.tolist()

    adata = adata[inds, :].copy()

    adata.write(ofile)

def run_scrublet(infile,
                 n_pcs,
                 min_counts,
                 min_cells,
                 min_gene_variability_pctl,
                 ofile):
    adata = sc.read(infile)

    # if number of cells is very low, don't call doublets, fill in
    if adata.X.shape[0] <= n_pcs:
        adata.obs['doublet_scores'] = 0

    # number of cells has to be more than number of PCs
    elif adata.X.shape[0] > 30:
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=min_counts,
                                                                  min_cells=min_cells,
                                                                  min_gene_variability_pctl=min_gene_variability_pctl,
                                                                  n_prin_comps=n_pcs)
        adata.obs['doublet_scores'] = doublet_scores

    adata.write(ofile)

def concat_adatas(adatas, ofile):
    for i, f in enumerate(adatas):
        if i == 0:
            adata = sc.read(f)

        else:
            temp = sc.read(f)

            temp.obs.reset_index(inplace=True)
            # if 'A10_A7_C3_Sublibrary_2_igvf_013' in temp.obs.cellID.tolist():
            #     import pdb; pdb.set_trace()
            temp.obs.set_index('cellID', inplace=True)
            adata = adata.concatenate(temp,
                        join='outer',
                        index_unique=None)
            adata.obs.reset_index(inplace=True)
            # if len(adata.obs.index) != len(adata.obs.cellID.unique().tolist()):
            #     import pdb; pdb.set_trace()
            adata.obs.set_index('cellID', inplace=True)

    adata.write(ofile)

def get_genotypes():
    g = ['WSBJ','NZOJ',
         'B6J','NODJ','129S1J',
         'CASTJ','AJ','PWKJ',
         'B6129S1F1J',
         'B6AF1J','B6PWKF1J',
         'B6NODF1J', 'B6WSBF1J',
         'B6CASTF1J', 'B6NZOF1J']
    return g

# def get_f1_genotype_pieces():
#     f1_genotypes = ['B6129S1F1J',
#         'B6AF1J','B6PWKF1J',
#         'B6NODF1J', 'B6WSBF1J',
#         'B6CASTF1J', 'B6NZOF1J']

def assign_demux_genotype(df):
    """
    Assigns a cell the genotype w/ the maximum of counts
    between the two genotypes that were loaded in the well
    the cell is from.

    Parameters:
        df (pandas DataFrame): DF of obs table for each cell w/
            klue counts for each genotype and multiplexed genotype
            columns
    """
    genotype_cols = get_genotypes()

    # restrict to nuclei w/ genetic multiplexing
    df = df.loc[df.well_type=='Multiplexed'].copy(deep=True)

    # fill nans once again
    df[genotype_cols] = df[genotype_cols].fillna(0)

    # loop through each multiplexed genotype combo
    # use those genotypes to determine which,
    # between the two, has the highest counts
    keep_cols = ['mult_genotype',
                 'mult_genotype_1',
                 'mult_genotype_2']+genotype_cols
    df = df[keep_cols]
    temp2 = pd.DataFrame()
    for g in df.mult_genotype.unique().tolist():
        # print(g)
        temp = df.loc[df.mult_genotype==g].copy(deep=True)

        g1 = temp.mult_genotype_1.unique().tolist()
        assert len(g1) == 1
        g1 = g1[0]

        g2 = temp.mult_genotype_2.unique().tolist()
        assert len(g2) == 1
        g2 = g2[0]

        # find the best match and report ties if the
        # values are the same
        temp['new_genotype'] = temp[[g1,g2]].idxmax(axis=1)
        temp.loc[temp[g1]==temp[g2], 'new_genotype'] = 'tie'

        temp2 = pd.concat([temp2, temp], axis=0)

    df = df.merge(temp2['new_genotype'], how='left',
                  left_index=True, right_index=True)
    df = df[['new_genotype']]

    assert len(df.loc[df.new_genotype.isnull()].index) == 0

    return df

def merge_kallisto_klue(f, genotypes, ofile):
    """
    Merge in the klue results with the kallisto results. Use
    a heuristic (which should be changeable / is subject to change)
    to determine the genotype assignment to each cell
    """

    adata = sc.read(f)
    df = pd.read_csv(genotypes, sep='\t')
    df.set_index('cellID', inplace=True)

    # make sure we won't dupe any cols
    assert len(set(df.columns.tolist())&set(adata.obs.columns.tolist())) == 0

    # merge in first; this way we have access to the genotypes that should
    # be in each well
    adata.obs = adata.obs.merge(df,
                                how='left',
                                left_index=True,
                                right_index=True)

    # assign genotype for multiplexed wells
    df = adata.obs.copy(deep=True)
    df = assign_demux_genotype(df)

    # merge in w/ adata and replace old values in "Genotype"
    # column for multiplexed wells with the klue results
    adata.obs = adata.obs.merge(df,
                                how='left',
                                left_index=True,
                                right_index=True)
    inds = adata.obs.loc[adata.obs.well_type=='Multiplexed'].index
    adata.obs.Genotype = adata.obs.Genotype.astype('str')
    adata.obs.loc[inds, 'Genotype'] = adata.obs.loc[inds, 'new_genotype']
    adata.obs.drop('new_genotype', axis=1, inplace=True)

    adata.write(ofile)
