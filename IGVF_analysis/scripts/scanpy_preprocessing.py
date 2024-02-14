import scanpy as sc
import pandas as pd
import leidenalg
import anndata as ad
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('-f', '--file', type=str, help='input file path')

    args = parser.parse_args()

    if args.file:
        process_file(args.file)
    else:
        print("Please provide a file using -f option.")

def process_file(file_path):
    print("Processing file:", file_path)

    adata = sc.read_h5ad('preprocessed_tissues/' + file_path)
    adata.layers['raw_counts'] = adata.X.copy()
    adata.var_names_make_unique()


### filtering data ###

    sc.pp.filter_cells(adata, min_counts = 500) # UMI cutoff
    sc.pp.filter_cells(adata, max_counts = 150000)
    sc.pp.filter_cells(adata, min_genes = 250) # min number of genes per cell

# filtering cells by mt content 
    adata.var['mt'] = adata.var['gene_name'].str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 1, :]

# filter cells by doublet score (scrublet)
    adata = adata[adata.obs.doublet_scores < 0.25, :]

### normalize the data ###
    sc.pp.normalize_total(adata, target_sum=1e4) # Counts per 10k
#adata.layers['norm_10k'] = adata.X
    sc.pp.log1p(adata)
#adata.layers['norm_log1p'] = adata.X

### calculate highly variable genes ###
# highly variable genes are used to compute the clustering 
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adatas = adata[:, adata.var.highly_variable]



#saving 'raw' (counts per 10k and log normalized)
#adata.raw = adata
#adata.write_h5ad('Hypothalamus_pituitarty_processed_newt.h5ad')

    sc.pp.regress_out(adatas, ['pct_counts_mt','n_genes_by_counts'])
    sc.pp.scale(adatas, max_value=10)
#adata.layers['scale10'] = adata.X

    sc.tl.pca(adatas, svd_solver='arpack', use_highly_variable = True)
    sc.pp.neighbors(adatas, n_neighbors=20, n_pcs=30)

    sc.tl.leiden(adatas, resolution = 1)


    sc.tl.umap(adatas, min_dist = 0.5, spread =2.0)

    adata.uns['neighbors'] = adatas.uns['neighbors']
    adata.uns['leiden'] = adatas.uns['leiden']
    adata.uns['umap'] = adatas.uns['umap']
    adata.obs['leiden'] = adatas.obs['leiden']
    adata.obsm = adatas.obsm
#adata.varm = adatas.varm
    adata.obsp = adatas.obsp

#    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    filename = file_path
    filename2 = filename.replace("preprocessed.h5ad", "")

    adata.write_h5ad(filename2 +'_processed.h5ad')

if __name__ == "__main__":
    main()
