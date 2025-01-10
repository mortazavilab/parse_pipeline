import scanpy as sc
import pandas as pd
import leidenalg
import anndata as ad
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process tissue adata from pipeline with cellbender.')
    parser.add_argument('-f', '--file', type=str, help='input file path')

    args = parser.parse_args()

    if args.file:
        process_file(args.file)
    else:
        print("Please provide a file using -f option.")

def process_file(file_path):
    print("Processing file:", file_path)

    adata = sc.read_h5ad('../cellbender_tissues/' + file_path)
    filename = file_path
    tissue = filename.replace(".h5ad", "")

    print("Unfiltered:", adata.shape, flush=True)
    
    adata = adata[adata.obs['Genotype'] != 'tie']
    print("No ties:", adata.shape, flush=True)

    ### filtering data ###
    adata = adata[adata.obs['total_counts_raw'] > 500, :] # UMI cutoff
    print(">500 UMIs:", adata.shape, flush=True)
    
    adata = adata[adata.obs['total_counts_raw'] < 150000, :]
    print("<150000 UMIs:", adata.shape, flush=True)
    
    adata = adata[adata.obs['n_genes_by_counts_cb'] > 250, :] # min number of genes per cell
    print(">250 genes:", adata.shape, flush=True)

    if tissue == 'PBMC':
        adata = adata[adata.obs['pct_counts_mt_cb'] < 20, :]
        print("<20% mt:", adata.shape, flush=True)
    else:
        adata = adata[adata.obs['pct_counts_mt_cb'] < 1, :]
        print("<1% mt:", adata.shape, flush=True)

    # filter cells by doublet score (scrublet)
    adata = adata[adata.obs['doublet_score'] < 0.25, :]
    print("<0.25 doublet score:", adata.shape, flush=True)

    ### normalize the data ###
    sc.pp.normalize_total(adata, target_sum=1e4, layers=None, inplace=True) # Counts per 10k
    sc.pp.log1p(adata, layer=None)

    ### calculate highly variable genes ###
    # highly variable genes are used to compute the clustering 
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    adatas = adata[:, adata.var.highly_variable]

    print("Regressing....")
    sc.pp.regress_out(adatas, ['pct_counts_mt_cb','n_genes_by_counts_cb'])
    sc.pp.scale(adatas, max_value=10)

    print("PCA....")
    sc.tl.pca(adatas, svd_solver='arpack', mask_var="highly_variable")
    sc.pp.neighbors(adatas, n_neighbors=20, n_pcs=30)

    print("Clustering....")
    sc.tl.leiden(adatas, resolution = 1)
    sc.tl.umap(adatas, min_dist = 0.5, spread =2.0)

    adata.uns['neighbors'] = adatas.uns['neighbors']
    adata.uns['leiden'] = adatas.uns['leiden']
    adata.uns['umap'] = adatas.uns['umap']
    adata.obs['leiden'] = adatas.obs['leiden']
    adata.obsm = adatas.obsm
    adata.obsp = adatas.obsp

    adata.write_h5ad("../cellbender_tissues/" + tissue  +'_processed.h5ad')

if __name__ == "__main__":
    main()
