import scanpy as sc
import pandas as pd
import anndata as ad


tissues = ['Kidney']

for tis in tissues:

    paths = pd.read_table(f'../tissue_paths_cellbender/{tis}_paths.txt', header = None)

    ad_list = []

    # Loop through the list of files
    for file_name in paths[0]:
        print(file_name)
        adata = sc.read_h5ad(file_name + '/adata.h5ad')
        ad_list.append(adata)

    # # Concatenate the list of Anndata objects
    concatenated_ad = ad.concat(ad_list, merge = 'same', label = None)
    
    concatenated_ad.write_h5ad(f'../cellbender_tissues/{tis}.h5ad')
