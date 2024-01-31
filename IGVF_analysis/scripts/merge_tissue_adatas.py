import scanpy as sc
import pandas as pd
import anndata as ad


tissues = ['Adrenal','CortexHippocampus','Gastrocnemius','GonadsFemale','GonadsMale','Heart','HypothalamusPituitary','Kidney','Liver']

for tis in tissues:

    paths = pd.read_table('/share/crsp/lab/seyedam/weberrl/IGVF_analysis/tissue_paths/' + tis + '_paths.txt', header = None)

    ad_list = []

    # Loop through the list of files
    for file_name in paths[0]:
        print(file_name)
        adata = sc.read_h5ad(file_name + '/adata.h5ad')
        ad_list.append(adata)

    # # Concatenate the list of Anndata objects
    concatenated_ad = ad.concat(ad_list, merge = 'same')
    
    concatenated_ad.write_h5ad('/share/crsp/lab/seyedam/weberrl/IGVF_analysis/preprocessed_tissues/' + tis + 'preprocessed.h5ad')
