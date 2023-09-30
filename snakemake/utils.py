import pandas as pd
import anndata
import scanpy as sc

def get_bc1(text):
    return text[-8:]

def get_bc2(text):
    return text[8:16]

def get_bc3(text):
    return text[:8]

def make_first_adata(adata, cgg, min_counts, ofile):
    """
    Add gene names to the kallisto output.
    Filter
    """
    adata = sc.read(adata)
    adata.obs.reset_index(inplace=True)
    adata.obs.columns = ['bc']

    adata.obs['bc1_sequence'] = adata.obs['bc'].apply(get_bc1)
    adata.obs['bc2_sequence'] = adata.obs['bc'].apply(get_bc2)
    adata.obs['bc3_sequence'] = adata.obs['bc'].apply(get_bc3)

    adata.var.reset_index(inplace=True)
    adata.var.columns = ['gene_name']
    names = pd.read_csv(cgg, sep="\t", header = None)
    names.columns = ['gene_id']
    adata.var['gene_id'] = names['gene_id']

    sc.pp.filter_cells(adata, min_counts = min_counts, inplace=True)
    adata.write(ofile)
