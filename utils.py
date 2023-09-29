import pandas as pd
import anndata
import scanpy as sc

def get_bc1(text):
    return text[-8:]

def get_bc2(text):
    return text[8:16]

def get_bc3(text):
    return text[:8]

def make_first_adata(cgb, cggn, cgg, cgn, min_counts, ofile):
    obs = pd.read_csv(cgb, sep="\t", header = None)
    obs.columns = ['bc']

    obs['bc1_sequence'] = obs['bc'].apply(get_bc1)
    obs['bc2_sequence'] = obs['bc'].apply(get_bc2)
    obs['bc3_sequence'] = obs['bc'].apply(get_bc3)

    var = pd.read_csv(cgg, sep="\t", header = None)
    var.columns = ['gene_id']
    names = pd.read_csv(cggn, sep="\t", header = None)
    names.columns = ['gene_name']
    var['gene_name'] = names['gene_name']

    adata = sc.read_mtx(cgn)
    X = adata.X
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    sc.pp.filter_cells(adata, min_counts = min_counts, inplace=True) # filter by 500 UMIs
    adata.write(ofile)
