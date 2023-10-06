import pandas as pd
import anndata
import scanpy as sc
from bc_utils import *

def get_bc1(text):
    return text[-8:]

def get_bc2(text):
    return text[8:16]

def get_bc3(text):
    return text[:8]

def make_subpool_adata(adata,
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

    adata.obs['subpool'] = wc.subpool

    sc.pp.filter_cells(adata, min_counts=min_counts, inplace=True)

    # merge in bc metadata
    temp = bc_df.copy(deep=True)
    temp = temp[['bc1_dt', 'well']].rename({'bc1_dt': 'bc1_sequence',
                                             'well': 'bc1_well'}, axis=1)
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_sequence')

    # merge in other bc information
    for bc in [2,3]:
        seq_col = f'bc{bc}_sequence'
        well_col = f'bc{bc}_well'
        bc_df = get_bcs(bc, kit, chemistry)
        bc_df.rename({'well': well_col,
                      'sequence': seq_col}, axis=1, inplace=True)
        adata.obs = adata.obs.merge(bc_df, how='left', on=seq_col)

    # merge in w/ sample-level metadata
    temp = sample_df.copy(deep=True)
    temp = temp.loc[(sample_df.plate==wildcards.plate)]
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_well')
    adata.var.set_index('gene_id', inplace=True)

    # make all object columns string columns
    for c in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[c].dtype):
            adata.obs[c] = adata.obs[c].fillna('NA')

    adata.write(ofile)

def make_subpool_sample_adata(infile, wc, ofile):
    adata = sc.read(infile)
    inds = []

    # filter on sample
    # TODO sample will be combo of {tissue}_{genotype}_{sex}_{rep}
    inds += adata.obs.loc[adata.obs['Mouse_Tissue_ID']==wc.sample].index.tolist()

    adata = adata[inds, :].copy()

    adata.write(ofile)
