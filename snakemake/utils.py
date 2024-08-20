import pandas as pd
import numpy as np
from snakemake.io import expand
import anndata
import scanpy as sc
import scrublet as scr
import shutil
import os
from cellbender.remove_background.downstream import anndata_from_h5

############################################################################################################
############################################### Barcode stuff ##############################################
############################################################################################################
    
# from the parse pipeline
def load_bc_dict(fname, verb=False):
    """ Load barcode edit dict
    """
    with open(fname, 'r') as INFILE:
        bc_dict = json.load(INFILE)
        if verb:
            print(f"Loaded {fname}")

    # Top level has int keys and holds default dicts
    new_dict = {}
    for k, v in bc_dict.items():
        new_dict[int(k)] = defaultdict(list, bc_dict[k])

    return new_dict

def get_bc_round_set(kit, chemistry):
    KIT_INT_DICT = {'custom_1': 1, 'WT': 48, 'WT_mini': 12, 'WT_mega': 96}
    kit_n = KIT_INT_DICT[kit]
    if kit_n == 12:
        bc_round_set = [['bc1','n24_v4'], ['bc2','v1'], ['bc3','v1']]
    if kit_n == 96:
        bc_round_set = [['bc1','n192_v4'], ['bc2','v1'], ['bc3','v1']]
    if kit_n == 48:
        bc_round_set = [['bc1','v2'], ['bc2','v1'], ['bc3','v1']]

    if kit == 'WT' and chemistry == 'v2':
        bc_round_set = [['bc1', 'n96_v4'], ['bc2', 'v1'], ['bc3', 'v1']]

    return bc_round_set

def load_barcodes(kit, chemistry):
    """
    Load the barcodes. Adapted from the Parse biosciences pipeline.

    Returns:
        edit_dict_set (dict): Dict for barcode<1,2,3> with
            key: query bc
            item: corrected bc
    """
    pkg_path = os.path.dirname(__file__)
    bc_path = '/'.join(pkg_path.split('/')[:-1])+'/barcodes/'

    bc_round_set = get_bc_round_set(kit, chemistry)

    edit_dict_set = {}
    for entry in bc_round_set:
        bc = entry[0]
        ver = entry[1]
        fname = bc_path + 'bc_dict_{}.json'.format(ver)
        edit_dict = load_bc_dict(fname)
        edit_dict_set[bc] = edit_dict

    return edit_dict_set

# From the Parse biosciences pipeline
def load_barcodes_set(kit, chemistry):
    """
    Load the barcodes. Adapted from the Parse biosciences pipeline.
    """
    pkg_path = os.path.dirname(__file__)
    bc_path = '/'.join(pkg_path.split('/')[:-1])+'/barcodes/'

    bc_round_set = get_bc_round_set(kit, chemistry)

    bc_set = {}
    for entry in bc_round_set:
        bc = entry[0]
        ver = entry[1]
        fname = bc_path + '/bc_data_{}.csv'.format(ver)
        bc_df = pd.read_csv(fname)
        bcs = set(bc_df.sequence.tolist())
        bc_set[bc] = bcs

    return bc_set

def get_bc1_matches(kit, chemistry):
    pkg_path = os.path.dirname(os.path.dirname(__file__))
    # pkg_path = '/'.join(pkg_path.split('/')[:-1])
    bc_round_set = get_bc_round_set(kit, chemistry)

    # determine file to use
    for entry in bc_round_set:
        if entry[0] == 'bc1':
            ver = entry[1]

    # read in and restructure such that each dt bc is
    # matched with its randhex partner from the same well
    fname = pkg_path+'/barcodes/bc_data_{}.csv'.format(ver)
    df = pd.read_csv(fname)
    df.loc[df.well.duplicated(keep=False)].sort_values(by='well')
    drop_cols = ['bci', 'uid', 'type']
    bc1_dt = df.loc[df['type'] == 'T'].drop(drop_cols, axis=1)
    bc1_dt.rename({'sequence': 'bc1_dt'}, axis=1, inplace=True)
    bc1_randhex = df.loc[df['type'] == 'R'].drop(drop_cols, axis=1)
    bc1_randhex.rename({'sequence': 'bc1_randhex'}, axis=1, inplace=True)
    bc_df = bc1_dt.merge(bc1_randhex, on='well')

    return bc_df

def get_bcs(bc, kit, chemistry):
    """
    Parameters:
        bc (int): {1,2,3}
    """
    pkg_path = os.path.dirname(os.path.dirname(__file__))
    bc_round_set = get_bc_round_set(kit, chemistry)

    bc_name = f'bc{bc}'

    # determine file to use
    for entry in bc_round_set:
        if entry[0] == bc_name:
            ver = entry[1]

    fname = pkg_path+'/barcodes/bc_data_{}.csv'.format(ver)
    df = pd.read_csv(fname)
    df.loc[df.well.duplicated(keep=False)].sort_values(by='well')

    if bc == 2 or bc == 3:
        assert len(df['type'].unique().tolist()) == 1

    drop_cols = ['bci', 'uid', 'type']
    df.drop(drop_cols, axis=1, inplace=True)
    df.rename({'sequence': bc_name}, axis=1, inplace=True)

    return df

############################################################################################################
############################################## Snakemake stuff #############################################
############################################################################################################
 
def parse_config(fname):
    df = pd.read_csv(fname, sep='\t')
    df['path'] = df.fastq.str.rsplit('/', n=2, expand=True)[0]+'/'
    df['path2'] = df.fastq.str.rsplit('/', n=1, expand=True)[0]+'/'
    df['fastqs'] = df.apply(lambda x: [x.fastq, x.fastq_r2] \
                                    if not pd.isna(x.fastq_r2) else [x.fastq],
                                    axis=1)
    df['fastq_pairs'] = df.fastqs

    return df

def get_df_info(wc, df, col):
    temp = df.copy(deep=True)
    temp = temp.loc[(temp.plate==wc.plate)&\
                    (temp.subpool==wc.subpool)&\
                    (temp.lane==wc.lane)&\
                    (temp.run==int(wc.run))]
    assert len(temp.index) == 1
    return temp[col].values[0]

def parse_sample_df(fname):
    df = pd.read_csv(fname)
    df.rename({'Experiment': 'plate'}, axis=1, inplace=True)

    # add multiplexed genotypes if relevant
    g_cols = ['mult_genotype_1', 'mult_genotype_2']
    df[g_cols] = df.Genotype.str.split('/', expand=True)

    # adjust single-genotype wells
    df.loc[df.well_type=='Single', g_cols] = np.nan

    # add a multiplexed genotype column
    inds = df.loc[df.well_type=='Multiplexed'].index
    df['mult_genotype'] = np.nan
    df.loc[inds, 'mult_genotype'] = df.loc[inds, g_cols[0]].astype(str)+'_'+\
                                   df.loc[inds, g_cols[1]].astype(str)

    # checks
    for g in g_cols:
        assert len(df.loc[(df.well_type=='Single')&(df[g].notnull())].index) == 0

    return df
    
    
############################################################################################################
#################################### Add metadata to cellbender output #####################################
############################################################################################################

def get_bc1(text):
    return text[-8:]

def get_bc2(text):
    return text[8:16]

def get_bc3(text):
    return text[:8]

def add_meta_filter(filt_h5,
                    unfilt_adata,
                    wc,
                    bc_df,
                    kit,
                    chemistry,
                    sample_df,
                    ofile):

    """
    Format cellbender output adata.
    """

    adata = anndata_from_h5(filt_h5)
   
    adata.obs['bc'] = adata.obs.index
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
    adata.var_names = adata.var['gene_id']
    
    # get gene names back....
    adata_kallisto = sc.read_h5ad(unfilt_adata)
    adata_kallisto.var_names = adata_kallisto.var['gene_id']
    adata.var['gene_name'] = adata_kallisto.var['gene_name']

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
    adata = adata[~adata.obs['Genotype'].isin(['WSBJ/CASTJ', 'AJ/129S1J', 'PWKJ/CASTJ'])].copy()
    
    # add raw counts layer
    adata.layers['raw_counts'] = adata.X.copy()
    
        
    adata.write_h5ad(ofile)
    


    
############################################################################################################
################################## Merge klue results, merge final adata ###################################
############################################################################################################

    
def touch_dummy(ofile):
    """
    Touch a dummy output file
    """
    open(ofile, 'a').close()

def is_dummy(fname):
    """
    Check if a file is a dummy file (ie whether is empty)
    Returns True if it's a dummy, False if not
    """
    return os.stat(fname).st_size == 0

# TODO maybe splitting these up will make things easier in the future
def get_founder_genotypes():
    g = ['WSBJ','NZOJ',
         'B6J','NODJ','129S1J',
         'CASTJ','AJ','PWKJ']
    return g

def get_f1_genotypes():
    g = ['B6129S1F1J',
    'B6AF1J','B6PWKF1J',
    'B6NODF1J', 'B6WSBF1J',
    'B6CASTF1J', 'B6NZOF1J']
    return g

def get_genotypes():
    g = get_founder_genotypes()
    g += get_f1_genotypes()

    return g

def get_f1_founder_genotype_dict():
    d = {'B6129S1F1J':'129S1J',
    'B6AF1J':'AJ',
    'B6PWKF1J':'PWKJ',
    'B6NODF1J':'NODJ',
    'B6WSBF1J':'WSBJ',
    'B6CASTF1J':'CASTJ',
    'B6NZOF1J':'NZOJ'}

    # make sure all keys are f1 genotypes
    # and all items (sans B6J) are founder genotypes
    assert set(list(d.keys())) == set(get_f1_genotypes())
    assert set(list(d.values())) == set(get_founder_genotypes())-set(['B6J'])

    return d

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
    print(genotype_cols)
    
    genotype_cols = [c for c in genotype_cols if c in df.columns.tolist()]

    # restrict to nuclei w/ genetic multiplexing
    df = df.loc[df.well_type=='Multiplexed'].copy(deep=True)

    # fill nans once again
    df[genotype_cols] = df[genotype_cols].fillna(0)

    # loop through each multiplexed genotype combo
    # use those genotypes to determine which, between the two, has the highest counts
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

        # find the best match and report ties if the values are the same
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
    # if we didn't run klue
    if is_dummy(genotypes):
        shutil.copy(f, ofile)
    # otherwise, decide which genotype each cell is
    else:
        adata = sc.read_h5ad(f)
        print("HELLO?")
        
        df = pd.read_csv(genotypes, sep='\t')
        df.set_index('cellID', inplace=True)

        # make sure we won't dupe any cols
        assert len(set(df.columns.tolist())&set(adata.obs.columns.tolist())) == 0

        # merge in first; this way we have access to the genotypes that should be in each well
        adata.obs = adata.obs.merge(df,how='left',left_index=True,right_index=True)

        # assign genotype for multiplexed wells
        df = adata.obs.copy(deep=True)
        df = assign_demux_genotype(df)

        # merge in w/ adata and replace old values in "Genotype" column for multiplexed wells with the klue results
        adata.obs = adata.obs.merge(df,how='left',left_index=True,right_index=True)
        
        inds = adata.obs.loc[adata.obs.well_type=='Multiplexed'].index
        adata.obs.Genotype = adata.obs.Genotype.astype('str')
        adata.obs.loc[inds, 'Genotype'] = adata.obs.loc[inds, 'new_genotype']
        adata.obs.drop('new_genotype', axis=1, inplace=True)  
        
        adata = adata[adata.obs['Genotype'] != 'tie'].copy() # aaaahh remove TIE's
        
        print(adata.obs['Mouse_Tissue_ID'].value_counts())
        print(adata.obs['Original_Mouse_Tissue_ID'].value_counts())
        
        # adjust mouse_tissue_id
        ms1 = ['B6J','AJ','WSBJ','129S1J']
        ms2 = ['NODJ','PWKJ','NZOJ','CASTJ']

        # Define a function to update 'Mouse_Tissue_ID' based on conditions
        def update_mouse_tissue_id(row):
            if row['well_type'] == "Multiplexed":
                if row['Genotype'] in ms1:
                    return row['Multiplexed_sample1']
                elif row['Genotype'] in ms2:
                    return row['Multiplexed_sample2']
            return row['Mouse_Tissue_ID']

        adata.obs['Mouse_Tissue_ID'] = adata.obs['Mouse_Tissue_ID'].astype(str)
        meta = adata.obs

        # Apply the function to update the 'Mouse_Tissue_ID' column
        meta['Mouse_Tissue_ID'] = meta.apply(update_mouse_tissue_id, axis=1)
        
        adata.obs['Original_Mouse_Tissue_ID'] = adata.obs['Mouse_Tissue_ID'] 
        

        
        adata.obs['Mouse_Tissue_ID'] = meta['Mouse_Tissue_ID']
        
        adata.write_h5ad(ofile)
        
############################################################################################################
################################################# Scrublet #################################################
############################################################################################################

def make_subpool_sample_adata(infile, wc, ofile):
    adata = sc.read_h5ad(infile)
    inds = []
    inds += adata.obs.loc[adata.obs['Original_Mouse_Tissue_ID']==wc.sample].index.tolist()

    adata = adata[inds, :].copy()

    adata.write_h5ad(ofile)
    
def run_scrublet(infile,
                 n_pcs,
                 min_cells,
                 min_gene_variability_pctl,
                 ofile):
    adata = sc.read_h5ad(infile)

    # if number of cells is very low, don't call doublets, fill in
    if adata.X.shape[0] < n_pcs*20:
        adata.obs['doublet_scores'] = 0

    # number of cells has to be more than number of PCs
    elif adata.X.shape[0] >= n_pcs*20:
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_cells=min_cells,
                                                                  min_gene_variability_pctl=min_gene_variability_pctl,
                                                                  n_prin_comps=n_pcs)
        adata.obs['doublet_scores'] = doublet_scores

    adata.write_h5ad(ofile)
    
############################################################################################################
############################################## Final concat ################################################
############################################################################################################

def concat_adatas(adatas, ofile):
    var_dfs = []
    print(adatas)
    
    for i, f in enumerate(adatas):
        temp = sc.read_h5ad(f)
        
        temp.var.drop(columns=['feature_type', 'genome', 'ambient_expression', 'cellbender_analyzed'], inplace=True)
        
#         plate = temp.obs['plate'][0]
#         sample = temp.obs['Mouse_Tissue_ID'][0]
#         subpool = temp.obs['subpool'][0]
        
#         temp.var.rename(columns={
#             'ambient_expression': f'ambient_expression_{sample}_{subpool}_{plate}',
#             'cellbender_analyzed': f'cellbender_analyzed_{sample}_{subpool}_{plate}'
#         }, inplace=True)
        
        var_dfs.append(temp.var)
        
        if i == 0:
            # Initialize the combined adata object
            adata = temp
        else:
            # Concatenate obs and X
            temp.obs.reset_index(inplace=True)
            temp.obs.set_index('cellID', inplace=True)
            adata = anndata.concat([adata, temp],
                                   join='outer',
                                   index_unique=None)
            adata.obs.reset_index(inplace=True)
            adata.obs.set_index('cellID', inplace=True)
    
    # Concatenate all var DataFrames and drop duplicates
    combined_var = pd.concat(var_dfs, axis=1, join='outer')
    combined_var = combined_var.loc[:, ~combined_var.columns.duplicated()]
    adata.var = combined_var
    
    adata.write(ofile)

