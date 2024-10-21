import pandas as pd
import numpy as np
from snakemake.io import expand
import anndata
import scanpy as sc
import shutil
import os
import json
from cellbender.remove_background.downstream import anndata_from_h5

############################################################################################################
############################################## Barcode parsing #############################################
############################################################################################################
def get_bc1(text):
    return text[-8:]

def get_bc2(text):
    return text[8:16]

def get_bc3(text):
    return text[:8]

def get_bc_round_set(kit, chemistry):
    valid_kits = {'custom_1', 'WT', 'WT_mini', 'WT_mega'}
    valid_chemistries = {'v1', 'v2', 'v3'}
    
    # Check if the provided kit and chemistry are valid
    if kit not in valid_kits:
        raise ValueError(f"Invalid kit: {kit}. Valid options are: {', '.join(valid_kits)}.")
    if chemistry not in valid_chemistries:
        raise ValueError(f"Invalid chemistry: {chemistry}. Valid options are: {', '.join(valid_chemistries)}.")

        
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
        
    # support for v3 chemistry
    if kit == 'WT_mini' and chemistry == 'v3':
        bc_round_set = [['bc1', 'n26_R1_v3_3'], ['bc2', 'v1'], ['bc3', 'R3_v3']]
        
    if kit == 'WT' and chemistry == 'v3':
        bc_round_set = [['bc1', 'n102_R1_v3_3'], ['bc2', 'v1'], ['bc3', 'R3_v3']]
        
    if kit == 'WT_mega' and chemistry == 'v3':
        bc_round_set = [['bc1', 'n208_R1_v3_3'], ['bc2', 'v1'], ['bc3', 'R3_v3']]

    return bc_round_set


def get_bc1_matches(kit, chemistry):
    pkg_path = os.path.dirname(os.path.dirname(__file__))
    bc_round_set = get_bc_round_set(kit, chemistry)

    # determine file to use
    for entry in bc_round_set:
        if entry[0] == 'bc1':
            ver = entry[1]

    # read in and restructure such that each dt bc is
    # matched with its randhex partner from the same well
    fname = pkg_path+'/barcodes/bc_data_{}.csv'.format(ver)
    df = pd.read_csv(fname)
    
    if chemistry == "v3":
        df = df[df['stype'] != "X"]
        df.rename({'stype': 'type'}, axis=1, inplace=True)
    
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
    
    if chemistry == "v3":
        if 'stype' in df.columns:
            df = df[df['stype'] != "X"]
            df.rename({'stype': 'type'}, axis=1, inplace=True)
        else:
            df = df[df['type'] != "X"]
        
    df.loc[df.well.duplicated(keep=False)].sort_values(by='well')

    if bc == 2 or bc == 3:
        assert len(df['type'].unique().tolist()) == 1

    df = df[['sequence', 'well']]
    df.rename({'sequence': bc_name}, axis=1, inplace=True)

    return df


def create_r1_RT_replace(kit, chemistry, ofile):
    # Use the get_bc_round_set function to get the appropriate bc_round_set
    bc_round_set = get_bc_round_set(kit, chemistry)
    bc, file_suffix = bc_round_set[0]  # Only use the first entry, which corresponds to bc1
    
    pkg_path = os.path.dirname(__file__)
    bc_path = '/'.join(pkg_path.split('/')[:-1])+'/barcodes/'
    filename = bc_path + f"bc_data_{file_suffix}.csv"
    
    df = pd.read_csv(filename)
    
    if chemistry == "v3":
        df = df[df['stype'] != "X"]
        df.rename({'stype': 'type'}, axis=1, inplace=True)
    
    # Filter rows based on 'type' column
    df_R = df[df['type'] == 'R']['sequence']
    df_T = df[df['type'] == 'T']['sequence']
    
    # Format
    formatted_sequences = [f"{r} *{t}" for r, t in zip(df_R, df_T)]
    
    output_dir = os.path.dirname(ofile)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Write
    output_filename = ofile
    with open(output_filename, 'w') as f:
        f.write("\n".join(formatted_sequences))


def create_r1r2r3_file(kit, chemistry, ofile):
    def generate_wells():
        return [f"{row}{col}" for row in "ABCDEFGH" for col in range(1, 13)]
    
    # Get barcode files based on kit and chemistry
    bc_round_set = get_bc_round_set(kit, chemistry)
    pkg_path = os.path.dirname(__file__)
    bc_path = os.path.join(pkg_path, '..', 'barcodes')
    bc_files = {bc[0]: os.path.join(bc_path, f"bc_data_{bc[1]}.csv") for bc in bc_round_set}
        
    # Read barcode data
    bc_data = {}
    for bc in bc_files:
        df = pd.read_csv(bc_files[bc])
        
        if chemistry == "v3":
            df.rename({'stype': 'type'}, axis=1, inplace=True)
            
        bc_data[bc] = df
            
    # Generate all well identifiers
    all_wells = generate_wells()

    # Create empty dfs for barcode types L/T and type R
    df_l_t = pd.DataFrame({'well': all_wells})
    df_r = pd.DataFrame({'well': all_wells})
    
    # Filter and merge DataFrames for type L/T
    for bc in ['bc1', 'bc2', 'bc3']:
        df = bc_data[bc]
        df_l_t_filtered = df[df['type'].isin(['L', 'T'])][['well', 'sequence']]
        df_l_t_filtered = df_l_t_filtered.rename(columns={'sequence': f'{bc}_sequence'})
        df_l_t = df_l_t.merge(df_l_t_filtered, on='well', how='left')
        
    # Format  
    df_l_t.fillna('-', inplace=True)
    
    # Filter and merge dfs for type R
    for bc in ['bc1', 'bc2', 'bc3']:
        df = bc_data[bc]
        df_r_filtered = df[df['type'] == 'R'][['well', 'sequence']]
        df_r_filtered = df_r_filtered.rename(columns={'sequence': f'{bc}_sequence'})
        df_r = df_r.merge(df_r_filtered, on='well', how='left')
        
    # Format
    df_r.fillna('-', inplace=True)
    
    # combine
    combined_df = pd.concat([df_l_t, df_r], ignore_index=True)
    
    sequence_columns = ['bc3_sequence','bc2_sequence','bc1_sequence']
    sequences_df = combined_df[sequence_columns]
    sequences_df = sequences_df[~(sequences_df[sequence_columns] == '-').all(axis=1)]

    sequences_df.to_csv(ofile, sep=' ', index=False, header=False)


############################################################################################################
############################################## Snakemake stuff #############################################
############################################################################################################
def get_fa_link(wc, config):
    genotype = wc.genotype
    
    if genotype in get_f1_genotypes():
        d = get_f1_founder_genotype_dict()
        genotype = d[genotype]
        
    link = config['ref']['genome']['link'][genotype]
    
    return link

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
    df['mult_genotype'] = '' # maybe change to this but needs to be tested
    df.loc[inds, 'mult_genotype'] = df.loc[inds, g_cols[0]].astype(str)+'_'+\
                                   df.loc[inds, g_cols[1]].astype(str)

    # checks
    for g in g_cols:
        assert len(df.loc[(df.well_type=='Single')&(df[g].notnull())].index) == 0

    return df

def get_subpool_fastqs(wc, df, config, how, read=None):
    """
    Get list of fastqs from the same subpool. Can
    either return as a Python list of strings or a
    formatted string list read to pass to a shell cmd.

    Parameters:
        how (str): {'str', 'list'}. 'list' will return
            Python list of str var. 'str' will return
            Python string
    """
    temp = df.copy(deep=True)
    temp = temp.loc[(temp.plate==wc.plate)&\
                    (temp.subpool==wc.subpool)]

    if how == 'list':
        reads = [read for i in range(len(temp.index))]
        if read == 'R1':
            entry = 'fastq_r1'
        elif read == 'R2':
            entry = 'fastq_r2'
        return expand(expand(config['raw'][entry],
                        zip,
                        lane=temp['lane'].tolist(),
                        run=temp['run'].tolist(),
                        allow_missing=True),
                        plate=wc.plate,
                        subpool=wc.subpool)

    elif how == 'str':
        r1s = expand(expand(config['raw']['fastq_r1'],
                        zip,
                        lane=temp['lane'].tolist(),
                        run=temp['run'].tolist(),
                        allow_missing=True),
                        plate=wc.plate,
                        subpool=wc.subpool)
        r2s = expand(expand(config['raw']['fastq_r2'],
                        zip,
                        lane=temp['lane'].tolist(),
                        run=temp['run'].tolist(),
                        allow_missing=True),
                        plate=wc.plate,
                        subpool=wc.subpool)
        fastq_str = ''
        for r1, r2 in zip(r1s, r2s):
            fastq_str+=f' {r1} {r2}'
        return fastq_str
    
    
def get_subpool_adatas(df, sample_df, wc, cfg_entry):
    """
    Get adatas that belong to the same subpool across the
    different pairs of genotypes ({mult_genotype_1},{mult_genotype_2})
    that klue was run on
    """
    # restrict to the plates in our input sample set
    temp = df.merge(sample_df, on='plate', how='inner')

    # restrict to this subpool / plate
    temp = temp.loc[temp.plate==wc.plate]
    temp = temp.loc[temp.subpool==wc.subpool]

    # restrict to multiplexed wells
    temp = temp.loc[temp.mult_genotype.notnull()]

    files = expand(cfg_entry,
                   zip,
                   plate=temp.plate.tolist(),
                   subpool=temp.subpool.tolist(),
                   mult_genotype_1=temp.mult_genotype_1.tolist(),
                   mult_genotype_2=temp.mult_genotype_2.tolist())
    files = list(set(files))

    return files

def get_subset_tissues(df, sample_df):
    temp = df.merge(sample_df, on='plate', how='inner')
    tissues = temp.Tissue.unique().tolist()
    return tissues
    
############################################################################################################
########################################## Format kallisto output ##########################################
############################################################################################################

def make_adata_from_kallisto(mtx,
                             cgb,
                             cggn,
                             cgg,
                             bc_df,
                             kit,
                             chemistry,
                             run_info,
                             bc_info,
                             config_tsv,
                             wc,
                             ofile,
                             knee_file,
                             report_df,
                             cb_df
                            ):

    """
    Make adata from unfiltered kallisto output using total mtx
    """

    print(ofile)
    data = sc.read_mtx(mtx)
    obs = pd.read_csv(cgb, header = None, sep="\t")
    obs.columns = ["bc"]
    var = pd.read_csv(cgg, header = None, sep="\t")
    var.columns = ["gene_id"]
    genes = pd.read_csv(cggn, header = None, sep="\t")
    genes.columns = ["gene_name"]
    var['gene_name'] = genes['gene_name']

    X = data.X
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.var_names = adata.var['gene_name']
    adata.var_names_make_unique()
    
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    
    adata.obs['bc1_sequence'] = adata.obs['bc'].apply(get_bc1)
    adata.obs['bc2_sequence'] = adata.obs['bc'].apply(get_bc2)
    adata.obs['bc3_sequence'] = adata.obs['bc'].apply(get_bc3)

    # merge in bc metadata
    temp = bc_df.copy(deep=True)
    temp = temp[['bc1_dt', 'well']].rename({'bc1_dt': 'bc1_sequence','well': 'bc1_well'}, axis=1)

    adata.obs = adata.obs.merge(temp, how='left', on='bc1_sequence')

    for bc in [2,3]:
        bc_name = f'bc{bc}'
        seq_col = f'bc{bc}_sequence'
        well_col = f'bc{bc}_well'
        bc_df = get_bcs(bc, kit, chemistry)
        bc_df.rename({'well': well_col,
                      bc_name: seq_col}, axis=1, inplace=True)
        adata.obs = adata.obs.merge(bc_df, how='left', on=seq_col)
        
    adata.obs.index = obs["bc"]

    ################## extra stuff for report ################## 
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    ### knee plot df
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    knee_df = pd.DataFrame({
        'Barcode rank': range(len(knee)),
        'UMIs': knee})
    
    knee_df.to_csv(knee_file, index=False)
    
    ### pseudoaligment stats
    json_file_path = run_info
    with open(json_file_path, 'r') as file:
        data_dict = json.load(file)

    df = pd.DataFrame([data_dict])
    df = df.drop(columns=['call', 'start_time', 'index_version', 'n_bootstraps', 'n_targets'])
    
    ### barcode stats
    json_file_path = bc_info
    with open(json_file_path, 'r') as file:
        data_dict = json.load(file)

    bc_info_df = pd.DataFrame([data_dict])
    
    
    #total_umis = np.sum(adata.X)
    df['Total UMIs'] = bc_info_df['numBarcodeUMIs']
    df['Sequencing saturation'] = 1 - (df['Total UMIs'] / df['n_pseudoaligned'])
    df['Sequencing saturation'] = round(df['Sequencing saturation'], 3)
    
    
    df['Total UMIs'] = format(int(df['Total UMIs']),',')
    
    df['n_processed'] = format(int(df['n_processed']),',')
    df['n_pseudoaligned'] = format(int(df['n_pseudoaligned']),',')
    df['n_unique'] = format(int(df['n_unique']),',')
    

    # Rename the remaining columns
    df = df.rename(columns={
        'n_processed': 'Total reads',
        'n_pseudoaligned': 'Total pseudoaligned',
        'n_unique': 'Total unique',
        'p_pseudoaligned': 'Pct. pseudoaligned',
        'p_unique': 'Pct. uniquely pseudoaligned',
        'kallisto_version': 'kallisto version'
    })
    
    df['Subpool'] = wc.subpool
    

    ### other QC stats
    adata_200 = adata[adata.obs['total_counts'] >= 200]
    df['Num. cells >= 200 UMI'] = f"{int(adata_200.shape[0]):,}"
    df['Median UMIs/cell (200 UMIs)'] = f"{int(adata_200.obs['total_counts'].median()):,}"
    df['Mean UMIs/cell (200 UMIs)'] = f"{round(adata_200.obs['total_counts'].mean(), 1):,.1f}"
    df['Median genes/cell (200 UMIs)'] = f"{int(adata_200.obs['n_genes_by_counts'].median()):,}"
    df['Mean genes/cell (200 UMIs)'] = f"{round(adata_200.obs['n_genes_by_counts'].mean(), 1):,.1f}"

    adata_500 = adata[adata.obs['total_counts'] >= 500]
    df['Num. cells >= 500 UMI'] = f"{int(adata_500.shape[0]):,}"
    df['Median UMIs/cell (500 UMIs)'] = f"{int(adata_500.obs['total_counts'].median()):,}"
    df['Mean UMIs/cell (500 UMIs)'] = f"{round(adata_500.obs['total_counts'].mean(), 1):,.1f}"
    df['Median genes/cell (500 UMIs)'] = f"{int(adata_500.obs['n_genes_by_counts'].median()):,}"
    df['Mean genes/cell (500 UMIs)'] = f"{round(adata_500.obs['n_genes_by_counts'].mean(), 1):,.1f}"
 
    df.to_csv(report_df, index=False)
    
    
    ### cb settings -- don't need to run CB first, get from config file
    cb_settings = pd.read_csv(config_tsv, sep='\t')
    
    columns_of_interest = ['plate', 'subpool', 'droplets_included', 'learning_rate', 'expected_cells']
    available_columns = [column for column in columns_of_interest if column in cb_settings.columns]
    if available_columns:
        cb_settings = cb_settings[available_columns]
        cb_settings = cb_settings.drop_duplicates().reset_index(drop=True)
        
    cb_settings.to_csv(cb_df, index=False)
    
    # WRITE
    adata.write(ofile)


############################################################################################################
################################################### Klue ###################################################
############################################################################################################

def get_genotype_counts(files, ofile):
    """
    Given a list of klue anndata objects, merge
    the genotype counts together for each cell and output
    as a tsv.
    """
    # If we have no genetic demultiplexing
    if not files:
        touch_dummy(ofile)
    # Otherwise, merge and output merged summary
    else:
        for i, f in enumerate(files):
            # Exclude specific genotype combinations
            if all(exclude not in f for exclude in ["WSBJ_CASTJ", "AJ_129S1J", "PWKJ_CASTJ"]):
                adata = sc.read_h5ad(f)
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
    adata.var['genotype'] = adata.var['gene_name'].str.split('/').str[-1].str.split('.').str[0]
    adata.var.set_index('genotype', inplace=True)
    adata.var.drop(adata.var.columns, axis=1, inplace=True)
    return adata

def add_meta_klue(mtx,
                    cgb,
                    cggn,
                    cgg,
                    wc,
                    bc_df,
                    kit,
                    chemistry,
                    sample_df,
                    ofile):

    """
    Add gene names to the klue output.
    """

    data = sc.read_mtx(mtx)
    obs = pd.read_csv(cgb, header = None, sep="\t")
    obs.columns = ["bc"]
    var = pd.read_csv(cgg, header = None, sep="\t")
    var.columns = ["gene_id"]
    genes = pd.read_csv(cggn, header = None, sep="\t")
    genes.columns = ["gene_name"]
    var['gene_name'] = genes['gene_name']

    X = data.X
    adata = anndata.AnnData(X=X, obs=obs, var=var)
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
    print(adata.obs)

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

    inds = adata.obs.loc[adata.obs.well_type=='Multiplexed'].index
    adata = adata[inds, :].copy()
    adata = rename_klue_genotype_cols(adata)
    
    # filter out sample swaps with wrong multiplexed genotype
    adata = adata[~adata.obs['Genotype'].isin(['WSBJ/CASTJ', 'AJ/129S1J', 'PWKJ/CASTJ'])].copy()
 
    adata.write(ofile)


############################################################################################################
#################################### Add metadata to cellbender output #####################################
############################################################################################################

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

def is_dummy(fname):
    """
    Check if a file is a dummy file (ie whether is empty)
    Returns True if it's a dummy, False if not
    """
    return os.stat(fname).st_size == 0

def touch_dummy(ofile):
    """
    Touch a dummy output file
    """
    open(ofile, 'a').close()
    
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


def add_meta_filter(filt_h5,
                    unfilt_adata,
                    genotypes,
                    wc,
                    bc_df,
                    kit,
                    chemistry,
                    sample_df,
                    ofile):

    """
    Format cellbender output adata:
    1. Add barcode and sample metadata
    2. Add cellbender counts as a layer
    3. Add raw counts as a layer
    4. Run scrublet within scanpy for each bc1_well using raw counts
    5. Create a cell ID to ensure uniqueness across subpools and experiments.
    6. If relevant, merge in the klue results. Update Genotype and Mouse_Tissue_ID for multiplexed wells. 
    7. Clean up adata obs
    """
    ############# 1. Add barcode and sample metadata #############
    print('Adding barcode and sample information...')
    adata = anndata_from_h5(filt_h5)
   
    adata.obs['bc'] = adata.obs.index
    adata.obs['bc1_sequence'] = adata.obs['bc'].apply(get_bc1)
    adata.obs['bc2_sequence'] = adata.obs['bc'].apply(get_bc2)
    adata.obs['bc3_sequence'] = adata.obs['bc'].apply(get_bc3)

    adata.obs['subpool'] = wc.subpool

    # merge in bc1 information
    temp = bc_df.copy(deep=True)
    temp = temp[['bc1_dt', 'well']].rename({'bc1_dt': 'bc1_sequence', 'well': 'bc1_well'}, axis=1)
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_sequence')

    # merge in bc2 and bc3 information
    for bc in [2,3]:
        bc_name = f'bc{bc}'
        seq_col = f'bc{bc}_sequence'
        well_col = f'bc{bc}_well'
        bc_df = get_bcs(bc, kit, chemistry)
        bc_df.rename({'well': well_col, bc_name: seq_col}, axis=1, inplace=True)
        adata.obs = adata.obs.merge(bc_df, how='left', on=seq_col)

    # merge in w/ sample-level metadata
    temp = sample_df.copy(deep=True)
    temp = temp.loc[(sample_df.plate==wc.plate)]
    adata.obs = adata.obs.merge(temp, how='left', on='bc1_well')
    
    # for now (to transfer info to/from the raw kallisto output)
    adata.obs_names = adata.obs['bc']
    
    # store gene names as a column in var
    adata.var['gene_name'] = adata.var_names
    
    # set gene ID as var names for further processing (avoid gene name duplication issues)
    adata.var_names = adata.var['gene_id']   
    
    ############# 2. Add cellbender counts as a layer #############
    # store cellbender normalized data as a layer
    adata.layers['cellbender_counts'] = adata.X.copy()
    
    ############# 3. Add raw counts as a layer #############
    # store raw data as a layer
    adata_raw = sc.read_h5ad(unfilt_adata)
    adata_raw.var_names = adata_raw.var['gene_id']   
    adata_raw_filtered = adata_raw[adata_raw.obs['bc'].isin(adata.obs['bc']), :].copy()
    adata_raw_filtered = adata_raw_filtered[adata.obs['bc'], :].copy()
    
    adata.layers['raw_counts'] = adata_raw_filtered.X.copy()
   
    ############# 4. Run scrublet within scanpy for each bc1_well #############
    # https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html
    # using raw counts
    print('Running scrublet...')
    
    adata_raw_filtered.obs['bc1_well'] = adata.obs['bc1_well']
    adata_raw_filtered.obs['bc1_well'] = adata_raw_filtered.obs['bc1_well'].astype(str)
    
    sc.pp.scrublet(adata_raw_filtered, batch_key="bc1_well", n_prin_comps=30)
    adata.obs['doublet_score'] = adata_raw_filtered.obs['doublet_score']
    adata.obs['predicted_doublet'] = adata_raw_filtered.obs['predicted_doublet']


    ############# 5. Create new cell IDs #############
    adata.obs['cellID'] = adata.obs['bc1_well']+'_'+\
        adata.obs['bc2_well']+'_'+\
        adata.obs['bc3_well']+'_'+\
        adata.obs['subpool']+'_'+\
        adata.obs['plate']
    adata.obs.reset_index(drop=True)

    # make sure these are unique + set as index
    assert len(adata.obs.index) == len(adata.obs.cellID.unique().tolist())
    adata.obs.set_index('cellID', inplace=True)
            
    # Add QC metrics from scanpy
    adata.var_names = adata.var['gene_name']
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], layer = 'raw_counts', percent_top=None, log1p=False, inplace=True)
    adata.obs.rename(columns={
        'n_genes_by_counts': 'n_genes_by_counts_raw',
        'total_counts': 'total_counts_raw',
        'total_counts_mt': 'total_counts_mt_raw',
        'pct_counts_mt': 'pct_counts_mt_raw'
    }, inplace=True)
    
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], layer = 'cellbender_counts', percent_top=None, log1p=False, inplace=True)
    adata.obs.rename(columns={
        'n_genes_by_counts': 'n_genes_by_counts_cb',
        'total_counts': 'total_counts_cb',
        'total_counts_mt': 'total_counts_mt_cb',
        'pct_counts_mt': 'pct_counts_mt_cb'
    }, inplace=True)

    
    ############# 6. Merge in the klue results #############
    # assign genotype for multiplexed wells
    if not is_dummy(genotypes):
        # filter out sample swaps with wrong multiplexed genotype
        adata = adata[~adata.obs['Genotype'].isin(['WSBJ/CASTJ', 'AJ/129S1J', 'PWKJ/CASTJ'])].copy()
        
        klue_counts_df = pd.read_csv(genotypes, sep='\t')
        klue_counts_df.set_index('cellID', inplace=True)

        # make sure we won't dupe any cols
        assert len(set(klue_counts_df.columns.tolist())&set(adata.obs.columns.tolist())) == 0

        # merge in first; this way we have access to the genotypes that should be in each well
        adata.obs = adata.obs.merge(klue_counts_df,how='left',left_index=True,right_index=True)

        # assign genotype for multiplexed wells
        df = adata.obs.copy(deep=True)
        df = assign_demux_genotype(df)

        print('Updating genotype for multiplexed wells...')
        # merge in w/ adata and replace old values in "Genotype" column for multiplexed wells with the klue results
        adata.obs = adata.obs.merge(df,how='left',left_index=True,right_index=True)
        
        # append "_klue_counts" to the genotype counts columns for clarity
        for col in klue_counts_df.columns:
            if col in adata.obs.columns:
                adata.obs.rename(columns={col: col + '_klue_counts'}, inplace=True)
        
        inds = adata.obs.loc[adata.obs.well_type=='Multiplexed'].index
        adata.obs.Genotype = adata.obs.Genotype.astype('str')
        adata.obs.loc[inds, 'Genotype'] = adata.obs.loc[inds, 'new_genotype']
        adata.obs.drop('new_genotype', axis=1, inplace=True)  
        
        genotype_list_founders = ['AJ', 'WSBJ', '129S1J', 'NODJ', 'PWKJ', 'NZOJ', 'CASTJ']
        genotype_list_f1s = ['B6AF1J', 'B6129S1F1J', 'B6NZOF1J', 'B6WSBF1J', 'B6NODF1J', 'B6PWKF1J', 'B6CASTF1J']

        if any(genotype in adata.obs['Genotype'].values for genotype in genotype_list_founders):
            ms1 = ['B6J', 'AJ', 'WSBJ', '129S1J']
            ms2 = ['NODJ', 'PWKJ', 'NZOJ', 'CASTJ']
        elif any(genotype in adata.obs['Genotype'].values for genotype in genotype_list_f1s):
            ms1 = ['B6J', 'B6AF1J', 'B6WSBF1J', 'B6129S1F1J']
            ms2 = ['B6NODF1J', 'B6PWKF1J', 'B6NZOF1J', 'B6CASTF1J']

        def update_mouse_tissue_id(row):
            if row['well_type'] == "Multiplexed" and row['Genotype'] != "tie":
                if row['Genotype'] in ms1:
                    return row['Multiplexed_sample1']
                elif row['Genotype'] in ms2:
                    return row['Multiplexed_sample2']
            return row['Mouse_Tissue_ID']
        
        obs = adata.obs
        obs['Mouse_Tissue_ID'] = obs.apply(update_mouse_tissue_id, axis=1)
        adata.obs['Mouse_Tissue_ID'] = obs['Mouse_Tissue_ID']
        
        # Update other sample-specific metadata columns
        multiplexed_obs = adata.obs[(adata.obs['well_type'] == "Multiplexed") & (adata.obs['Genotype'] != "tie")]

        sample_df_unique = sample_df.drop_duplicates(subset=['Mouse_Tissue_ID', \
                                                             'Tissue_weight_mg', \
                                                             'ZT', \
                                                             'Dissection_time', \
                                                             'Dissection_date', \
                                                             'Estrus_cycle', \
                                                             'Body_weight_g', \
                                                             'Age_days', \
                                                             'DOB'])
        
        merged = multiplexed_obs[['Mouse_Tissue_ID']].merge(sample_df_unique[['Mouse_Tissue_ID', \
                                                                              'Tissue_weight_mg', \
                                                                              'ZT', \
                                                                              'Dissection_time',\
                                                                              'Dissection_date', \
                                                                              'Estrus_cycle', \
                                                                              'Body_weight_g', \
                                                                              'Age_days', \
                                                                              'DOB']],
                                                            on='Mouse_Tissue_ID', how='left')
        
        merged.set_index(multiplexed_obs.index)
        adata.obs.update(merged.set_index(multiplexed_obs.index))
        
        
    ############# 7. Clean up obs #############
    adata.obs['lab_sample_id'] = adata.obs['Mouse_Tissue_ID']

    adata.obs.drop(columns=['mult_genotype_1', 'mult_genotype_2', 'mult_genotype'], inplace=True)
    adata.var = adata.var.drop(columns=['feature_type', 'genome'])

    
    if wc.plate == 'igvf_b01':
        adata.obs['subpool_type'] = np.where(adata.obs['subpool'].isin(['Subpool_7', 'Subpool_EX']), 'EX', 'NO')
    else:
        adata.obs['subpool_type'] = np.where(adata.obs['subpool'] == 'Subpool_EX', 'EX', 'NO')

    # add accessions if the file exists
    file_path = f"configs/{wc.plate.upper()}_barcode_sample_map.tsv"
    if os.path.exists(file_path):
        acc_map = pd.read_csv(file_path, sep="\t")
        acc_map = acc_map[acc_map['parse_barcode_type'] == "T"]

        acc_map = acc_map.drop(columns=['barcode', 'seqspec_region_id', 'parse_barcode_type'])
        adata_obs_df = adata.obs.reset_index()
        merged_df = adata_obs_df.merge(acc_map, left_on='bc1_well', right_on='well', how='left')
        merged_df.set_index('cellID', inplace=True)

        # Handle potential multiple values in 'sample' per well by collapsing them, comma-separated
        merged_df['sample'] = merged_df.groupby('cellID')['sample'].transform(lambda x: ','.join(x.dropna().unique()))

        adata.obs['sample'] = merged_df['sample']
    else:
        adata.obs['sample'] = "NA"
        

    # make all object columns string columns
    for c in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[c].dtype):
            adata.obs[c] = adata.obs[c].fillna('NA')
            
    keys = ["lab_sample_id","sample","plate",\
            "subpool","SampleType","Tissue","Sex",\
            "Age","Genotype","subpool_type","Protocol","Chemistry",\
            "bc","bc1_sequence","bc2_sequence","bc3_sequence","bc1_well","bc2_well","bc3_well",\
            "Row","Column","well_type",\
            "Mouse_Tissue_ID","Multiplexed_sample1", "Multiplexed_sample2","DOB",\
            "Age_days","Body_weight_g","Estrus_cycle",\
            "Dissection_date","Dissection_time","ZT",\
            "Dissector","Tissue_weight_mg",\
            "Notes","n_genes_by_counts_raw","total_counts_raw","total_counts_mt_raw",\
            "pct_counts_mt_raw","n_genes_by_counts_cb","total_counts_cb","total_counts_mt_cb",\
            "pct_counts_mt_cb","doublet_score","predicted_doublet","background_fraction",\
            "cell_probability","cell_size","droplet_efficiency"]
        
    if not is_dummy(genotypes):
        adata.obs = adata.obs[keys]
        print(adata.obs.head())
        print(adata.var.head())
        print(f'Saving file {ofile}')
        adata.write_h5ad(ofile)
    else:
        # for non-multiplexed experiments, remove the irrelevant multiplexed obs columns
        keys = [key for key in keys if key not in ["Multiplexed_sample1", "Multiplexed_sample2"]]

        adata.obs.drop(columns=['Multiplexed_sample1', 'Multiplexed_sample2'], inplace=True)
        
        adata.obs = adata.obs[keys]
        print(adata.obs.head())
        print(adata.var.head())
        
        print(f'Saving file {ofile}')
        adata.write_h5ad(ofile)
        

############################################################################################################
############################################## Final concat ################################################
############################################################################################################

def get_tissue_adatas(df, sample_df, wc, cfg_entry):

    # limit to input tissue
    temp_sample = sample_df.copy(deep=True)
    temp_sample = temp_sample.loc[temp_sample.Tissue==wc.tissue]

    # merge this stuff in with the fastq df
    fastq_df = df.copy(deep=True)
    temp = fastq_df.merge(temp_sample, on='plate', how='inner')

    # get the plate / subpool / sample info for this tissue
    plates = temp.plate.tolist()
    subpools = temp.subpool.tolist()
    samples = temp.Mouse_Tissue_ID.tolist()

    files = expand(cfg_entry,
           zip,
           plate=plates,
           sample=samples,
           subpool=subpools)
    files = list(set(files))

    return files


def concat_adatas(adatas, ofile):
    var_dfs = []
    print(adatas)
    
    for i, f in enumerate(adatas):
        temp = sc.read_h5ad(f)
        
        temp.var.drop(columns=['ambient_expression', 'cellbender_analyzed'], inplace=True)
        
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

def concat_adatas_raw_counts(adatas, ofile): 
    adata = None  # Initialize adata

    for i, f in enumerate(adatas):
        temp = sc.read_h5ad(f)
        del temp.layers['cellbender_counts']    
        temp.X = temp.layers['raw_counts'].copy()
        
        if adata is None:
            adata = temp
        else:
            adata = anndata.concat([adata, temp], join='outer', label = None, index_unique=None)

    adata.var = temp.var[['gene_id', 'gene_name', 'mt']]
    
    # Final save
    adata.write_h5ad(ofile)
        
        
def concat_adatas_cb_counts(adatas, ofile, temp_dir="cb_temp_files"):
    adata = None  # Initialize adata

    for i, f in enumerate(adatas):
        print(f)
        temp = sc.read_h5ad(f)
        del temp.layers['raw_counts']
        temp.X = temp.layers['cellbender_counts']

        if adata is None:
            adata = temp
        else:
            adata = anndata.concat([adata, temp], join='outer', label = None, index_unique=None)
            
    adata.var = temp.var[['gene_id', 'gene_name', 'mt']]
    
    # Final save
    adata.write_h5ad(ofile)