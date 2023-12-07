import pandas as pd
import numpy as np
from snakemake.io import expand

def parse_config(fname):
    df = pd.read_csv(fname, sep='\t')
    df['path'] = df.fastq.str.rsplit('/', n=2, expand=True)[0]+'/'
    df['path2'] = df.fastq.str.rsplit('/', n=1, expand=True)[0]+'/'
    # df['r2_fastq'] = df.fastq.str.replace('_R1_', '_R2_')
    # df['fastq_pairs'] = df.apply(lambda x: (x.fastq, x.fastq_r2) \
    #                                 if not pd.isna(x.fastq_r2) else (x.fastq),
    #                                 axis=1)
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
    df.loc[inds, 'mult_genotype'] = df.loc[inds, g_cols[0]]+'_'+\
                                   df.loc[inds, g_cols[1]]

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
            entry = 'r1_fastq'
        elif read == 'R2':
            entry = 'r2_fastq'
        return expand(expand(config['raw'][entry],
                        zip,
                        lane=temp['lane'].tolist(),
                        run=temp['run'].tolist(),
                        allow_missing=True),
                        plate=wc.plate,
                        subpool=wc.subpool)

    elif how == 'str':
        r1s = expand(expand(config['raw']['r1_fastq'],
                        zip,
                        lane=temp['lane'].tolist(),
                        run=temp['run'].tolist(),
                        allow_missing=True),
                        plate=wc.plate,
                        subpool=wc.subpool)
        r2s = expand(expand(config['raw']['r2_fastq'],
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
