import pandas as pd
from snakemake.io import expand

def parse_config(fname):
    df = pd.read_csv(fname, sep='\t')
    df['path'] = df.fastq.str.rsplit('/', n=2, expand=True)[0]+'/'
    df['path2'] = df.fastq.str.rsplit('/', n=1, expand=True)[0]+'/'
    df['r2_fastq'] = df.fastq.str.replace('_R1_', '_R2_')
    return df

def get_df_info(wc, df, col):
    temp = df.copy(deep=True)
    temp = temp.loc[(temp.plate==wc.plate)&\
                    (temp.subpool==wc.subpool)&\
                    (temp.lane==wc.lane)&\
                    (temp.run==int(wc.run))]
    import pdb; pdb.set_trace()
    assert len(temp.index) == 1
    return temp[col].values[0]

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
                        allow_missing=True),
                        plate=wc.plate,
                        subpool=wc.subpool)
        r2s = expand(expand(config['raw']['r2_fastq'],
                        zip,
                        lane=temp['lane'].tolist(),
                        allow_missing=True),
                        plate=wc.plate,
                        subpool=wc.subpool)
        fastq_str = ''
        for r1, r2 in zip(r1s, r2s):
            fastq_str+=f' {r1} {r2}'
        return fastq_str
