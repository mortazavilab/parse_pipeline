import scanpy as sc
import pandas as pd 

igvf = 'igvf_010'

sub_path = ['../' + igvf + '/aj_pwk/aj_pwk_concat.csv','../' + igvf + '/b6_nod/b6_nod_concat.csv','../' + igvf + '/wsb_nzo/wsb_nzo_concat.csv','../' + igvf + '/129_cast/129_cast_concat.csv']

merged_df = pd.DataFrame()

for sub in sub_path:
    df = pd.read_csv(sub)
    dfa=df.iloc[:,[1,2]]
    dfa['code'] = df['sublibrary'] + '-' + df['Barcode']
    if merged_df.empty:
        merged_df = dfa
    else:
        merged_df = merged_df.merge(dfa, left_on='code', right_on='code', how = 'outer')




merged_df.columns = merged_df.columns.str.replace(r'\.fa\.gz', '')
merged_df.columns = merged_df.columns.str.replace(r'/home/delaney/igvf/references/PRJNA923323/Mus_musculus_', '')




merged_df.to_csv('../'+igvf+'/demultiplex_matrix.csv')
