import scanpy as sc
import pandas as pd


igvf = '010'
pairs = ['aj_pwk','b6_nod','wsb_nzo','129_cast']
subs = ['2','3','4','5','6','7','8','9','10','11','12','13','14','15','16']

for pair in pairs:

    alldf = []

    for sub in subs:
        path = '../igvf_' + igvf + '/' + pair + '/' + 'igvf_'  +  igvf + '_sub' +  sub + '/counts_unfiltered_modified/adata.h5ad'
        print(path)
        a = sc.read_h5ad(path)
        adf = a.to_df()
        adf['sublibrary'] = 'sub' + sub
        adf['Barcode']=adf.index
        print('check')

        #if pair == 'aj_pwk':
        #    adf['pred_genotype'] = adf.apply(lambda row: "AJ" if row[0] > row[1] else "PWKJ", axis=1)
        #elif pair == 'b6_nod':
        #    adf['pred_genotype'] = adf.apply(lambda row: "B6J" if row[0] > row[1] else "NODJ", axis=1)
        #elif pair == 'wsb_nzo':
        #    adf['pred_genotype'] = adf.apply(lambda row: "WSBJ" if row[0] > row[1] else "NZOJ", axis=1)
        #elif pair == '129_cast':
        #    adf['pred_genotype'] = adf.apply(lambda row: "129S1J" if row[0] > row[1] else "CASTJ", axis=1)

        alldf.append(adf)
    
    concat_df = pd.concat(alldf, ignore_index=True)

    concat_df.to_csv('../igvf_' + igvf + '/' + pair + '/'  + pair+ '_concat.csv')
