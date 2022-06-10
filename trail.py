from bioinfokit import analys, visuz



# load dataset as pandas dataframe
df = analys.get_data('volcano').data
df.head(2)
          GeneNames  value1  value2    log2FC       p-value
0  LOC_Os09g01000.1    8862   32767 -1.886539  1.250000e-55
1  LOC_Os12g42876.1    1099     117  3.231611  1.050000e-55

visuz.GeneExpression.volcano(df=df, lfc='log2FC', pv='p-value')


df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] < pv_thr[0]), 'color_add_axy'] = color[0]  # upregulated
df.loc[(df[lfc] <= -lfc_thr[1]) & (df[pv] < pv_thr[1]), 'color_add_axy'] = color[2]  # downregulated
df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate
df['logpv_add_axy'] = -(np.log10(df[pv]))