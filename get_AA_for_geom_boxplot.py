#!/usr/bin/python3

import pandas as pd
import numpy as np

data = pd.read_table('17aa_2019_12_17_results_row1_T.txt')
data.drop(columns=['Compound'],inplace=True)
data = np.log10(data+1)

warm = pd.DataFrame(data[:245])
warm['treat'] = 'Warm'
warm_index = list(warm.keys()) #获取列索引，warm.index()获取行索引

cold = pd.DataFrame(data[245:490])
cold['treat'] = 'Cold'
cold = cold.reset_index(drop=True)
cold.columns = warm_index


all_df1 = {}

for i in warm_index:

    aa = warm[[i,'treat']]
    aa1 = cold[[i,'treat']]

    all_df1['%s_df' % i] = pd.concat([aa,aa1])
    all_df1['%s_df' % i]['AA'] = i    
    all_df1['%s_df' % i].columns = ['Value','treat','AA_type']


data_last = pd.concat([all_df1['Asn_df'],all_df1['Asp_df'],all_df1['Gln_df'],all_df1['Glu_df']])
data_last.to_csv('4aa.csv',index=None)

