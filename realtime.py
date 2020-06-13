#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

#读取realtime结果excel文件，跳过前六行，把第一行为列名，并选取1，2，9列
rawdata = pd.read_excel(r'D:\yanglab_data\real_time\realtime.xls',skiprows=6,header=1,usecols=[1,2,9],names=['sample_name','target_name','cycle'])

#删除倒数5行
rawdata.drop(rawdata.index[-5:],inplace=True)

#获取所用的引物名字
target = rawdata['target_name'].unique()

#根据引物将数据分为几个不同的dataframe
data0 =  rawdata[rawdata['target_name'].isin([target[0]])]

#重命名列名cycle为actin_cycle
data0 = data0.rename({'cycle':'actin_cycle'},axis='columns')

#删除target_name一列
data0.drop(['target_name'],axis=1,inplace=True)

#重建索引
data0 = data0.reset_index(drop=True)

#得到引物为cbf的dataframe
data1 =  rawdata[rawdata['target_name'].isin([target[1]])]
data1 = data1.rename({'cycle':'cbf_cycle'},axis='columns')
data1.drop(['target_name'],axis=1,inplace=True)
data1 = data1.reset_index(drop=True)

#计算actin与cbf的差的cycle数
diff = pd.Series(data1['cbf_cycle'] - data0['actin_cycle'])

#获取dataframe的样本名
index = pd.Series(data0['sample_name'])

#将样本名和cycle差值合并为一个新的dataframe
data_used = pd.concat([index,diff],axis=1)

#重命名列名0为diffcyc
data_used.rename({0:'diffcyc'},axis='columns',inplace=True)

#获取不同重复的与对照相比后的cycle差
rep1 = []
for i in range(0,len(data_used),3):
    data = data_used['diffcyc'][i] - data_used['diffcyc'][0]
    rep1.append(data)

rep2 = []
for i in range(1,len(data_used),3):
    data = data_used['diffcyc'][i] - data_used['diffcyc'][1]
    rep2.append(data)

rep3 = []
for i in range(2,len(data_used),3):
    data = data_used['diffcyc'][i] - data_used['diffcyc'][2]
    rep3.append(data)

#将三个重复合并为一个dataframe
rep_all = pd.DataFrame({'rep1':rep1,'rep2':rep2,'rep3':rep3})

#获取索引
name_all = []
for i in range(0,len(data_used),3):
    data = data_used['sample_name'][i]
    name_all.append(data)

#修改rep_all的索引为样本名并排序
rep_all.index = pd.Series(name_all)
rep_all.sort_index(inplace=True)
#对每个样本的不同重复计算 2^(-n)，得到基因表达量
pow = 2 ** (-rep_all)
print(pow)

#pow.to_csv(r'D:\yanglab_data\qpcr_results.csv')