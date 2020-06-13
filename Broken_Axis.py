import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

################################ 整理数据 ##############################

data = pd.read_csv(r'D:\yanglab_data\real_time\qpcr_results.csv',index_col=0) # 读入数据时把第一列作为索引
mean = data.mean(axis='columns') # 求平均值
std = data.std(axis='columns') # 求标准差
data_plt = pd.concat([mean,std],axis=1) # 合并
data_plt.columns=['mean','std'] # 重命名列名


ref = np.array(data_plt.loc[['Ref0','Ref2','Ref4','Ref8']]['mean']) # 读取ref数据并存为列表
ref_err = np.array(data_plt.loc[['Ref0','Ref2','Ref4','Ref8']]['std']) # 读取ref数据的靶值
mut1 = np.array(data_plt.loc[['mut1-0','mut1-2','mut1-4','mut1-8']]['mean'])
mut1_err = np.array(data_plt.loc[['mut1-0','mut1-2','mut1-4','mut1-8']]['std'])
mut2 = np.array(data_plt.loc[['mut2-0','mut2-2','mut2-4','mut2-8']]['mean'])
mut2_err = np.array(data_plt.loc[['mut2-0','mut2-2','mut2-4','mut2-8']]['std'])
mut3 = np.array(data_plt.loc[['mut3-0','mut3-2','mut3-4','mut3-8']]['mean'])
mut3_err = np.array(data_plt.loc[['mut3-0','mut3-2','mut3-4','mut3-8']]['std'])
print(mut3_err)
print(mut3)
################################ 设置图形格式 #################################

n = 3; m = 1 # 设置图形比例的数值
gs = gridspec.GridSpec(2,1, height_ratios = [n,m],hspace=0.1) # 设置两个子图比例为3：1，间距为0.1

plt.figure()

ax = plt.subplot(gs[0,0:]) 
ax2 = plt.subplot(gs[1,0:], sharex = ax)

ax.spines['top'].set_visible(False) # 使上边框不可见
ax.spines['right'].set_visible(False)
ax2.spines['right'].set_visible(False)

ax.set_ylim(100, 6500)  # 设置第一个子图的范围
ax2.set_ylim(0, 6)  # 设置第二个子图的范围

ax.spines['bottom'].set_visible(False) # 截断处的横线不显示
ax2.spines['top'].set_visible(False)
ax.axes.get_xaxis().set_visible(False) # 使第一个子图的x轴不可见
ax.tick_params(labeltop=False)  # 最上面不显示刻度
ax2.xaxis.tick_bottom() # Move ticks and ticklabels (if present) to the bottom of the axes


d = .015 # 设置截断线段的长度
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)

on = (n+m)/n; om = (n+m)/m
ax.plot((-d,d),(-d*on,d*on), **kwargs) # top-left diagonal
#ax.plot((1-d*on,1+d*on),(-d*on,d*on), **kwargs) # top-right diagonal
kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
ax2.plot((-d,d),(1-d*om,1+d*om), **kwargs) # bottom-left diagonal
#ax2.plot((1-d*on,1+d*on),(1-d*om,1+d*om), **kwargs) # bottom-right diagonal

############################## 绘图 ###############################

n_groups = 4 # 样本的数目（ref，mut1-3）
index = np.arange(n_groups) # 转换为列表
bar_width = 0.2 # 设置柱子的宽度
opacity = 0.8 # 设置柱子的透明度
error_config = {'ecolor': 'black','capsize' :2} # 设置误差线为黑色，宽度为2

############ 画第一个子图的柱状图 ############

# 画第一个样本的柱子
rects1 = ax.bar(index, ref, bar_width,
                alpha=opacity, color='peachpuff',
                yerr=ref_err, error_kw=error_config,
                label='ref',edgecolor = 'white')
# 画第二个样的柱子
rects2 = ax.bar(index + bar_width, mut1, bar_width,
                alpha=opacity, color='tomato',
                yerr=mut1_err, error_kw=error_config,
                label='mut1',edgecolor = 'white')

rects3 = ax.bar(index + bar_width*2, mut2, bar_width,
                alpha=opacity, color='orangered',
                yerr=mut2_err, error_kw=error_config,
                label='mut2',edgecolor = 'white')

rects4 = ax.bar(index + bar_width*3, mut3, bar_width,
                alpha=opacity, color='red',
                yerr=mut3_err, error_kw=error_config,
                label='mut3',edgecolor = 'white')

############ 画第二个子图的柱状图 ############

rects1 = ax2.bar(index, ref, bar_width,
                alpha=opacity, color='peachpuff',
                yerr=ref_err, error_kw=error_config,
                label='ref',edgecolor = 'white')

rects2 = ax2.bar(index + bar_width, mut1, bar_width,
                alpha=opacity, color='tomato',
                yerr=mut1_err, error_kw=error_config,
                label='mut1',edgecolor = 'white')

rects3 = ax2.bar(index + bar_width*2, mut2, bar_width,
                alpha=opacity, color='orangered',
                yerr=mut2_err, error_kw=error_config,
                label='mut2',edgecolor = 'white')

rects4 = ax2.bar(index + bar_width*3, mut3, bar_width,
                alpha=opacity, color='red',
                yerr=mut3_err, error_kw=error_config,
                label='mut3',edgecolor = 'white')


ax2.set_xlabel('Hours in 4 ℃')
ax.set_ylabel('Relative gene expression')
#ax.set_title('Creat by python')
ax.set_xticks(index + bar_width*1.5)
ax.set_xticklabels(('0', '2', '4', '8'))
ax.legend(frameon=False,loc='upper left')


plt.show()