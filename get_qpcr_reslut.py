import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def get_qpcr_result(sample_name=[],treat=[],color_set=[],file_name='realtime.xls',target1='mock',target2='condition',cut='no',ax_setylim=100,ax2_setylim=10):

    rawdata = pd.read_excel(file_name,skiprows=6,header=1,usecols=[1,2,9],names=['sample_name','target_name','cycle']) # 读取数据
    rawdata.drop(rawdata.index[-5:],inplace=True) # 去掉最后没用的行

    target = sorted(list(rawdata['target_name'].unique())) # 获取引物的名字并排序，得到的列表中第一个必须是对照，若不是，将原始数据中的对照名字改成最小即可
    #print(type(target))

    if target1 == target[0] and target2 in target: # 判断列表第一个引物是否是对照，以及引物二是否在列表中
        actin_df = rawdata[rawdata['target_name'].isin([target[0]])] # 选取'target_name'为对照的数据存为新的dataframe
        actin_df = actin_df.rename({'cycle':'actin_cycle'},axis='columns') # 修改列名
        actin_df.drop(['target_name'],axis=1,inplace=True) # 删除'target_name'一列
        actin_df = actin_df.reset_index(drop=True) # 重置索引
        #print(actin_df)

        target2_df = rawdata[rawdata['target_name'].isin([target2])]
        target2_df = target2_df.rename({'cycle':target2},axis='columns')
        target2_df.drop(['target_name'],axis=1,inplace=True)
        target2_df = target2_df.reset_index(drop=True)
        #print(target2_df)    

        diff = pd.Series(target2_df[target2] - actin_df['actin_cycle']) # 与引物对照相比得到差异cycle数
        #print(diff)
        index1 = pd.Series(actin_df['sample_name']) # 获取样本名
        data_used = pd.concat([index1,diff],axis=1) # 合并样本名与差异cycle两个Series
        data_used.rename({0:'diffcyc'},axis='columns',inplace=True) # 重命名列名
        #print(data_used)

        rep1 = []  # 获取三个重复各自对应于样本对照的差异cycle数
        for i in range(0,len(data_used),3):
            data = data_used['diffcyc'][i] - data_used['diffcyc'][0] # 与样本的对照相比得到差异cycle，第一个重复
            rep1.append(data)

        rep2 = []
        for i in range(1,len(data_used),3):
            data = data_used['diffcyc'][i] - data_used['diffcyc'][1]
            rep2.append(data)

        rep3 = []
        for i in range(2,len(data_used),3):
            data = data_used['diffcyc'][i] - data_used['diffcyc'][2]
            rep3.append(data)
        
        rep_all = pd.DataFrame({'rep1':rep1,'rep2':rep2,'rep3':rep3})

        name_all = [] # 获取样本名字
        for i in range(0,len(data_used),3):
            data = data_used['sample_name'][i]
            name_all.append(data)

        rep_all.index = pd.Series(name_all) # 合并样本名与差异cycle
        rep_all.sort_index(inplace=True) # 按样本名重新排序

        pow_all = 2 ** (-rep_all) 
        #print(pow_all)

        mean = pow_all.mean(axis='columns') # 计算平均值与方差
        std = pow_all.std(axis='columns')
        data_plt = pd.concat([mean,std],axis=1)
        data_plt.columns=['mean','std']
        #print(data_plt)

        pow_all.to_csv('result.csv')


        mean = [] # 获得各个样本对应不同处理的基因表达平均值，存为列表
        error = [] # 方差

        for i in sample_name:
            
            l = []
            for j in treat:        
                l.append(str(i)+'-'+str(j))

            locals()[str(i)+'_mean'] = np.array(data_plt.loc[l]['mean']) # 批量命名变量并赋值
            locals()[str(i)+'_err'] = np.array(data_plt.loc[l]['std'])
            
            mean.append(locals()[str(i)+'_mean']) # 将变量对应的值添加到对应的列表中
            error.append(locals()[str(i)+'_err'])

            #print(l)
        #print(mean);print(error)
        #print(mut3_mean)

        n_groups = len(sample_name) # 样本数目
        total_width = 0.8 # 设置图的总宽度
        index = np.arange(len(treat)) # 柱子的横坐标位置
        bar_width = total_width/n_groups # 柱子的宽度
        opacity = 0.8 # 柱子透明度
        error_config = {'ecolor': 'black','capsize' :2} # 设置误差线的形式
        #print(index)


        all_bar = [] # 存放样本对应的柱子位置
        all_bar.append(index)
        for i in range(1,len(sample_name)): # 有多少个样本就有多少个index
            bar = index + bar_width*i # 每个柱子往右迁移一个样本宽度
            all_bar.append(bar)
        #print(all_bar)


        if cut == 'no': # 是否进行纵坐标截断

            fig, ax = plt.subplots()

            ax.spines['top'].set_visible(False) 
            ax.spines['right'].set_visible(False)

            for i in range(len(sample_name)): # 批量依次画样本的柱子

                ax.bar(all_bar[i], mean[i], bar_width,
                    alpha=opacity, color=color_set[i],
                    yerr=error[i], error_kw=error_config,
                    label=sample_name[i],edgecolor = 'white')            
        
            ax.set_xlabel('Conditon')
            ax.set_ylabel('Relative gene expression')
            ax.set_title(target2.upper(),loc='left')
            ax.set_xticks(index+bar_width*0.5*(n_groups-1))
            ax.set_xticklabels(treat)
            ax.legend(frameon=False,loc='upper left')
            plt.savefig('result.png',dpi=400)
            plt.show() 

        elif cut == 'yes':

            n = 3; m = 1 # 设置图形比例的数值
            gs = gridspec.GridSpec(2,1, height_ratios = [n,m],hspace=0.1) # 设置两个子图比例为3：1，间距为0.1

            plt.figure()

            ax = plt.subplot(gs[0,0:]) # 画第一个子图
            ax2 = plt.subplot(gs[1,0:], sharex = ax) # 画第一个子图，共享x轴

            ax.spines['top'].set_visible(False) # 使上边框不可见
            ax.spines['right'].set_visible(False)
            ax2.spines['right'].set_visible(False)

            ax_setylim2 = max(data_plt['mean']) + max(data_plt['std'])

            ax.set_ylim(ax_setylim, ax_setylim2)  # 设置第一个子图的范围
            ax2.set_ylim(0, ax2_setylim)  # 设置第二个子图的范围

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


            for i in range(len(sample_name)):

                ax.bar(all_bar[i], mean[i], bar_width,
                    alpha=opacity, color=color_set[i],
                    yerr=error[i], error_kw=error_config,
                    label=sample_name[i],edgecolor = 'white')

                ax2.bar(all_bar[i], mean[i], bar_width,
                    alpha=opacity, color=color_set[i],
                    yerr=error[i], error_kw=error_config,
                    label=sample_name[i],edgecolor = 'white')                


            ax2.set_xlabel('Condition')
            ax.set_ylabel('Relative gene expression')
            ax.set_title(target2.upper(),loc='left')
            ax.set_xticks(index + bar_width*0.5*(n_groups-1))
            ax.set_xticklabels(treat)
            ax.legend(frameon=False,loc='upper left')
            plt.savefig('result_cut.png',dpi=400)
            plt.show()

        else:

            print('please enter yes or no')


    else:
        print('The first target should be actin and target2 should in the xls')

