#### Figure 文件说明

```
bar_plot.py 为画柱状图的脚本

Broken_Axis.py 为对纵坐标进行截断的脚本

get_qpcr_reslut.py 包含一个自定义函数，可对qpcr数据直接分析出图
```

#### get_qpcr_reslut.py 参数说明

```

file_name='realtime.xls'： qpcr文件，文件中样本名称需设置的格式为 sample_name-treat,样本名在前，处理条件在后，否则该模块无法识别，

sample_name=[]： 样本的名称，用逗号隔开，需和qpcr的xls文件中的样本名称保持一致，如 ['col','m1','m2','m3']

treat=[]： 处理的条件名称，用逗号隔开,如 [0,1,2]

color_set=[]： 自定义设置颜色，颜色数目需和样本数目保持一致

target1='mock'： target1 必须为对照的引物，和xls文件保持一致

target2='condition'： 另一个引物名称

cut='no'： 默认为不对柱状图进行截断，修改为 'yes' 可进行截断

ax_setylim=100： 截断是上面部分默认为100

ax2_setylim=10： 截断是下面部分默认为10

```
