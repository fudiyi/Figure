
library("plyr")
library("ggplot2")

setwd("E:/")
cabbage_exp <- read.table("S_for_geombar.txt",header = T)
cabbage_exp

### 文件格式 ###
treat	group	num	sample
h3	AS	1309	col
h24	AS	3373	c
h3	AS	2807	c
h24	AS	1800	col
h3	TF	1747	col
h24	TF	837	c
h3	TF	139	c
h24	TF	5311	col

ce <- ddply(cabbage_exp, "treat", transform,
            percent = num / sum(num) * 100) 

### 改变x轴标签的顺序
cabbage_exp$treat <- factor(cabbage_exp$treat,levels = c("h3","h24")) 

### 改变因子的名称
levels(cabbage_exp$group)[levels(cabbage_exp$group)=="AS"] <- "Alter"
levels(cabbage_exp$group)[levels(cabbage_exp$group)=="TF"] <- "Trans"

levels(cabbage_exp$sample)
cabbage_exp$sample = factor(cabbage_exp$sample, levels=c('col','c'))

ggplot(cabbage_exp,aes(x=treat,y=num,fill=sample)) +
  geom_bar(stat="identity",position = "stack")+theme_classic()+
  theme(
    strip.background = element_rect(
      color = "white", fill = "white"),
    panel.grid = element_blank(),
    panel.border = element_blank(),axis.line = element_line(size=1),
    axis.text = element_text(size=20),
    axis.ticks = element_line(size=1))

### 相关参数 ###

geom_bar() 函数用法如下：
geom_bar(mapping = NULL, data = NULL, stat = "count", width=0.9, position="stack")

stat：设置统计方法，有效值是count（默认值） 和 identity，其中，count表示条形的高度是变量的数量，identity表示条形的高度是变量的值；
position：位置调整，有效值是stack、dodge和fill，默认值是stack（堆叠），是指两个条形图堆叠摆放，dodge是指两个条形图并行摆放，fill是指按照比例来堆叠条形图，每个条形图的高度都相等，但是高度表示的数量是不尽相同的
width：条形图的宽度，是个比值，默认值是0.9
color：条形图的线条颜色
fill：条形图的填充色

theme() # 自定义图片背景

