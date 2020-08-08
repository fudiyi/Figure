#!/usr/bin/R
#分组箱线图，读入的文件可由 get_AA_for_geom_boxplot.py 产生

library(ggplot2)
library(ggsignif)
library(ggpubr)

setwd("E:/")
data <-read.csv("4aa.csv",header=T)
data=data.frame(data)

# 修改柱子/x轴顺序
data$treat = factor(data$treat, levels=c('Warm','Cold'))
data$AA_type = factor(data$AA_type, levels=c('Asn','Asp','Gln','Glu'))

# 修改背景
mytheme<-theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major.y=element_line(color="grey",linetype=1),panel.grid.minor.y=element_line(color="grey",linetype=2))

pdf(file="4aa.pdf",width=8,height=6)

# 进行显著性分析
stat.test <- compare_means(Value ~ treat, data = data, group.by = "AA_type",method = "t.test", ref.group = "Warm")


ggplot(data,aes(x=AA_type,y=Value))+  
  geom_boxplot(aes(fill=treat),outlier.size=0.1)+
  scale_fill_manual(values=c("orange2","dodgerblue2"))+theme_classic() + theme(legend.position = c(0.9,0.3))+
  stat_pvalue_manual(stat.test, x = "AA_type",label = "p.signif",y.position = 7.5)

dev.off()


## Ref: [添加显著性标识的方法] (http://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html)
