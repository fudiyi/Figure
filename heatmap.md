```shell
#!/usr/bin/R

setwd("E:/")
library(pheatmap)

# 读取数据
data<-read.table("16AA_content_for_heatmap.txt",header=T,row.names=1,sep="\t")
data1<-log10(data+1)

# 重命名行名为 'NIL1'
rownames(data) <- paste0("NIL", seq(nrow(data)))

# 设置色度条
bk = unique(c(seq(0,10, length=100)))

#pdf(file="file.pdf",width=4,height=6)

# 选取某几列数据
#data.cut1=data1[,c('Gln','Asp','Asn')]

# 构建行注释
ann_row = data.frame(
  treat = factor(rep(c("Warm", "Cold"), c(245,245)))
)
rownames(ann_row) = paste("NIL", 1:490, sep = "")

# 注释行名是否相同
all.equal(rownames(data1), rownames(ann_row))

pheatmap(data1,cluster_col=TRUE,cluster_row=FALSE,show_row=F,fontsize=8,color=colorRampPalette(c("blue","yellow","red"))(100),border_color=FALSE,
         annotation_row = ann_row, breaks=bk)

#dev.off()
```

Ref:

1. 超链接：[pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap)
