```shell
#!/usr/bin/R

setwd("E:/")
library(pheatmap)

# 读取数据
data<-read.table("16AA_content_for_heatmap.txt",header=T,row.names=1,sep="\t")
data1<-log10(data+1)

# 格式如下

              col_0_1     col_0_2     col_0_3     hy5_0_1     hy5_0_2     hy5_0_3   bbx78_0_1    bbx78_0_2   bbx78_0_3    col_8_1    col_8_2    col_8_3
AT5G08640 137.4938945 128.7763029 124.9079163  32.5018732  29.7941045  31.2528709  24.1846376  25.11270342  22.0918924 207.290228 197.301405 197.494511
AT3G51240 283.3414534 293.0110004 291.5439683  80.1877812  79.0644451  80.2063292  41.7151096  50.64761638  46.3740204 359.633906 369.037178 374.488689
AT2G22590  12.6206693  12.1731615  12.2334189   0.5862999   0.4531913   0.2337643   0.3485756   0.09405986   0.2283223  67.324538  71.937155  69.650895
AT3G55120 122.9868211 118.6099185 108.2246962  48.0817801  50.9628248  47.4541527  38.3075761  42.44671061  38.2799624 180.486677 176.795525 176.866845

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


# 格式如下
           Type Time
col_0_1     Col   t0
col_0_2     Col   t0
col_0_3     Col   t0
hy5_0_1     HY5   t0
hy5_0_2     HY5   t0
hy5_0_3     HY5   t0
bbx78_0_1   BBX   t0
bbx78_0_2   BBX   t0
bbx78_0_3   BBX   t0
col_8_1     Col   t8
col_8_2     Col   t8
col_8_3     Col   t8
hy5_8_1     HY5   t8
hy5_8_2     HY5   t8
hy5_8_3     HY5   t8
bbx78_8_1   BBX   t8
bbx78_8_2   BBX   t8
bbx78_8_3   BBX   t8
col_24_1    Col  t24
col_24_2    Col  t24
col_24_3    Col  t24
hy5_24_1    HY5  t24
hy5_24_2    HY5  t24
hy5_24_3    HY5  t24
bbx78_24_1  BBX  t24
bbx78_24_2  BBX  t24
bbx78_24_3  BBX  t24

# 注释行名是否相同
all.equal(rownames(data1), rownames(ann_row))

pheatmap(data1,cluster_col=TRUE,cluster_row=FALSE,show_row=F,fontsize=8,color=colorRampPalette(c("blue","yellow","red"))(100),border_color=FALSE,
         annotation_row = ann_row, breaks=bk,cutree_rows =5)

# 修改排序

annotation_col$Time = factor(annotation_col$Time, levels=c('t0','t8','t24'))

# 输出 cluster 文件
list=pheatmap(data1,scale = "row")
row_cluster=cutree(list$tree_row,k=10)

newOrder=data1[list$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder),names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
write.table(newOrder,file ="cluster.csv",sep =",",quote=FALSE)


#dev.off()
```

Ref:

1. 超链接：[pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap)
