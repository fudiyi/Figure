# R作图脚本汇总

#### 以前使用过的一些R作图脚本，未进行细致整理，以后绘制相应图时逐一更新

#### ggviolin.r

```R
#########ggviolin#######

setwd("E:/R/SYT")  #设置文件存放路径
library(ggpubr)   #载入包
data<-read.table("2019_12_9.txt",header=T) #读取文件
data$variable <- factor(data$variable,levels = c("warm","cold")) #改变x轴数据顺序
table(data$variable)  #查看顺序

#data form
type variable       value
min     warm  0.00000000
min     warm  0.63454628
min     warm -1.30174159
min     warm -0.02093827
min     warm  1.81735627
min     warm  0.32769434


pdf=(file="file.pdf",width=8,height=6) #导出pdf
ggviolin(data, "variable","value", fill = "type",palette = c("cyan", "orange"), add = "boxplot", add.params = list(color = "black"))  #作图
dev.off() 
```

#### GO.r

```R
##barplot

data_bar <- read.table("col_up_mac.txt",header = T)
data$Pathway <- factor(data$Pathway,levels = c("Molecular_Function","Biological_Process"))
ggbarplot(data_bar,x="Name",y="Score",fill="Pathway",rotate=T,palette = c("#ADC9D7","#E39E3E"),xlab = "",ylab="-log10(FDR)",sort.val = "asc",group="Pahway",sort.by.groups = TRUE)+theme_base()

```

#### CMplot.r

```R
 #CMplot
 
BiocManager::install("CMplot")
library(CMplot)
setwd("E:/mGWAS/2019_11_7/cold")
ser_chr3 <- read.table("ser_chr3.txt")
ser_chr3[1:4,1:4]

##############data#############
     V1            V2 V3     V4
1 pheno  chr3.S_34582  3  34582
2 pheno  chr3.S_34593  3  34593
3 pheno chr3.S_223342  3 223342
4 pheno chr3.S_223360  3 223360

ser_chr3_cmplot <- ser_chr3[,c(2,3,4,7)]

##############cut_col##########
             V2 V3     V4      V7
1  chr3.S_34582  3  34582 0.25171
2  chr3.S_34593  3  34593 0.44181
3 chr3.S_223342  3 223342 0.44337
4 chr3.S_223360  3 223360 0.52940
5 chr3.S_224171  3 224171 0.34006
6 chr3.S_224270  3 224270 0.20247

colnames(ser_chr3_cmplot) <- c( 'SNP', 'chr','pos', 'P')

###############################
            SNP chr    pos       P
1  chr3.S_34582   3  34582 0.25171
2  chr3.S_34593   3  34593 0.44181
3 chr3.S_223342   3 223342 0.44337
4 chr3.S_223360   3 223360 0.52940
5 chr3.S_224171   3 224171 0.34006
6 chr3.S_224270   3 224270 0.20247

CMplot(all_chr3,chr.den.col=c("darkgreen", "yellow","red"),chr.labels=c("Lys","Arg","Asn","Ser"),threshold=c(1e-5,1e-6),threshold.col=c("red","orange"),threshold.lwd=1, threshold.lty=c(1,2), amplify= TRUE,signal.col="red")


```

#### corrplot.r

```R
#corrplot

library(corrplot)
M <- cor(mtcars)
res1 <- cor.mtest(mtcars, conf.level = .95)
corrplot(M,method=c("color"),tl.cex=0.1,tl.col="black")
corrplot(M, p.mat = res1$p, insig = "label_sig",sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")
```

#### heatmap.r

```R
setwd("E:/RNA-Seq/XSW")
library(pheatmap)
data<-read.table("349_gene.txt",header=T,row.names=1,sep="\t")
pdf(file="file.pdf",width=4,height=6)
pheatmap(data,scale="row",cluster_col=TRUE/FALSE,cluster_row=TRUE/FALSE,show_row/colnames=T/F,fontsize=x,breaks= bk,color=colorRampPalette(c("blue","yellow","red"))(100),border_color=FALSE)
dev.off()

# 构建列注释信息 
annotation_col = data.frame( CellType = factor(rep(c("CT1", "CT2"), 5)), Time = 1:5 ) 
rownames(annotation_col) = paste("Test", 1:10, sep = "") 
head(annotation_col)

# 设置色度条在0-5
bk = unique(c(seq(0,5, length=100)))

```

#### PCA.r

```R
setwd("D:\\Users\\PCA")
library(ggplot2)
cpd=read.table("compound2.txt",sep="\t",header=T,row.names=1)
cpd2=log2(cpd+1)
cpd2=data.frame(cpd2)
class.info=read.table("class.txt",sep="\t",header=T)
class.info=class.info[class.info$class!="",]

index=match(as.character(class.info$ID),names(cpd2))
index2=index[!is.na(index)]
cpd3=cpd2[,index2]

index3=match(names(cpd3),class.info$ID)
class.info=class.info[index3,]


cpd_pr=princomp(cpd3,cor=TRUE,scale=T)
loadings=cpd_pr$loadings

PCs=loadings[,1:2]
PCs=data.frame(PCs)
PCs$class2=class.info$class2
PCs$class=class.info$class

ggplot(PCs,aes(x=Comp.1,y=Comp.2,color=class,shape=class2)) + 
  #geom_point(size=4,col="white") +
  geom_point(size=4) +
  theme(axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
		axis.text.x=element_text(size=8),
		axis.text.y=element_text(size=8)) +
  #scale_color_manual(values=c("red","green","blue"))+
  #xlim(-0.37,0.38) +
  #ylim(-0.45,0.45) +
  #xlab("PC1 (16.3%)") +
  #ylab("PC2 (9.2%)") +
  theme(#panel.grid.major = element_line(size=1,linetype =1,color="w"), 
        panel.grid.minor = element_blank())
        #panel.background = element_rect(fill = "lightblue"))
		#axis.line = element_line(colour = "black"))

```

#### qqman.r

```R
#qqman
asp_warm<-read.table("5_chr3.txt",header=T,sep="\t")
hig = read.table("152908_genebody.txt")
> pdf(file="asp_warm_chr3.pdf",width=6,height=4)
> manhattan(asp_warm,genomewideline = F,suggestiveline = -log10(1e-06),col = c("gray"), highlight = hig$V1,ylim = c(0, 10))
> dev.off()

# hig文件格式              
 V1
1 chr9.S_122220227
2 chr9.S_122220238
3 chr9.S_122220286
4 chr9.S_122220287
5 chr9.S_122220393
6 chr9.S_122220419


```

#### upset.r

```R
between <- function(row, min, max){
     newData <- (row["RI_in_mac"] < max) & (row["RI_in_mac"] > min)
}

	upset(data_all_RI,sets=c("RI_0_3","SE_0_3","A3SS_0_3","A5SS_0_3","MXE_0_3","RI_3_24","SE_3_24","A3SS_3_24","A5SS_3_24","MXE_3_24","RI_0_24","SE_0_24","A3SS_0_24","A5SS_0_24","MXE_0_24"),order.by = "freq",keep.order = TRUE,mainbar.y.label = "Gene Intersections", sets.x.label = "Splicing Form", mb.ratio = c(0.6, 0.4), 
		queries = list(list(query = intersects, params = list("RI_0_3")),list(query = between, params=list(1994,1996), color="red", active=T),
		list(query = intersects, params = list( "RI_0_24")),list(query = between, params=list(1995,1997), color="red", active=T),
		list(query = intersects, params = list("RI_3_24")),list(query = between, params=list(1996,1998), color="red", active=T),
		list(query = intersects, params = list("RI_0_3", "RI_0_24")),list(query = between, params=list(1997,1999), color="red", active=T),
		list(query = intersects, params = list("RI_0_3", "RI_3_24")),list(query = between, params=list(1998,2000), color="red", active=T),
		list(query = intersects, params = list("RI_0_24", "RI_3_24")),list(query = between, params=list(1999,2001), color="red", active=T),
		list(query = intersects, params = list("RI_0_3", "RI_0_24","RI_3_24")),list(query = between, params=list(2000,2002), color="red", active=T)))

for skip

upset(data,sets=c("RI_0_3","SE_0_3","A3SS_0_3","A5SS_0_3","MXE_0_3","RI_3_24","SE_3_24","A3SS_3_24","A5SS_3_24","MXE_3_24","RI_0_24","SE_0_24","A3SS_0_24","A5SS_0_24","MXE_0_24"),order.by = "freq",keep.order = TRUE,mainbar.y.label = "Gene Intersections", sets.x.label = "Splicing Form", mb.ratio = c(0.6, 0.4), 
+       queries = list(list(query = intersects, params = list("RI_0_3")),list(query = between, params=list(1999,2001), color="red", active=T),
+                      list(query = intersects, params = list( "RI_0_24")),list(query = between, params=list(1997,1999), color="red", active=T),
+                      list(query = intersects, params = list("RI_3_24")),list(query = between, params=list(2000,2002), color="red", active=T),
+                      list(query = intersects, params = list("RI_0_3", "RI_0_24")),list(query = between, params=list(1994,1996), color="red", active=T),
+                      list(query = intersects, params = list("RI_0_3", "RI_3_24")),list(query = between, params=list(1998,2000), color="red", active=T),
+                      list(query = intersects, params = list("RI_0_24", "RI_3_24")),list(query = between, params=list(1995,1997), color="red", active=T),
+                      list(query = intersects, params = list("RI_0_3", "RI_0_24","RI_3_24")),list(query = between, params=list(1996,1998), color="red", active=T)))


##upsetR##
upset(movies, queries = list(list(query = intersects, params = list("Drama", 
    "Comedy", "Action"), color = "orange", active = T), list(query = intersects, 
    params = list("Drama"), color = "red", active = F), list(query = intersects, 
    params = list("Action", "Drama"), active = T)))
```

#### box.r

```R
setwd("E:/R/fdy/realtime")
	
library(ggplot2)
library(ggsignif)
pdf(file="warm.pdf",width=4,height=4)
my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
ggplot(data,aes(x=type,y=content))+geom_boxplot(fill=c("blue","red"))+geom_signif(comparisons=list(c("cold","warm")),test=t.test,map_signif_level=TRUE)
dev.off()
#单个组
#color cyan

> data<-read.table("all1.txt",header=T)
> pdf(file="fdy6.pdf",width=8,height=4)
> ggplot(data,aes(x=type,y=value,color=variable))+geom_boxplot(aes(color=variable),outlier.size=0.1)+scale_color_manual(values=c("blue","red"))
> dev.off()
#多个组

> data<-read.table("all1.txt",header=T)
> pdf(file="last2.pdf",width=8,height=4)
> ggplot(data,aes(x=type,y=value))+geom_boxplot(aes(fill=variable),outlier.size=0.1)+scale_fill_manual(values=c("blue","red"))
> dev.off()
#全部填充

ggplot(data=data,aes(x=few,y=value,fill=type))+geom_bar(position=position_dodge(), stat="identity")+geom_errorbar(aes(ymax = value+std, ymin = value-std),position = position_dodge(0.9), width = 0.2)+facet_grid(.~variable)+theme_wsj()+scale_fill_wsj("rgby","")
#柱状图
data$few= factor(data$few, levels=c('col','a_myc','c_myc','ko_1')) 
#x轴排序
dt$obj = factor(dt$obj, levels=c('D','B','C','A','E')) 
## 设置柱条的顺序

> pdf(file="AA_content.pdf",width=8,height=4,bg="white")
> mytheme<-theme(panel.background=element_rect(fill="white",color="black"),panel.grid.major.y=element_line(color="grey",linetype=1),panel.grid.minor.y=element_line(color="grey",linetype=2))
> ggplot(data,aes(x=type,y=value,color=variable))+geom_boxplot(aes(color=variable),outlier.size=0.1)+scale_color_manual(values=c("blue","orange"))+mytheme
> dev.off()
#修改背景


# Grouped bar plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggpubr)
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
# Comparisons against reference
stat.test <- compare_means(len ~ dose, data = ToothGrowth, group.by = "supp",method = "t.test", ref.group = "0.5")
#> Warning: `cols` is now required.
#> Please use `cols = c(p)`
stat.test
#> # A tibble: 4 x 9
#>   supp  .y.   group1 group2            p      p.adj p.format p.signif method
#>   <fct> <chr> <chr>  <chr>         <dbl>      <dbl> <chr>    <chr>    <chr> 
#> 1 VC    len   0.5    1      0.000000681  0.000002   6.8e-07  ****     T-test
#> 2 VC    len   0.5    2      0.0000000468 0.00000019 4.7e-08  ****     T-test
#> 3 OJ    len   0.5    1      0.0000878    0.000088   8.8e-05  ****     T-test
#> 4 OJ    len   0.5    2      0.00000132   0.0000026  1.3e-06  ****     T-test
# Plot
bp <- ggbarplot(ToothGrowth, x = "supp", y = "len",fill = "dose", palette = "jco",add = "mean_sd", add.params = list(group = "dose"),position = position_dodge(0.8))
bp + stat_pvalue_manual(stat.test, x = "supp", y.position = 33,label = "p.signif",position = position_dodge(0.8))

hap_data$hap <- as.factor(hap_data$hap)
stat.test <- compare_means(value ~ hap, data = hap_data, group.by = "AA_type",method = "t.test", ref.group = "A")
ggplot(data,aes(x=AA_type,y=value))+geom_boxplot(aes(fill=hap),outlier.size=0.1)+scale_fill_manual(values=c("cyan","orange","green"))+mytheme+ stat_pvalue_manual(stat.test, x = "AA_type", y.position = 8,label = "p.signif",position = position_dodge(0.8))


#小提琴tu
ggviolin(data, x = "hap", y = "value", fill = "hap",palette = c("cyan", "orange", "green"),add = "boxplot", add.params = list(fill = "white"))
ggviolin(data, x = "hap", y = "value", fill = "hap",palette = c("cyan", "orange", "green"),add = "boxplot", add.params = list(fill = "white"),ylim=c(-3,4.5))+stat_compare_means(comparisons = my_comparisons, label = "p.signif",label.y = c(3,3.5,4))

```

