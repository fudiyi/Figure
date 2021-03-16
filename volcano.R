library(gdata)
library(ggplot2)
setwd('E:/转录组/RNA-Seq/YH/hisat2_featurecounts_DEseq2/figure/火山图')
nrDEG = read.table("col24_NLP24.txt",row.names = 1,header = T)
nrDEG[is.na(nrDEG)] <- 1
nrDEG[nrDEG$pvalue <0.05 & nrDEG$log2FoldChange >1,ncol(nrDEG)+1]="Down"
nrDEG[nrDEG$pvalue <0.05 & nrDEG$log2FoldChange < -1,ncol(nrDEG)]="Up"
nrDEG[nrDEG$pvalue>=0.05 | 1 > abs(nrDEG$log2FoldChange),ncol(nrDEG)]="Normal"
colnames(nrDEG)[ncol(nrDEG)]="Regulate"
nrDEG$Regulate=factor(nrDEG$Regulate,levels = c("Up","Down","Normal"),order=T)
head(nrDEG)

#统计差异基因数目
length(which((nrDEG[,-1]=="Down")))
length(which((nrDEG[,-1]=="Up")))

pdf(file="col24_nlp24.pdf",width=10,height=10) 
col=c("red","lightseagreen", "black")
p_volcano = ggplot(nrDEG,aes(x=log2FoldChange,y=-log10(pvalue)))+
  geom_point(aes(color=nrDEG$Regulate),alpha=0.5)+scale_color_manual(values =col)+
  geom_hline(yintercept=c(-log10(0.05)),linetype=4)+
  geom_vline(xintercept=c(-log2(2),log2(2)),linetype=4)+
  theme_classic()+theme(legend.justification=c(1,0), legend.position=c(1,0.5))
plot(p_volcano)
dev.off()
