
# DEG analysis -DESeq2 #
rm(list=ls(all=TRUE))
argsv = commandArgs(trailingOnly=TRUE)
library("gplots") 
library("ggplot2")
library("Rtsne")
require("lattice")
library("GGally")
library('ggrepel')
library("DESeq2")
library('umap')
library('calibrate')
library('RColorBrewer')

# Load sample information # 
colData<- read.table(argsv[2], header=TRUE, sep="\t", row.names=NULL)
#colData<-colData[sort(colData$X),]
names<-colData$Name
Genotype<-colData$Genotype
count.matrix <- read.table(argsv[1], header=TRUE, sep="\t", row.names=1,quote = "",check.names=FALSE)
#count.matrix<-count.matrix[,sort(colnames(count.matrix))]
dds = DESeqDataSetFromMatrix(countData =round(count.matrix), colData = colData, design = ~ Genotype)
dds <- estimateSizeFactors(dds)
dat <- counts(dds, normalized = TRUE)
dds <- DESeq(dds, betaPrior = FALSE)
# Counts Log transform #
tt = counts(dds,normalized = TRUE)
tt<-log(tt+1)
tt<-t(tt)
names<-row.names(tt)


# PCA plot #
pca<-prcomp(tt)
PCi<-data.frame(pca$x,S=names)
# proportion of variance
variance = pca$sdev^2 / sum(pca$sdev^2)
variance = round(variance, 3) * 100
allpc=paste(as.character(variance),collapse=", ")
tiff(paste0(argsv[3],".PCA.tif"),width=8, height=6, units='in',res=300)
ggplot(PCi,aes(x=PC1,y=PC2))+
  geom_point(aes(color=Genotype),alpha=1,size=3)+
  #geom_text(aes(label=names),hjust=0.5, vjust=-1)+
  theme_bw()+
  scale_x_continuous(name = paste(c("PC1 :",variance[1],"%" ), collapse = ""))+
  scale_y_continuous(name = paste(c("PC2 :",variance[2],"%" ), collapse = ""))+
  geom_text_repel(label=names)+
  ggtitle(paste(argsv[3],c(" PCs:",allpc),collapse=" "))+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0))
dev.off()

# t-SNE #
# Correlatin matrix# 
cm<-cor(t(tt), method = c("spearman"))
set.seed(2020)
tsne_out <- Rtsne(tt,  perplexity = 1,  max_iter=500,is_distance = F) # Run TSNE 
scores = data.frame('Y1' = tsne_out$Y[,1], 'Y2' = tsne_out$Y[,2])
tiff(paste0(argsv[3],".t-SNE.tif"),width=8, height=6, units='in',res=300)
ggplot(scores, aes(x=Y1, y=Y2)) + 
  geom_point(size=2,aes(color=Genotype),alpha=0.5)+
  #geom_text(aes(label=names),hjust=0, vjust=0)+
  #theme_classic()+
  theme_bw()+
  theme(axis.title=element_text(size=12),text= element_text(size=12))+
  scale_x_continuous(name = "t-sne")+
  scale_y_continuous(name = "t-sne")+
  geom_text_repel(label=names)+
  ggtitle(paste0(argsv[3],": t-SNE plot"))+
  theme(plot.title = element_text(hjust = 0.5,vjust = 0))
dev.off()