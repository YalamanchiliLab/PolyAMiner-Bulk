

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
res.geno <- results(dds, contrast = c("Genotype", "TR", "CR"))
res.geno$FC <- 2^res.geno$log2FoldChange

# log transformation #
rld <- rlogTransformation(dds)

## Order by adjusted p-value
res <- res.geno[order(res.geno$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
# print(resdata, nrow = 10)
# names(resdata)[1] <- "Symbol"

#Write results
# write.csv(resdata, file=paste0(argsv[3],"_DEG-results.txt"),row.names=F)
write.table(resdata, file = paste0(argsv[3],"_DEG-results.txt"), sep = "\t", row.names = FALSE)


