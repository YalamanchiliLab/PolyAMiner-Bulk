library(EnhancedVolcano)

args<-commandArgs(TRUE)

plab=c(unlist(strsplit(args[7], '_')),unlist(strsplit(args[8], '_')))



frame<-read.table(args[1], header=TRUE, sep="\t")
frame<-frame[is.finite(frame$padj),]
frame = frame[!duplicated(frame$Gene),]
# row.names(frame)=frame$Symbol
rownames(frame) = make.names(frame$Symbol, unique=TRUE)

keyvals <- ifelse(
  frame$padj >0.05, 'black',
  ifelse(frame$log2FoldChange > 0.263, 'red',
    ifelse(frame$log2FoldChange < -0.263, 'royalblue','black')))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'black'] <- 'Un'
names(keyvals)[keyvals == 'royalblue'] <- 'Dn'

tiff(paste0(args[2]),width=5, height=6, units='in',res=300)
  EnhancedVolcano(frame,
                lab = rownames(frame),
                # lab = as.character(frame$Symbol),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = plab,
                #selectLab = c(""),
                xlab = bquote(~Log[2]~ 'fold change'),
                FCcutoff = 0.263,
                pCutoff = 0.05,
                pointSize = 3,
                labSize = 3,
                labCol = 'black',
                #labFace = 'bold',
                boxedLabels = FALSE,
                colAlpha = 0.5,
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                colCustom=keyvals,
                ylim=c(0,25),
                xlim=c(-8,8),
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                title=args[9],subtitle=paste0("DEGs:",args[3]," FC 20%:",args[4]," UP:",args[5]," DN:",args[6]),caption='',ylab=bquote('-'~Log[10]~ 'Adj Pvalue')
                )
dev.off()