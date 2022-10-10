library(EnhancedVolcano)
library(plyr)

args<-commandArgs(TRUE)

plab=c(unlist(strsplit(args[7], '_')),unlist(strsplit(args[8], '_')))



frame<-read.table(args[1], header=TRUE, sep="\t")
frame<-frame[is.finite(frame$PolyAIndex),]
frame = frame[!duplicated(frame$Gene),]
# row.names(frame)=frame$Symbol
rownames(frame) = make.names(frame$Symbol, unique=TRUE)
xlim_mag = round_any(max(abs(max(frame$PolyAIndex)), abs(min(frame$PolyAIndex))), 10) + 5

keyvals <- ifelse(
  frame$AdjG_Pval >0.05, 'grey',
  ifelse(frame$PolyAIndex >= 0.5, 'red',
    ifelse(frame$PolyAIndex <= -0.5, 'royalblue',
       ifelse((frame$PolyAIndex) > 0 & (frame$PolyAIndex < 0.5), 'salmon',
        ifelse((frame$PolyAIndex) < 0 & (frame$PolyAIndex > -0.5), 'skyblue','black')))))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'red'] <- 'SigEl'
names(keyvals)[keyvals == 'salmon'] <- 'El'
names(keyvals)[keyvals == 'grey'] <- 'NoSig'
names(keyvals)[keyvals == 'royalblue'] <- 'SigSh'
names(keyvals)[keyvals == 'skyblue'] <- 'Sh'


print(colnames(frame))
tiff(paste0(args[2]),width=6, height=7, units='in',res=200)
#svg(paste0(args[2],"_Vol_0.1.svg"),width = 6, height = 7)
  EnhancedVolcano(frame,
                #lab = rownames(frame),
                lab = as.character(frame$Symbol),
                x = 'PolyAIndex',
                y = 'AdjG_Pval',
                #selectLab = plab,
                # selectLab = c("Ctnnb1","Fat1","Scoc","Usp12","Zswim6","Mob1b","Eaf1","Sde2","Ndufa6","Lactb2","Araf","Itga5","Crebzf","Golga5","Ldb1"),
                #xlab = bquote(~Log[2]~ 'fold change'),
                xlab = "PolyAIndex",
                FCcutoff = 0.5,
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
                widthConnectors = 0.5,
                colConnectors = 'black',
                colCustom=keyvals,
                ylim=c(0,15),
                xlim=c(-xlim_mag,xlim_mag),
                #title=args[2],
                title=args[9],
                subtitle=paste0("APAs:",args[3]," PAIndex:",args[4]," El:",args[5]," St:",args[6]),
                caption='',
                gridlines.major=FALSE,
                gridlines.minor=FALSE,
                ylab=bquote('-'~Log[10]~ 'Adj Pvalue')
                )
dev.off()