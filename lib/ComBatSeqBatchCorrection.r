args<-commandArgs(TRUE)

myData <-read.delim(args[1], header=TRUE, sep="\t")
countMatrix <- myData[,c(3:ncol(myData))]
featureMatrix <- myData[,c(1:2)]

batch <- as.integer(c(strsplit(args[2], ",")[[1]]))

adjusted <- sva::ComBat_seq(countMatrix, batch=batch, group=NULL)

adjustedData <- merge(featureMatrix, adjusted, by=0)
adjustedData <- subset (adjustedData, select = -Row.names)
write.table(adjustedData, file = args[1], sep = "\t", row.names = FALSE)