######## TCGA dataset preprocessing
raw<-read.csv("TCGA-BRCA-TPM.csv",row.names = 1)
raw<-log2(raw+1)
library(limma)
raw<-normalizeBetweenArrays(raw)
LNCS<-read.csv("Human-lncRNAs-all.csv")

r<-t(data.frame(strsplit(rownames(raw),"[.]")))
rownames(raw)<-r[,1]
# change ESEMBLE id into SYMBOL id
library(clusterProfiler)
ids<-bitr(rownames(raw),fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db',drop = TRUE)
ids<-ids[!duplicated(ids$SYMBOL),]
ids<-ids[!duplicated(ids$ENSEMBL),]
exprSet<-raw
exprSet<-exprSet[which(rownames(exprSet) %in% ids$ENSEMBL),]
rownames(exprSet)<-ids[match(rownames(exprSet),ids$ENSEMBL),"SYMBOL"]
write.csv(exprSet,file = "TCGA-BRCA-normSYMBOL.csv")
rm(list=ls())

