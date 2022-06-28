############# Reload the RNA expression matrix dataset ######################

raw<-read.csv("TCGA-BRCA-normSYMBOL.csv",row.names = 1)
LNC<-read.table('../Human-lncRNAs-all.txt')
exprSet<-raw

normal<-data.frame(exprSet[,grep('.11A.',colnames(exprSet))])
tumor<-exprSet[,grep('.01',colnames(exprSet))]

group_list<-c(rep('T',ncol(tumor)),rep('N',ncol(normal)))
new_exprSet<-cbind(tumor,normal)

modType=c(group_list) 
library(limma)
design <- model.matrix(~0+factor(modType))
colnames(design)=levels(factor(modType))
rownames(design)=colnames(new_exprSet)
design
#contrast.matrix<-makeContrasts(paste0(unique(modType),collapse = "-"), levels = design)
contrast.matrix<-makeContrasts("T-N", levels = design) 
contrast.matrix 
deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)  ## default no trend !!!
  tempOutput = topTable(fit2, coef=1, n=Inf)
  #nrDEG = na.omit(tempOutput) #È¥µô¿ÕÖµ
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  #head(nrDEG)
  return(tempOutput)}
re = deg(new_exprSet,design,contrast.matrix)
nrDEG<-re
sgnfc<-nrDEG[abs(nrDEG$logFC)>2 & nrDEG$adj.P.Val<0.01,]
write.csv(sgnfc,file = "3.DE_2_0.01TvN.csv")
