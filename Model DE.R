####### Different Expression analysis in Riskscore high/low groups ##########
data=clid_ps
data$exp<-ifelse(data$Riskscore>median(data$Riskscore),"H","L")
H_list<-rownames(data[which(data$exp=='H'),])
L_list<-rownames(data[which(data$exp=='L'),])
colnames(tumor)<-substr(colnames(tumor),9,12)

exp_H<-tumor[,which(colnames(tumor)%in%H_list)]
exp_L<-tumor[,which(colnames(tumor)%in%L_list)]
RS_exp<-cbind(exp_H,exp_L)
Description<-c(rep("NONE",nrow(RS_exp)))
RS_frame<-cbind(Description,RS_exp)
write.table(RS_frame,file = "RS_exp.txt",quote = F,sep = '\t') # then edit into GSEA software input format(gct)
modType=c(rep("High",ncol(exp_H)),rep("Low",ncol(exp_L))) # differential
write.table(t(data.frame(modType)),file = "modType.txt")

library(limma) # limma
design <- model.matrix(~0+factor(modType))
colnames(design)=levels(factor(modType))
rownames(design)=colnames(RS_exp)
design
contrast.matrix<-makeContrasts("High-Low", levels = design) 
contrast.matrix 
deg = function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  fit2 <- eBayes(fit2)  ## default no trend !!!
  tempOutput = topTable(fit2, coef=1, n=Inf)
  return(tempOutput)
}
re = deg(RS_exp,design,contrast.matrix)
write.csv(re,file = "RiskScore_related_DEGs.csv")

### clusterProfiler: KEGG/GO enrichment of DEs ###
DE<-DEG
genelist_down<-DE[which(DE$logFC<(-1.5) & DE$adj.P.Val<0.05),]
genelist_up<-DE[which(DE$logFC>0 & DE$adj.P.Val<0.05),]
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)
genelist_up_l<-bitr(rownames(genelist_up),fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db",drop=TRUE)
eKo <- enrichKEGG(
  gene          = genelist_up_l$ENTREZID,
  keyType     = "kegg",
  organism   = 'human',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.1
)
eGo <- enrichGO(
  gene          = genelist_up_l$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.1,
  readable      = TRUE)
p1<-barplot(eGo)+ggtitle('up-regulated GO')
print(p1)
#ggsave(g,filename = "riskscore_upGO.pdf",width = 14,height = 8)
p2<-barplot(eKo)+ ggtitle('up-regulated KEGG')
#ggsave(g,filename = "riskscore_upKEGG.pdf",width = 10,height = 8)
library(gridExtra)
library(cowplot)
pdf('enrichment_up.pdf',height = 8)
plot_grid(p1, p2, nrow = 2, align = "v")
dev.off()
