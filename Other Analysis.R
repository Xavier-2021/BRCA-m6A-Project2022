
################### Genes of interests ##################

# correlation
data<-data.frame(t(tumor))
data2<-data[,c('LINC01198','ZEB1','MYC')]
g1<-ggplot(data,aes(x=ZEB1,y=LINC01198))+stat_smooth(method="lm")+
  geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson")+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())
g2<-ggplot(data,aes(x=MYC,y=LINC01198))+stat_smooth(method="lm")+
  geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson")+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())
ggsave(g2,filename = 'LINC01198 MYC cor.pdf',width = 4,height = 4)
ggarrange(g1,g2,ncol=2)

# Tumor vs Normal Different Expresssion
LGexprSet<-new_exprSet[c(Genes,'IGF2BP2'),]
group_list<-c(rep('T',ncol(tumor)),rep('N',ncol(normal)))
library(reshape2)
LGexprSet_L<-melt(as.matrix(LGexprSet))
LGexprSet_L$group<-rep(group_list,each=4)
colnames(LGexprSet_L) = c('Gene','sample', 'expression','group')                      
LGexprSet_L<-LGexprSet_L[LGexprSet_L$Gene=='IGF2BP2',]

ggplot(LGexprSet_L,aes(x=Gene,y=expression,color=group))+geom_boxplot()+
  stat_compare_means(method = 'wilcox')+theme_set(theme_bw())+
  ggtitle('IGF2BP2 TCGA-BRCA')


# survival
library(survival)
library(survminer)
library(SummarizedExperiment)
data<-data.frame(t(tumor))
rownames(data)
data1<-clid_ps
data1$IGF2BP2<-data[rownames(data1),'IGF2BP2']
data1$IGF2BP2
data1$exp<-ifelse(data1$IGF2BP2>median(data1$IGF2BP2),"H","L")
fit <- survfit(Surv(time,logic_vital) ~ exp , data = data1)
res.sum <- surv_summary(fit,data=data)
pdf('IGF2BP2 survival.pdf',width = 4,height = 5,onefile = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           data=data1
)+ggtitle("IGF2BP2 Survival")
dev.off()

# relavent genes enrichment


High2Low<-function(Gene){
  gexp<-data.frame(t(tumor[Gene,]))
  gexp$modType<-ifelse(gexp[,1]>median(gexp[,1]),'H','L')
  tumor<-tumor[,rownames(gexp)]
  
  modType<-gexp$modType
  library(limma)
  design <- model.matrix(~0+factor(modType))
  colnames(design)=levels(factor(modType))
  rownames(design)=colnames(tumor)
  design
  contrast.matrix<-makeContrasts('H-L', levels = design)
  contrast.matrix 
  deg = function(exprSet,design,contrast.matrix){
    fit <- lmFit(exprSet,design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2, coef=1, n=Inf)
    return(tempOutput)
  }
  reGene = deg(tumor,design,contrast.matrix)
  write.csv(reGene,file = paste0('Highlow_DE_',Gene,'.csv'))
  return(reGene)
  
  ##功能富集分析##
  ##将差异基因列表转化为entrezid作为后续的富集分析#
  DEup<-rownames(reGene[which(reGene$logFC>0 & reGene$adj.P.Val<0.01),])
  DEdown<-rownames(reGene[which(reGene$logFC<0 & reGene$adj.P.Val<0.01),])
  library(clusterProfiler)
  library(org.Hs.eg.db)
  DEup<-bitr(DEup,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db',drop = TRUE)
  DEdown<-bitr(DEdown,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db',drop = TRUE)
  ##GO富集####
  e_up_go <- enrichGO(gene          = DEup$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
  p1<-dotplot(e_up_go)+ggtitle('upregulated GO')
  e_down_go <- enrichGO(gene          = DEdown$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)
  p2<-dotplot(e_down_go,showCategory=10)+ggtitle('downregulated GO')
  library(gridExtra)
  library(cowplot)
  pdf(paste0('GO_',Gene,'.pdf'),width = 8)
  print(plot_grid(p1, p2, nrow = 2, align = "v"))
  dev.off()
  
  ##KEGG富集####
  eKoup <- enrichKEGG(
    gene          = DEup$ENTREZID,
    keyType     = "kegg",
    organism   = 'hsa',
    pvalueCutoff      = 0.05,
    pAdjustMethod     = "BH",
    qvalueCutoff  = 0.1
  )
  eKodown <- enrichKEGG(
    gene          = DEdown$ENTREZID,
    keyType     = "kegg",
    organism   = 'hsa',
    pvalueCutoff      = 0.05,
    pAdjustMethod     = "BH",
    qvalueCutoff  = 0.1
  )
  k1<-dotplot(eKoup,showCategory=10)+ggtitle('upregulated KEGG')
  k2<-dotplot(eKodown,showCategory=10)+ggtitle('downregulated KEGG')
  library(gridExtra)
  library(cowplot)
  pdf(paste0('KEGG_',Gene,'.pdf'))
  print(plot_grid(k1, k2, nrow = 2, align = "v"))
  dev.off()
  
}

for (genename in c('LINC01198','IGF2BP2','')) {
  reGene<-High2Low(genename)
}


