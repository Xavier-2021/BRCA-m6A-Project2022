############ TME analysis using ESTIMATE algorithm #############
library(estimate)

Esti<-read.csv('BRCA_estimate.csv',row.names = 1)
Esti<-Esti[!duplicated(substr(rownames(Esti),9,12)),]
rownames(Esti)<-substr(rownames(Esti),9,12)
Esti<-Esti[rownames(clid_ps),]
Esti_clid<-cbind(clid_ps,Esti)

library(survival)
library(survminer)
library(SummarizedExperiment)
data<-Esti_clid
data<-na.omit(data)

data$exp<-ifelse(data$ImmuneScore>median(data$ImmuneScore),"H","L") # check Pvalue one by one
fit <- survfit(Surv(time,logic_vital) ~ exp , data = data)
res.sum <- surv_summary(fit,data=data)
pdf('ImmuneScore survival.pdf',width = 4,height = 6,onefile = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           data=data
)+ggtitle("ImmuneScore Survival")
dev.off()
library(ggplot2);library(ggpubr)
ggplot(data,aes(x=Riskscore,y=ImmuneScore))+stat_smooth(method="lm")+
  geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson", label.y =0)+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())+xlab('Risk score')

################### CIBERSORT: TME analysis ####################
MC.col<-read.table('CIBERSORT_SKCM.txt',header = T)

MC.col<-MC.col[!duplicated(substr(MC.col$SampleID,9,12)),]
rownames(MC.col)<-substr(MC.col$SampleID,9,12)

MC.col<-MC.col[rownames(clid_ps),]
MC.col$Riskscore<-clid_ps$Riskscore
library(ggExtra)
for (Cell in colnames(MC.col)[3:24]) {
  tb<-cor.test(MC.col[,Cell],MC.col$Riskscore) 
  if (abs(tb$estimate)>0.3) {
    g<-ggplot(MC.col,aes_string(x='Riskscore',y=Cell))+stat_smooth(method="lm")+
      geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson")+
      theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())+xlab('Risk score')
    pdf(paste0(Cell,'_RS_cor.pdf'))
    print(ggMarginal(g,fill = 'darkgreen'))
    dev.off()
  }
}

data<-MC.col[,c(3:24,28:31)]
data<-na.omit(data)
ccc<-cor(data,use = "everything",method = "pearson")
ccc
cccor<-ggcorrplot::cor_pmat(data)
pdf('cor_CIBERSORT.pdf',height = 8,width = 8)
corrplot(ccc,type = 'lower',p.mat = cccor,insig = 'label_sig',pch.cex = 0.4,
         sig.level = c(.001, .01, .05),
         method = 'shade',
         tl.cex = 0.6,tl.col = 'black')+
  theme(axis.text.y = element_text(angle = 0),axis.title.x = element_text(angle = 0))
dev.off()

# Immune Cell Markers
for (Gene in Genes) {
  for (Marker in c('CD8A','MRC1','TNF')) {
    data<-data.frame(t(tumor[c(Gene,Marker),]),check.names = F)
    g<-ggplot(data,aes_string(x=paste0('`',Gene,'`'),y=Marker))+stat_smooth(method="lm")+
      geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson")+
      theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())+
      xlab(Gene)+ylab(Marker)
    pdf(paste0(Gene,' ',Marker,' Cor.pdf'),width = 5,height = 5,onefile = F)
    print(g)
    dev.off()  
  }
}
