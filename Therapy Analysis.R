##################### TIDE ##########################
TIDE<-read.table('11.TIDE/TCGA.BRCA.RNASeq.norm_subtract.OS_base',sep = '\t',header = T)
rownames(TIDE)<-substr(rownames(TIDE),9,12)
TIDE$Riskscore<-clid_ps[rownames(TIDE),'Riskscore']
TIDE<-na.omit(TIDE)
TIDE$Riskscore_exp<-ifelse(TIDE$Riskscore>median(TIDE$Riskscore),'High','Low')
library(ggplot2);library(ggpubr)
ggplot(TIDE,aes(x=Riskscore_exp,y=M2,color=Riskscore_exp))+geom_violin()+geom_jitter(width = 0.2)+
  stat_compare_means()+ggtitle('TIDE in Riskscore High/Low group')

TIDE.prog<-read.csv('11.TIDE/TIDE output.csv',row.names = 1)
TIDE.prog<-TIDE.prog[rownames(clid_ps),]
TIDE.prog$Riskscore<-clid_ps[rownames(TIDE.prog),'Riskscore']
TIDE.prog<-na.omit(TIDE.prog)
TIDE.prog$Riskscore_exp<-ifelse(TIDE.prog$Riskscore>median(TIDE.prog$Riskscore),'High','Low')
library(ggplot2);library(ggpubr)
ggplot(TIDE.prog,aes(x=Riskscore_exp,y=TIDE,color=Riskscore_exp))+geom_violin()+geom_jitter(width = 0.2)+
  stat_compare_means()+ggtitle('TIDE in Riskscore High/Low group')

#################### immuno Check Point #############
ICPgenes<-c('SNCA','CD274','CTLA4','LAG3','HAVCR2','TJP1','IDO1')
data<-tumor[ICPgenes,]
data<-data.frame(t(data))
data<-data[rownames(clid_ps),]
data$Riskscore<-clid_ps$Riskscore

library(corrplot)
library(ggplot2)

data<-na.omit(data)
ccc<-cor(data,use = "everything",method = "pearson")
ccc
cccor<-ggcorrplot::cor_pmat(data)
pdf('corhallmark.pdf',height = 8,width = 8)
corrplot(ccc,type = 'upper',p.mat = cccor,insig = 'label_sig',pch.cex = 1,
         sig.level = c(.001, .01, .05),
         method = 'circle',
         tl.cex = 1,tl.col = 'black')+
  theme(axis.text.y = element_text(angle = 0),axis.title.x = element_text(angle = 0))
dev.off()

library(ggpubr)
g<-ggplot(data,aes(x=Riskscore,y=HAVCR2))+stat_smooth(method="lm")+
  geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson", label.y =0)+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())+xlab('Risk score')
library(ggExtra)
ggMarginal(g,fill='darkgreen')

##################### pRRophetic to predict drug sensitivity ##############
################# Run in another dir ############
write.csv(clid_ps,'clid_ps.csv') 

######## clid_ps.csv and tumor.csv was copy to working dictionary ###############
library(parallel)
library(pRRophetic)
library(ggplot2)
library(ggpubr)

exprData<-tumor
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
data(cgp2016ExprRma)
possibleDrugs2016 <- unique( drugData2016$Drug.name)
possibleDrugs2016

clid_ps<-read.csv('clid_ps.csv',row.names = 1)
tumor<-read.csv('tumor.csv',row.names = 1)

tumor<-tumor[,!duplicated(substr(colnames(tumor),9,12))]
exprData<-tumor
colnames(exprData)<-substr(colnames(exprData),9,12)
r<-intersect(colnames(exprData),rownames(clid_ps))
clid_ps<-clid_ps[r,]
exprData<-exprData[,r]
Data<-rbind(exprData,clid_ps$Riskscore)
Riskscore<-clid_ps$Riskscore

drugDraf<-c()
possibleDrugs2016
possibleDrugs2016 <- unique( drugData2016$Drug.name)
Druglist<-c()

DrugSig<-function(exprData,DDrug){
  predictedPtype<-pRRopheticPredict(
    testMatrix=as.matrix(exprData), # can be my exprData
    drug=DDrug,
    tissueType = "all", 
    batchCorrect = "eb",
    selection=1,
    dataset = "cgp2016")
  df <- data.frame(predictedPtype)
  return(df)
}

for (drug in c('681640','Erlotinib','Navitoclax')) {
  df<-tryCatch({
    DrugSig(exprData,drug)
  },error = function(e){
    3
    print("error")
  })
  if (df!=3) {
    try({
      df$Riskscore<-clid_ps[rownames(df),"Riskscore"]
      df<-na.omit(df)
      df$exp<-ifelse(df$Riskscore>median(df$Riskscore),'H','L')
      q<-wilcox.test(df[df$exp=='H',"predictedPtype"],df[df$exp=='L',"predictedPtype"])
      if (q$p.value<0.01) {
        pdf(paste0('Riskscore_',drug,'_predict_plot.pdf'),width = 5,height = 5)
        print(ggplot(df,aes(x=exp,y=predictedPtype,color=exp))+
                geom_boxplot(width=0.5)+geom_jitter(width=0.2)+
                xlab(NULL)+ggtitle(paste0(drug,'\nTCGA-BRCA Riskscore'))+
                stat_compare_means()+
                theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw()))
        dev.off()
        drugDraf<-rbind(drugDraf,
                        c(drug,q$p.value))
      }})
  }
}

write.csv(drugDraf,file = 'DrugsP.csv')

drugDraf<-rownames(clid_ps)
for (drug in c('681640','Erlotinib','Navitoclax')) {
  df<-tryCatch({
    DrugSig(exprData,drug)
  },error = function(e){
    3
    print("error")
  })
  if (df!='error') {
    colnames(df)<-drug
    drugDraf<-cbind(drugDraf,
                    df)
  }
  gc()
}

clid_compli<-cbind(clid_ps,drugDraf[,2:4])

########## Drug predictedPtype Clinical ##############
# KMplot #

library(survival);library(survminer);library(SummarizedExperiment)

data<-clid_compli
data<-na.omit(data)
data$exp<-ifelse(data$`681640`>median(data$`681640`),"H","L")
fit <- survfit(Surv(time,logic_vital) ~ exp , data = data)
res.sum <- surv_summary(fit,data=data)
pdf('681640 survival.pdf',width = 4,height = 5)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           data=data
)+ggtitle("681640")
dev.off()

table(data$Tumor_Stage)
data<-data[-which(data$Tumor_Stage %in% c('','Stage X')),]
library(ggplot2);library(ggpubr)

pdf('Erlotinib-Gender.pdf',width = 5,height = 3.5,onefile = F)
ggplot(data,aes(x=Gender,y=`Erlotinib`,color=Gender))+geom_boxplot()+
  stat_compare_means()+theme_set(theme_bw())+
  ggtitle('Erlotinib')
dev.off()




############## CIBERSORT: relationship of TME and Drugs ###################

CIBERSORT<-read.table('CIBERSORT_BRCA.txt',sep = '\t',header = T)
CIBERSORT<-CIBERSORT[!duplicated(CIBERSORT$SampleID),]
CIBERSORT<-CIBERSORT[-grep('11A',CIBERSORT$SampleID),]
CIBERSORT<-CIBERSORT[!duplicated(substr(CIBERSORT$SampleID,9,12)),]

rownames(CIBERSORT)<-substr(CIBERSORT$SampleID,9,12)
CIBERSORT<-CIBERSORT[rownames(clid_compli),]
clid_compli$Macrophage.M2<-CIBERSORT$Macrophages.M2
clid_compli$CD8T<-CIBERSORT$T.cells.CD8
clid_compli$NK.activated<-CIBERSORT$NK.cells.activated

library(ggplot2);library(ggpubr) ;   library(ggExtra)
for (cell in c('Macrophage.M2','CD8T','NK.activated')) {
  for (drug in c('681640','Erlotinib','Navitoclax')) {
    g<-ggplot(clid_compli,aes_string(x=paste0('`',cell,'`'),y=paste0('`',drug,'`')))+stat_smooth(method="lm")+
      geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson")+
      theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())
    pdf(paste0(cell,'_',drug,'_cor.pdf'),width = 5,height = 5)
    print(ggMarginal(g,fill='darkgreen'))
    dev.off()
  }
}

library(survival)
library(survminer)
library(SummarizedExperiment)

data<-clid_compli
data<-na.omit(data)
for (cell in c('Macrophage.M2','CD8T','NK.activated')){
  data$exp<-ifelse(data[,cell]>median(data[,cell]),"H","L")
  fit <- survfit(Surv(time,logic_vital) ~ exp , data = data)
  res.sum <- surv_summary(fit,data=data)
  pdf(paste0(cell,' survival.pdf'),width = 4,height = 5,onefile = F)
  show(ggsurvplot(fit,
                  pval = TRUE, conf.int = F,
                  risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "strata", # Change line type by groups
                  surv.median.line = "hv", # Specify median survival
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("#E7B800", "#2E9FDF"),
                  data=data
  )+ggtitle(cell))
  dev.off()
  
}
for (cell in c('Macrophage.M2','CD8T','NK.activated')){
  for (drug in c('681640','Erlotinib','Navitoclax')) {
    data$exp<-ifelse(data[,cell]>median(data[,cell]),"H","L")
    data$sense<-ifelse(data[,drug]>median(data[,drug]),'H','L')
    data$group<-ifelse(data$exp=='H',ifelse(data$sense=='H','High Cell High Sensitive','High Cell Low Sensitive'),
                       ifelse(data$sense=='H','Low Cell High Sensitive','Low Cell Low Sensitive'))
    fit <- survfit(Surv(time,logic_vital) ~ group , data = data)
    res.sum <- surv_summary(fit,data=data)
    pdf(paste0(cell,' ',drug,' survival.pdf'),width = 8.5,height = 6,onefile = F)
    show(ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "strata", # Change line type by groups
                    surv.median.line = "hv", # Specify median survival
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("purple1", "orange2",'violet','slateblue1'),
                    data=data,
                    legend='right',
    )+ggtitle(paste0(cell,' ',drug))
    )
    dev.off()
  }
  
}
