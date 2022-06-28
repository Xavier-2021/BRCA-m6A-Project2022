############## download clinicaldata #################
library(cgdsr)
library(DT)
#create a CGDS object
mycgds <- CGDS("http://www.cbioportal.org/")
#Specify a dataset
mycancerstudy = 'brca_tcga' 
getCaseLists(mycgds,mycancerstudy)[,1]
#Get available genetic data profiles for a specific cancer study
getGeneticProfiles(mycgds,mycancerstudy)[,1]
mycaselist ='brca_tcga_rna_seq_v2_mrna'  
mygeneticprofile = 'brca_tcga_rna_seq_v2_mrna'  
#Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)
cliD<-DT::datatable(myclinicaldata,
                    extensions = 'FixedColumns',
                    options = list(                    #dom = 't',
                      scrollX = TRUE,
                      fixedColumns = TRUE
                    ))
write.csv(myclinicaldata,file='cliD_LIHC.csv')

############## bind clinicaldata and expression dataset #################

cliD<-myclinicaldata
colnames(cliD)
cliD<-cliD[,c("OS_MONTHS","OS_STATUS")]
clid_p<-cliD[which(substr(rownames(cliD),9,12)%in%substr(colnames(tumor),9,12)),]
poser<-tumor[rownames(tumor) %in% rownames(sgnfc),]#
colnames(poser)<-substr(colnames(poser),9,12)

clid_p<-clid_p[!duplicated(substr(rownames(clid_p),9,12)),]
rownames(clid_p)<-substr(rownames(clid_p),9,12)
poser<-poser[,rownames(clid_p)]#´
ncol(poser)
crucial_cols<-t(poser)
clid_p<-cbind(clid_p,crucial_cols)
write.csv(clid_p,file = "Clinical_DE.csv")

write.csv(tumor,file='tumor.csv')

####### Find Prognostic Significant Genes #######
library(survival)
library(survminer)
library(SummarizedExperiment)
results<-c()
vector<-colnames(clid_p[,c(3:ncol(clid_p))])
data<-clid_p
data$logic_vital<-ifelse(data$OS_STATUS=='1:DECEASED',2,1)
data$time<-as.numeric(data$OS_MONTHS)
for (gene in vector) {
  data$exp<-ifelse(data[,gene]>median(data[,gene]),"H","L")
  if ("H"%in%data$exp & "L"%in%data$exp) {
    sdf<-survdiff(Surv(time,logic_vital) ~ exp , data = data)
    p.val<-1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    if (p.val<0.01) {
      results<-c(results,gene)
    }
  }
}
df<-c()
for (gene in results) {
  data$exp<-ifelse(data[,gene]>median(data[,gene]),"H","L")
  if ("H"%in%data$exp & "L"%in%data$exp) {
    outsur6<-coxph(Surv(time,logic_vital)~ data[,gene],data = data)
    tb<-summary(outsur6)[["coefficients"]]
    tb<-data.frame(tb)
    p.val<-tb$Pr...z..
    if (p.val<0.01) {
      df<-c(df,gene)
    }
  }
}
df=as.data.frame(df)
df<-na.omit(df)
df$df<-gsub('[.]','-',df$df)
sig_surv<-df$df 

sig_surv<-sig_surv[sig_surv %in% LNCS$V1] 
LNC_kick<-tumor[which(rownames(tumor) %in% sig_surv),]
colnames(LNC_kick)<-substr(colnames(LNC_kick),9,12)

######### TvN volcano plot ###############
library(ggplot2)
logFC_cutoff=2
DEG<-re[rownames(re) %in% LNCS$V1,]
DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.01 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)

g = ggplot(data=DEG, 
           aes(x=logFC, y=-log10(adj.P.Val), 
               color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 adj.p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','#A9A9A9','red'))+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())
print(g)
ggsave(g,filename = 'volcano_TvN.pdf',width = 6,height = 6)
