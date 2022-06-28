############### Riskscore in phenotypes ################
rownames(clid_ps)==rownames(cliD)
clid_ps$Tumor_Stage<-cliD$AJCC_PATHOLOGIC_TUMOR_STAGE
table(clid_ps$Tumor_Stage)
clid_ps$Tumor_Stage<-ifelse(clid_ps$Tumor_Stage %in% c('Stage I','Stage IA','Stage IB'),
                            'Stage I',clid_ps$Tumor_Stage)
clid_ps$Tumor_Stage<-ifelse(clid_ps$Tumor_Stage %in% c('Stage II','Stage IIA','Stage IIB','Stage IIC'),
                            'Stage II',clid_ps$Tumor_Stage)
clid_ps$Tumor_Stage<-ifelse(clid_ps$Tumor_Stage %in% c('Stage III','Stage IIIA','Stage IIIB','Stage IIIC'),
                            'Stage III',clid_ps$Tumor_Stage)
table(clid_ps$Tumor_Stage)

clid_ps$Pathologic_PT<-cliD$AJCC_TUMOR_PATHOLOGIC_PT
table(clid_ps$Pathologic_PT)
clid_ps$Pathologic_PT[grep('T1',clid_ps$Pathologic_PT)]<-'T1'
clid_ps$Pathologic_PT[grep('T2',clid_ps$Pathologic_PT)]<-'T2'
clid_ps$Pathologic_PT[grep('T3',clid_ps$Pathologic_PT)]<-'T3'
clid_ps$Pathologic_PT[grep('T4',clid_ps$Pathologic_PT)]<-'T4'
table(clid_ps$Pathologic_PT)

clid_ps$age<-as.numeric(cliD$AGE)
clid_ps$Gender<-cliD$SEX
table(clid_ps$Gender)
clid_ps$age
clid_ps$Age_subtype<-ifelse(clid_ps$age>65,'Elder','Younger')

library(ggplot2);library(ggpubr)
data<-clid_ps[which(clid_ps$Tumor_Stage!='' & clid_ps$Tumor_Stage!='Stage X'),]
data<-na.omit(clid_ps)
ggplot(data,aes(x=Tumor_Stage,y=Riskscore,color=Tumor_Stage))+geom_boxplot()+geom_jitter(width = 0.2)+
  stat_compare_means()+theme_set(theme_bw())

data<-clid_ps[clid_ps$Age_subtype=='Younger',]
data<-rbind(clid_ps[clid_ps$Age_subtype=='Elder',],data[sample(700,318),])
data<-na.omit(data)
ggplot(data,aes(x=Age_subtype,y=Riskscore,color=Age_subtype))+geom_boxplot()+geom_jitter(width = 0.2)+
  stat_compare_means()+theme_set(theme_bw())

g<-ggplot(data,aes(x=Riskscore,y=age))+stat_smooth(method="lm")+
  geom_point(cex=1)+geom_smooth(method=lm,se=F)+stat_cor(method = "pearson", label.y =0)+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())+xlab('Risk score')
library(ggExtra)
ggMarginal(g,fill='darkgreen')