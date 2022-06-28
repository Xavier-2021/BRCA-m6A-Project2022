Eraser<-c("ALKBH5","FTO")
Writer<-c("METTL3","METTL14","METTL16","WTAP","RBM15","RBM15B","VIRMA","CBLL1","ZC3H13")
Reader<-c("YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3",
          "IGF2BP1","IGF2BP2","IGF2BP3","PRRC2A")

Getit<-tumor[c(Eraser,Writer,Reader),]
Getit<-na.omit(Getit)
colnames(Getit)<-substr(colnames(Getit),9,12)

All<-rbind(Getit,LNC_kick)
### filt by significant correlation to m6A enzymes ###
library(ggcorrplot)
library(pheatmap)
library(SummarizedExperiment)
samples_Cor<-t(All)
data_Cor<-cor(samples_Cor,use = "everything",method = "pearson")
data_Cor<-data_Cor[sig_surv,rownames(Getit)]
data_Cor_pre<-data_Cor
muu<-t(All)
corr.p<-ggcorrplot::cor_pmat(muu)

data_Cor<-data_Cor[which(abs(rowMaxs(data_Cor))>0.3),]
corr.p<-corr.p[rownames(data_Cor_pre),rownames(Getit)]
data_Cor<-data_Cor[which(rownames(data_Cor) %in% rownames(corr.p)),]

col_anno<-data.frame(type=c(rep("E",2),rep("W",9),rep("R",9)),row.names = colnames(data_Cor))
ann_colors = list( Time = c("Blue", "Red", "Yellow"))
pheatmap(data_Cor_pre,cluster_cols = F,
         annotation_col = col_anno,
         annotation_colors = ann_colors,
         gaps_col = c(2,11),
         display_numbers = data.frame(ifelse(corr.p < 0.01, "**", ifelse(corr.p<0.05,"*",""))),
         fontsize_number = 10
)
write.csv(data_Cor,file="cor_m6A_matrix.csv")

################## Genes MultiCOX ######################
library(survival);library(survminer);library(SummarizedExperiment);library(forestplot)
data<-clid_ps
data$logic_vital<-ifelse(data$OS_STATUS=='1:DECEASED',2,1)
data$time<-as.numeric(data$OS_MONTHS)

outsur1 <-coxph(Surv(time,logic_vital)~ LINC01235+`P3H2-AS1`+LINC01198
                ,data = clid_ps)
tb1<-cbind(summary(outsur1)[["coefficients"]],summary(outsur1)[["conf.int"]])
tb1
write.csv(tb1,file="多因素univa.csv",quote = F)

tb1<-c()
data<-clid_p
data$logic_vital<-ifelse(data$OS_STATUS=='1:DECEASED',2,1)
data$time<-as.numeric(data$OS_MONTHS)
for (gene in sig_surv) {
  outsur2 <-coxph(as.formula(paste0('Surv(time,logic_vital)~ `',gene,'`'))
                  ,data = data)
  tb2<-cbind(summary(outsur2)[["coefficients"]],summary(outsur2)[["conf.int"]])
  tb1<-rbind(tb1,tb2)
}
tb1
write.csv(tb1,file="Univariate sig_surv.csv",quote = F) # Edit in Excel

tb<-read.csv("Univariate sig_surv.csv",header = F)
colnames(tb)<-tb[1,]
cochrane_from_rmeta<-rbind(c(NA,NA,NA),tb[-1,c(4:6)])
h<-data.frame(tb[,c(1:3)])
pdf('sigsurv_uniCOX.pdf',onefile = F)
forestplot(h,
           mean=cochrane_from_rmeta$mean,
           lower=cochrane_from_rmeta$lower,
           upper=cochrane_from_rmeta$upper,
           zero = 1, 
           clip = c(0,3),
           xticks = c(0.7,0.85,1,1.15,1.3,1.45),
           xlog=FALSE, 
           fn.ci_norm = fpDrawDiamondCI, 
           boxsize = 0.4, 
           col=fpColors(line = "#CC79A7", 
                        box="#D55E00"), 
           lty.ci = 7,   
           lwd.ci = 3,   
           ci.vertices.height = 0.05, 
           txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1), 
           lineheight = "auto", 
) 
dev.off()
gc()

####### glmnet: LASSO model build #######
clid_train<-clid_p[sample(466,300),]
clid_test<-clid_p[-which(rownames(clid_p) %in% rownames(clid_train)),]

# glmnet 
colnames(clid_p)<-gsub('[.]','-',colnames(clid_p))

write.csv(clid_p,file="clinical_for_glmnet_input.csv")

dat_use<-clid_train[,rownames(data_Cor)]
dat_model<-dat_use
dat_model$y<-data[match(rownames(dat_model), rownames(data)),"time"]
dat_model$z<-data[match(rownames(dat_model), rownames(data)),"logic_vital"]
dat_model<-na.omit(dat_model)
#set.seed(1) # 
N <- as.numeric(nrow(dat_use)) # split data
test_index <- sample(N,0.4*N)
train_index <- c(1:N)[-test_index]
test_data <- dat_model[test_index,]
train_data <- dat_model[train_index,]
# 
y_name <- c('y','z')
thedat <- na.omit(train_data)
thedat<-thedat[which(thedat$y>0),]
y <-  thedat[,y_name]
library(survival)
y$y<-as.double(y$y)
y$z<-as.double(y$z)
y <-data.matrix(Surv(y$y,y$z))
x <- thedat
x <- data.matrix(x)
x<-x[,-c(ncol(x),ncol(x)-1)]
#
library(glmnet)
# fit the model
fit <- glmnet(x, y,family = 'cox')
# area under the ROC curve, CV lambda
# set.seed
#set.seed(1)
fit_cv <- cv.glmnet(x, y,family='cox')
pdf("glmnet_ficv.pdf",width=5,height = 4)
plot(fit_cv,cex.axis=1,cex.lab=1)
dev.off()
pdf("glmnet_cofi.pdf",width=5,height = 4)
plot(fit,cex.axis=1,cex.lab=1,lwd=2)
dev.off()
Coefficients<-as.matrix(coef(fit_cv,s="lambda.min"))
fit_cv$lambda.min
rownames(Coefficients) <- gsub('[.]','-',rownames(Coefficients))
write.csv(Coefficients,file = "Coefficients.csv")
Coefficients<-read.csv("Coefficients.csv")
Genes<-Coefficients[which(Coefficients$X1!=0),"X"]

####### Calculate Riskscore of each patient #######
rownames(Coefficients)<-Coefficients$X
Coefficients<-Coefficients[Genes,]
riskscore<-rownames(clid_p)
clid_ps<-cbind(clid_p[,c(1:2)],clid_p[,Genes]) # 
for (gene in Genes) {
  r<-clid_ps[,gene]*Coefficients[gene,"X1"]
  riskscore<-cbind(riskscore,as.numeric(r))
}
rownames(riskscore)<-riskscore[,1]
riskscore<-riskscore[,-1]
write.csv(riskscore,file = "riskscore.csv") # Use Excel to edit

riskscore<-read.csv("riskscore.csv",header=F)
clid_ps$Riskscore<-riskscore$V6
write.csv(clid_ps,file = "clinical_RS.csv")

##### plot Riskscore in training and testing datasets ######
library(ggplot2)
mayo=clid_ps[rownames(clid_test),]
rownames(mayo)<-c(1:nrow(mayo))
mayo$num<-c(1:nrow(mayo))
mayo$num<-as.numeric(mayo$num)
mayo<-mayo[which(mayo$time!='NA'),]
p1<-ggplot(data=mayo,
           aes(x=num,y=time,color=OS_STATUS))+geom_point(alpha=0.4, size=1.75
           )+scale_colour_manual(values = c('green','red','grey')
           )+geom_vline(xintercept =0.5*nrow(mayo),lty=4,lwd=1)+labs(x=NULL,y="survival time")+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())+ggtitle('Testing Data')
may<-mayo[order(mayo$Riskscore),]
may$num<-c(1:nrow(may))
may$Risk_level<-ifelse(may$Riskscore>median(may$Riskscore),"High","Low")
p2<-ggplot(data=may,
           aes(x=num,y=Riskscore,color=`Risk_level`))+geom_point(alpha=0.4, size=1.75
           )+scale_colour_manual(values = c('red','green')
           )+geom_vline(xintercept =0.5*nrow(mayo),lty=4,lwd=1)+labs(x=NULL)+
  theme(panel.grid.major=element_line(colour=NA))+theme_set(theme_bw())

plots = list(p1,p2)
library(gridExtra)
library(cowplot)
pdf('distribution_RSplot_test.pdf')
plot_grid(p1, p2, nrow = 2, align = "v") #
dev.off()

# Model Genes DE
ModelexprSet<-new_exprSet[Genes,]
group_list<-c(rep('T',ncol(tumor)),rep('N',ncol(normal)))
library(reshape2)
Modelexp_L<-melt(as.matrix(ModelexprSet))
Modelexp_L$group<-rep(group_list,each=3)
colnames(Modelexp_L) = c('Gene','sample', 'expression','group')                      
ggplot(Modelexp_L,aes(x=Gene,y=expression,color=group))+geom_boxplot()+
  stat_compare_means(method = 'wilcox')+theme_set(theme_bw())+
  ggtitle('Risk lncRNAs DE')

### Riskscore survival KMplot ###
library(survival)
library(survminer)
library(SummarizedExperiment)
clid_ps$logic_vital<-ifelse(clid_ps$OS_STATUS=='1:DECEASED',2,1)
clid_ps$time<-as.numeric(clid_ps$OS_MONTHS)
data<-clid_ps[rownames(clid_test),]
data<-na.omit(data)

data$exp<-ifelse(data$Riskscore>median(data$Riskscore),"H","L")
fit <- survfit(Surv(time,logic_vital) ~ exp , data = data)
res.sum <- surv_summary(fit,data=data)
pdf('Risk score survival test.pdf',width = 4,height = 6,onefile = F)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           data=data
)+ggtitle("Risk Score Survival test")
dev.off()

### univariable Cox ###
tb1<-c()
data<-clid_ps
for (gene in c(Genes,'Riskscore')) {
  outsur2 <-coxph(as.formula(paste0('Surv(time,logic_vital)~ `',gene,'`'))
                  ,data = data)
  tb2<-cbind(summary(outsur2)[["coefficients"]],summary(outsur2)[["conf.int"]])
  tb1<-rbind(tb1,tb2)
}
write.csv(tb1,file="univaCOX.csv",quote = F) # edit using Excel 
tb<-read.csv("COXuniva.csv",header = F)
colnames(tb)<-tb[1,]
cochrane_from_rmeta<-tb[-1,c(4:6)]
h<-data.frame(tb[,c(1:3)])
pdf('uni_COX.pdf')
forestplot(h,  #
           mean=cochrane_from_rmeta$V4,
           lower=cochrane_from_rmeta$V5,
           upper=cochrane_from_rmeta$V6,
           zero = 1, #
           clip = c(0,3),
           xticks = c(0.5,0.75,1,1.25,1.5,1.75,2,2.5,3),
           xlog=FALSE, #
           fn.ci_norm = fpDrawDiamondCI, 
           boxsize = 0.1, 
           col=fpColors(line = "#CC79A7", #
                        box="#D55E00"), #
           lty.ci = 7,   # 
           lwd.ci = 3,   # 
           ci.vertices.height = 0.05, # 
           txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1), #
           lineheight = "auto", #
) 
dev.off()

######### Risk score ROC analysis at 1\2\3\5 years###
mayo<-clid_ps[rownames(clid_ps),]
daysROCprepare<-function(days){
  mayo$year<-ifelse(mayo$logic_vital==2,ifelse(mayo$time<days+1,0,1),ifelse(mayo$time>=days,1,NA))
  mame<-mayo[which(mayo$year!='NA'),]
  return(mame)
}
library(pROC)
library(ggplot2)
mame<-daysROCprepare(12)
ROC1 <- roc(response = mame$year, predictor=mame$Riskscore,levels = c(0,1))
mame<-daysROCprepare(24)
ROC2 <- roc(response = mame$year, predictor=mame$Riskscore,levels = c(0,1))
mame<-daysROCprepare(36)
ROC3 <- roc(response = mame$year, predictor=mame$Riskscore,levels = c(0,1))
mame<-daysROCprepare(60)
ROC4 <- roc(response = mame$year, predictor=mame$Riskscore,levels = c(0,1))
pdf('ROC.pdf')
plot(ROC1,method='density', 
     col = 'blue', 
     legacy.axes = TRUE, 
     xlab = '1-Specificity',
     print.auc =TRUE,
     print.auc.x = 0.4,print.auc.y = 0.5,
     auc.polygon = F
)
plot.roc(ROC2,
         add=T,  
         col="red", 
         print.auc=TRUE, print.auc.x=0.4,print.auc.y=0.4,
         smooth = F)  
plot.roc(ROC3,
         add=T,  
         col="black",  
         print.auc=TRUE, print.auc.x=0.4,print.auc.y=0.3,
         smooth = F) 
plot.roc(ROC4,
         add=T, 
         col="brown", 
         print.auc=TRUE, print.auc.x=0.4,print.auc.y=0.2,
         smooth = F) 
legend(1,1,  
       bty = "n",  #
       title="",   
       legend=c("1year","2year","3year","5year"),  
       col=c("blue","red",'black','brown'),  
       lwd=2)  
dev.off()


####### nomogram #######
cliD<-myclinicaldata
colnames(cliD)
cliD<-cliD[!duplicated(substr(rownames(cliD),9,12)),]
rownames(cliD)<-substr(rownames(cliD),9,12)
cliD<-cliD[rownames(clid_ps),]
clid_t<-clid_ps
clid_t$tumor_stage<-cliD$AJCC_PATHOLOGIC_TUMOR_STAGE
clid_t$gender<-cliD$SEX
clid_t$age<-cliD$AGE
clid_t$metastasis_pm<-cliD$AJCC_METASTASIS_PATHOLOGIC_PM
clid_t$weight<-cliD$SAMPLE_INITIAL_WEIGHT

library(Hmisc); library(grid); library(lattice);library(Formula); library(ggplot2)
library(rms)
data=clid_t
data<-data[which(data$time>0),]
data<-data[-which(data$tumor_stage==''),]
dd=datadist(data)
options(datadist="dd")#
#data<-data[-which(data$metastasis_pm==''),]
f2 <- psm(Surv(time,logic_vital) ~ age+tumor_stage+gender+Riskscore+metastasis_pm+weight, data =  data, dist='lognormal')
med <- Quantile(f2)  
surv <- Survival(f2) 
nom <- nomogram(f2, fun=function(x) med(lp=x),
                funlabel="Median Survival Time")
plot(nom)
nom <- nomogram(f2, fun=list(function(x) surv(365, x),
                             function(x) surv(730, x),
                             function(x) surv(1095,x),
                             function(x) surv(1825,x)),
                funlabel=c("1-year Survival Probability",
                           "2-year Survival Probability",
                           "3-year Survival Probability",
                           "5-year Survival Probability"))
pdf('nomogram.pdf',width=16)
plot(nom, xfrac=.2)
dev.off()
### Calibrate ###
rcorrcens(Surv(time,logic_vital) ~ predict(f2), data =  data)
f2 <- psm(Surv(time,logic_vital) ~ age+tumor_stage+gender+Riskscore+metastasis_pm+weight, 
          data =  data, x=T, y=T, dist='lognormal')
cal5 <- calibrate(f2, cmethod='KM', method="boot", u=60, m=300, B=1089)
cal3 <- calibrate(f2, cmethod='KM', method="boot", u=36, m=300, B=1089)
cal2 <- calibrate(f2, cmethod='KM', method="boot", u=24, m=300, B=1089)
cal1 <- calibrate(f2, cmethod='KM', method="boot", u=12, m=300, B=1089)
pdf("Calibrates1235.pdf")
plot(cal5,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.4,1),ylim=c(0.4,1),xlab = 'predicted',ylab = 'fraction surviving',
     col='brown')
par(new=T)
plot(cal3,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.4,1),ylim=c(0.4,1),xlab = 'predicted',ylab = 'fraction surviving',
     col='black')
par(new=T)
plot(cal2,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.4,1),ylim=c(0.4,1),xlab = 'predicted',ylab = 'fraction surviving',
     col='red')
par(new=T)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.4,1),ylim=c(0.4,1),xlab = 'predicted',ylab = 'fraction surviving',
     col='blue')
legend(0,1,  # 图例位置
       bty = "n",  # 图例样式，默认为 "o"
       title="",   # 引号内添加图例标题
       legend=c("1year","2year","3year","5year"),  # 添加分组
       col=c("blue","red",'black','brown'),  # 颜色跟前面一致
       lwd=2)  # 线条粗细
dev.off()

