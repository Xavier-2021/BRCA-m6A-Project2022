################### m6A enzymes DE ###########################
m6AexprSet<-new_exprSet[rownames(Getit),]
ModelexprSet<-new_exprSet[Genes,]
group_list<-c(rep('T',ncol(tumor)),rep('N',ncol(normal)))
library(reshape2)
m6Aexp_L<-melt(as.matrix(m6AexprSet))
m6Aexp_L$group<-rep(group_list,each=20)
colnames(m6Aexp_L) = c('Gene','sample', 'expression','group')    

ggplot(m6Aexp_L,aes(x=Gene,y=expression,color=group))+geom_boxplot()+
  stat_compare_means(aes(label =..p.signif..))+theme_set(theme_bw())+
  coord_flip()+ggtitle('m6A Enzymes DE')

library(scatterplot3d)
df<-m6AexprSet
df<-t(df)
dfGroup<-as.data.frame(group_list)
rownames(dfGroup)<-colnames(m6AexprSet)
pca_result <- prcomp(df,scale=T)
pca_result$x<-data.frame(pca_result$x)
tb<-pca_result$x
colors <- c(rep("red",1102),rep("green",99))
pVar <- pca_result$sdev^2/sum(pca_result$sdev^2)
pVar = round(pVar,digits = 2)
s3d <- scatterplot3d(pca_result$x[,1:3],
                     pch = 16,       # 点形状
                     color=colors,   # 点颜色
                     cex.symbols = 2 # 点大小
)
legend("right",
       legend = unique(dfGroup[,1]),
       col =  c("red","green"),
       pch = 16,
       inset = -0.1,
       xpd = TRUE,
       horiz = TRUE)
ggtitle('m6A PCA')

#m6A enzymes survival 
library(survival)
library(survminer)
library(SummarizedExperiment)
ncol(Getit)
Get_L<-data.frame(t(Getit))
Get_L<-Get_L[rownames(clid_ps),]
data<-cbind(Get_L,clid_ps[,c('logic_vital','time')])

for (gene in rownames(Getit)) {
  data$exp<-ifelse(data[,gene]>median(data[,gene]),"H","L")
  sdf<-survdiff(Surv(time,logic_vital) ~ exp , data = data)
  p.val<-1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  if (p.val<0.05) {
    fit <- survfit(Surv(time,logic_vital) ~ exp , data = data)
    res.sum <- surv_summary(fit,data=data)
    pdf(paste0(gene,' survival.pdf'),width = 4,height = 6,onefile = F)
    print(ggsurvplot(fit,
                     pval = TRUE, conf.int = F,
                     risk.table = TRUE, # Add risk table
                     risk.table.col = "strata", # Change risk table color by groups
                     linetype = "strata", # Change line type by groups
                     surv.median.line = "hv", # Specify median survival
                     ggtheme = theme_bw(), # Change ggplot2 theme
                     palette = c("#E7B800", "#2E9FDF"),
                     data=data
    )+ggtitle((paste0(gene,' survival'))))
    dev.off()
  }
}

