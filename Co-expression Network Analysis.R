###################### WGCNA ################################

library(WGCNA)
library(data.table)
library(stringr)
library(reshape2)
library(stringr)
library(AnnotationDbi)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=4
memory.limit(size = 20000)
corType='pearson'
robustY = ifelse(corType=="pearson",T,F)

dataExpr<-tumor
dataExpr<-na.omit(dataExpr)

m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))

traitData<-clid_t

exprSize = checkSets(fixDataStructure(dataExpr))
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;
gsg = goodSamplesGenesMS(fixDataStructure(dataExpr), verbose = 3)
gsg$allOK

dim(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
clust = cutreeStatic(sampleTree, cutHeight = 600, minSize = 10)
keepSamples = (clust== 1)
dataExpr = dataExpr[keepSamples, ]
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")


powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType='unsigned', verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red")
abline(h=0.9,col="red")
# Soft thresholdä
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
#power = sft$powerEstimate
power=sft$powerEstimate

##One-step network construction and module detection

cor <- WGCNA::cor
net = blockwiseModules(dataExpr, power = power, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=T, corType = "pearson", 
                       maxPOutliers=1, loadTOMs=TRUE,
                       saveTOMFileBase = paste0("1.tom"),
                       verbose = 3)
#cor<-stats::cor

table(net$colors)
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = F, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

merge_modules = mergeCloseModules(dataExpr, moduleColors, cutHeight = 0.1, verbose = 3)
mergedColors = merge_modules$colors
mergedMEs = merge_modules$newMEs
plotDendroAndColors(net$dendrograms[[1]], cbind(moduleColors[net$blockGenes[[1]]], mergedColors[net$blockGenes[[1]]]),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# module eigengene
MEs = net$MEs
MEs_col = mergedMEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# marDendro/marHeatmap
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
traitData<-traitData[rownames(MEs_col),]
Riskscore<-traitData$Riskscore

expcli<-cbind(MEs_col, Riskscore)
MEs_colpheno = orderMEs(expcli)
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap",
                      marHeatmap = c(6,6,2,3), 
                      plotDendrograms = F, 
                      xLabelsAngle = 90)

trait<-traitData[,c(6,9:13)]
r<-cbind(c('Stage I','Stage IA','Stage IB','Stage II','Stage IIA','Stage IIB','Stage III',
           'Stage IIIA','Stage IIIB','Stage IIIC','Stage IV','Stage X'),
         c(1:14))
r<-data.frame(r)
trait$Stages<-r[match(trait$tumor_stage,r$X1),'X2']
trait$age<-clid_t[rownames(trait),'age']
trait$Sex<-ifelse(trait$gender=='Male',1,0)
trait<-trait[,c(1,3,5,6,7)]
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, trait, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, trait, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

pdf('Module-trait relationships.pdf')
par(mar=c(4, 6, 3, 1))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(trait), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
# Select module
module = "green";
# Select module probes
probes = colnames(dataExpr) 
inModule = (moduleColors==module)
modProbes = probes[inModule]

########find different color modules' GO pathway ###############################
library(clusterProfiler);library(org.Hs.eg.db)
for (thecolor in c('green')) {
  
  module = thecolor
  # Select module probes
  probes = colnames(dataExpr) ## 
  inModule = (moduleColors==module)
  modProbes = probes[inModule]
  
  modProbes_choosen =mapIds(org.Hs.eg.db,keys = modProbes,keytype = 'SYMBOL',column = 'ENTREZID')
  go_ALL <- enrichGO(gene = modProbes_choosen,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pvalueCutoff = 0.5,
                     qvalueCutoff = 0.5)
  
  kegg_all<- enrichKEGG(
    gene          = modProbes_choosen,
    keyType     = "kegg",
    organism   = 'hsa',
    pvalueCutoff      = 0.05,
    pAdjustMethod     = "BH",
    qvalueCutoff  = 0.1
  )
  BP_df<- go_ALL@result %>% arrange(desc(Count)) %>% filter(ONTOLOGY=="BP") %>% head()
  CC_df<- go_ALL@result %>% arrange(desc(Count)) %>% filter(ONTOLOGY=="CC") %>% head()
  MF_df<- go_ALL@result %>% arrange(desc(Count)) %>% filter(ONTOLOGY=="MF") %>% head()
  KEGG_df<-kegg_all@result %>% arrange(desc(Count)) %>% head()
  KEGG_df$ONTOLOGY<-'KEGG'
  go_df <- rbind(BP_df,CC_df,MF_df,KEGG_df)
  BP_df<- go_ALL@result %>% arrange(desc(Count)) %>% filter(ONTOLOGY=="BP") %>% head()
  CC_df<- go_ALL@result %>% arrange(desc(Count)) %>% filter(ONTOLOGY=="CC") %>% head()
  MF_df<- go_ALL@result %>% arrange(desc(Count)) %>% filter(ONTOLOGY=="MF") %>% head()
  KEGG_df<-kegg_all@result %>% arrange(desc(Count)) %>% head()
  KEGG_df$ONTOLOGY<-'KEGG'
  
  go_df <- rbind(BP_df,CC_df,MF_df,KEGG_df)
  go_df$go_term_order=factor(x = c(1:nrow(go_df)),labels = go_df$Description)
  
  library(ggsci)
  library(ggplot2)
  g<-ggplot(data=go_df, aes(x=go_term_order,y=Count, fill=ONTOLOGY)) +
    geom_bar(stat="identity", width=0.8)  + 
    scale_fill_jama() + 
    theme_classic() +
    xlab(NULL) + ylab("Number of Genes") + labs(title = paste0("The Most Enriched terms of ME",thecolor))+ 
    theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 80,vjust = 1, hjust = 1 )) 
  ggsave(g,filename = paste0('ME',thecolor,' GO KEGG.pdf'))
}

print(g)
HubGenes <- chooseTopHubInEachModule(dataExpr,moduleColors)
write.csv(HubGenes,file='BLCA_FLPS_hub.csv')
