### Somatic SubClustering for Liger-merged 4 subjects After Removing Doublets Clusters 5 
# 1.8.2019 by Qianyi
### Comparing and contrasting heterogeneous single cell profiles using liger

### set paths and load library
home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/human"
setwd(home)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(cowplot)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
redblue100<-rgb(read.table('../redblue100.txt',sep='\t',row.names=1,header=T))
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette1=gg_color_hue(4)
datainfo=read.table("datainfo_human")
datainfo
dataset=as.character(datainfo[,2])
subject=unique(gsub("-.*","",dataset))
### Use Human Subject 1, 2, 3, and 5 - exclude 4 and 6
dataset=dataset[c(1:9,11)]
subject=unique(gsub("-.*","",dataset))
length(dataset) # 10
length(subject) # 4

library(liger)
library(cowplot)


# 1.8.2019 by Qianyi
### load object for directly merged 4 subjects after removing Cluster 5
load(file="HumanMerged4-1235-NoC5Somatic.Robj")
dge4=dge
dge # 35063 genes across 3765 samples.

### load object for CCA of 4 subjects after removing Cluster 5
load(file="Human4-1235-NoC5Somatic_CCA.Robj")
pbmc

dgefile=dgename="LIGER/LIGER_"

## Data Preprocessing
### The algorithm takes a list of two or more digital gene expression (DGE) matrices as input. Genes should be in rows and cells in columns. Before running the factorization, we need to normalize the data to account for different numbers of UMIs per cell, select variable genes, and scale the data. Note that we do not center the data because nonnegative matrix factorization accepts only positive values. The selectGenes function performs variable gene selection on each of the datasets separately, then takes the union. Note that coresponding genes in each dataset need to have the same names (though the genes do not need to be in the same order in each dataset). For cross-species analysis, it may be convenient to convert all gene names to uppercase; you can do this using the capitalize=T option of the selectGenes function.

setlist=list()
for(i in 1:length(subject)){
### Extract Somatic cells for each human subject
cells.use=names(dge4@ident)[which(dge4@meta.data$indiv == subject[i])]
dgedata=dge4@raw.data[,cells.use]
nCellperGene <- rowSums(dgedata>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgedata2=dgedata[genes.use,]
setlist[[i]]=dgedata2
}
names(setlist)=subject

## Create liger object, normalize, select highly-variable genes, scale not center
ligerex = createLiger(setlist) #Can also pass in more than 2 datasets
ligerex = normalize(ligerex)
pdf(file=paste0(dgefile,"selectGenes.pdf"),height=5,width=4.5)
ligerex = selectGenes(ligerex, var.thresh = 0.1)
dev.off()
ligerex = scaleNotCenter(ligerex)


## Selecting k and lambda
### The suggestK and suggestLambda functions can aid in selecting k and lambda. 
### We want to find the smallest k for which the increase in entropy metric begins to level off (an “elbow” in the plot). 
### Similarly, we want the smallest lambda for which the alignment metric stabilizes.
pdf(file=paste0(dgefile,"suggestK.pdf"))
k.suggest <- suggestK(ligerex, num.cores = 8, gen.new = T, return.results = T)
dev.off()
pdf(file=paste0(dgefile,"suggestLamda_k20.pdf"))
lamda.suggest <- suggestLambda(ligerex, k=20, num.cores = 8, gen.new = T, return.results = T)
dev.off()
pdf(file=paste0(dgefile,"suggestLamda_k10.pdf"))
lamda.suggest <- suggestLambda(ligerex, k=10, num.cores = 8, gen.new = T, return.results = T)
dev.off()


## Performing the Factorization
### Next we perform the factorization using an alternating least squares algorithm. 

## tSNE

## SNF clustering and quantile alignment
### After performing the factorization, we identify cells that load on corresponding cell factors and quantile normalize their factor loadings across datasets. 
### The key parameters here are the number of factors (k), the penalty parameter (lambda), and the clustering resolution. 
### In most cases, the default settings of lambda=5.0 and resolution=1.0 provide reasonable results.

ligerex = optimizeALS(ligerex, k = 20,lambda=5) 
ligerex = runTSNE(ligerex)
ligerex = quantileAlignSNF(ligerex,resolution=0.5) 
table(ligerex@clusters)
   0    1    2    3    4    5    6 
1090  561  544  516  473  292  289
table(pbmc@ident,ligerex@clusters)
       0    1    2    3    4    5    6
  1    9    3   26    0    1    4    0
  2   26    4    3    1  430    1    2
  3 1014  553  492   31   38   10    5
  4   17    1   19  482    1    0    1
  5    3    0    3    2    3    1  281
  6    2    0    0    0    0  231    0
  7   19    0    1    0    0   45    0
save(ligerex,file="Human4-1235-NoC5Somatic_LIGER_k20.Robj")

ligerex = optimizeALS(ligerex, k = 10,lambda=10) 
ligerex = runTSNE(ligerex)
ligerex = quantileAlignSNF(ligerex,resolution=0.5) 
table(ligerex@clusters)
  0   1   2   3   4   5   6 
913 716 537 523 496 294 286
table(pbmc@ident,ligerex@clusters)
      0   1   2   3   4   5   6
  1   3   4   4   2   0   0  30
  2   1   9 431  21   0   5   0
  3 905 683  85 429  20  11  10
  4   3  18   8  15 475   1   1
  5   1   2   9   1   1 277   2
  6   0   0   0   1   0   0 232
  7   0   0   0  54   0   0  11
save(ligerex,file="Human4-1235-NoC5Somatic_LIGER_k10.Robj")

table(ligerex@clusters,dge4@ident)
table(ligerex@clusters[names(dge5@ident)],dge5@meta.data$res.0.6order)
table(ligerex@clusters,pbmc@ident)


## Visualizing the results
pdf(file=paste0(dgefile,"word_clouds.pdf"),height=6)
plotWordClouds(ligerex)
dev.off()

p<-plotByDatasetAndCluster(ligerex,return.plots=T) #Can also pass in different set of cluster labels to plot

pdf(file=paste0(dgefile,"tSNE.pdf"),height=4,width=5.5)
print(p[[1]])
dev.off()

pdf(file=paste0(dgefile,"clusters.pdf"),height=4,width=5)
print(p[[2]]+scale_color_brewer(type='qual',palette='Paired'))
dev.off()


## Finding marker genes
### We can use the factorization to identify shared and dataset-specific markers. The function below returns a list, where the first element contains dataset-specific markers for dataset 1, the second element contains shared markers, the third element contains dataset-specific markers for dataset 2, and the last 2 elements indicate the number of factors in which each marker is found. This information allows the identification of ubiquitous vs. cell-type-specific dataset differences.
### note: can only get shared markers for two datasets
markers = getFactorMarkers(ligerex, dataset1="Human5",dataset2="Human2",num.genes = 10)
markers$shared[1:5,]

# Show top 10 highly loading dataset specific and shared genes
# Prevent text progress bar from printing
library(ggrepel)
library(grid)

word_clouds <- plotWordClouds(ligerex, dataset1="Human5",dataset2="Human2",num.genes = 10, do.spec.plot = F,
return.plots = T)
pdf(file=paste0(dgefile,"word_clouds_markers_k10_Human5.pdf"),width=5,height=8)
for(i in 1:length(word_clouds)){
	print(word_clouds[[i]])
}
dev.off()

plotGene(ligerex, gene = "CD163")

pdf(file=paste0(dgefile,"knownmarker.pdf"),height=10,width=22)
plotGeneViolin(ligerex, gene = "CD163")
plotGeneViolin(ligerex, gene = "VWF")
plotGeneViolin(ligerex, gene = "ACTA2")
plotGeneViolin(ligerex, gene = "PDGFRA")
plotGeneViolin(ligerex, gene = "DLK1")
dev.off()

