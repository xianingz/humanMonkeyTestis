### Monkey Pericyte subclustering
#  3.24.2019 by Qianyi
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering 
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing doublets
### Removed Cluster 7/7 of SCyte-Somatic doublets from somatic subclustering by CCA
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing C1 and C7 doublets
### Removed 239 cells with >=10% MT transcripts that we failed to exclude as we missed 7 MT gene symbols using Adrienne MT gene list
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7-No10pctMT cells 
### Extract Pericyte, directly-merged and do subclustering
#-> 2 subclusters (Pericytes + Subcluster2)
# note: here I do not use CCA for Pericyte subclustering as there are too few cells for each subject

### set paths and load library
home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/monkey"
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
myBrewerPalette1=gg_color_hue(5)
datainfo=read.table("datainfo_monkey")
datainfo
dataset=as.character(datainfo[,2])
subject=unique(gsub("-.*","",dataset))
length(dataset) # 9
length(subject) # 5


### load object for CCA of 5 subjects after NoC1-NoC7-No10pctMT
# CCA and somatic subclustering after removing Cluster 1/8 and 7/7 of spermatid-myoid doublets for somatic subclustering 
dgefile=dgename="CCA5Subjects-Somatic-NoC1-NoC7-No10pctMT/"
load(file="Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA-No10pctMT.Robj")
testis # 17251 genes across 2098 samples.
testis14=testis

library(Seurat)
# cowplot enables side-by-side ggplots
library(cowplot)

# setup Seurat objects 
### removed cells with >10% MT from previous CCA and somatic subclustering
table(testis14@ident)
#   1    2    3    4    5    6 
#  37  178  155  362 1324   42
testis=SubsetData(testis14,ident.use=3)
table(testis@ident)
#3 
#155 
table(testis@ident,testis@meta.data$protocol)
#     Monkey1 Monkey2 Monkey3 Monkey4 Monkey5
#      122       6      14      10       3
nCellperGene <- rowSums(testis@data>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # [1] 10454

datalist=list()
for(i in 1:length(subject)){
### Extract Pericyte for each monkey subject
cells.use=names(testis@ident)[which(testis@meta.data$protocol == subject[i])]
dgedata=testis@raw.data[genes.use,cells.use]
  dge <- CreateSeuratObject(raw.data = dgedata)
  dge <- NormalizeData(object = dge)
  dge <- ScaleData(object = dge)
  dge <- FindVariableGenes(object = dge, do.plot = FALSE)
  dge@meta.data[,"protocol"] <- subject[i]
  datalist[[i]]=dge
}
datalist
# note: too few cells for each subject, so cannot use CCA


dgedata2=testis@raw.data[genes.use,names(testis@ident)]
dim(dgedata2) # [1] 10454   155

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project="MonkeyPericyte_DirectMerged5Subjects", min.cells=1, min.genes=1)
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("Monkey_Pericyte_DirectMerged5Subjects_HVG0.2.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # [1] 2133
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 10761 genes across 174 samples     
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
#Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
#    122       6      14      10       3 
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
save(dge,file="MonkeyPericyte_DirectMerged5Subjects-No10pctMT.Robj")
table(dge@meta.data$indiv)
#Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
#    122       6      14      10       3 
table(gsub("_.*","",names(dge@ident)))
#Monkey1.1 Monkey1.2   Monkey2 Monkey3.1 Monkey3.2 Monkey4.1 Monkey4.2 Monkey5.2 
#       68        54         6         9         5         7         3         3 

### PCA
print(Sys.time())  #[1] "2019-03-24 12:00:42 EDT"
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)

### PCA Plot
pdf("PCA_Pericyte_DirectMerged5Subjects.pdf",width=12,height=10)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
### Scree Plot for PCA
numPCs=3;i=1 # HVG
pdf("dge_PCA_Variablel_variation.pdf",height=4.5,width=4.5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@dr$pca@sdev,type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
eigenvalue=dge@dr$pca@sdev[numPCs[i]]
text(numPCs[i]+1,eigenvalue+2,col="red",paste(numPCs[i],"PCs"))
#legend("topright",legend=infile[i,2],cex=1.5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
eigenvalue=dge@dr$pca@sdev[numPCs[i]]
print(eigenvalue)
plot(density(dge@dr$pca@sdev),col="red",lwd=2,xlim=c(0,4),xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@dr$pca@sdev),col="black")
lines(density(dge@dr$pca@sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=(eigenvalue+dge@dr$pca@sdev[numPCs[i]+1])/2,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,0.3,col="red",paste(numPCs[i],"PCs"))
plot(density(dge@dr$pca@sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@dr$pca@sdev),col="black")
lines(density(dge@dr$pca@sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=(eigenvalue+dge@dr$pca@sdev[numPCs[i]+1])/2,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,0.3,col="red",paste(numPCs[i],"PCs"))
#legend("topright",legend=infile[i,2],cex=1.5)
dev.off()

### Louvain-Jaccard clustering, tSNE, UMAP using top PCs
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,1,by=0.1), save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])

save(dge,file="MonkeyPericyte_DirectMerged5Subjects-No10pctMT.Robj")


### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
# 1 2 2 2 2 2 3 3 3 3

table(dge@meta.data$res.0.1,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.2,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.3,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.4,dge@meta.data$res.0.5)
# exactly the same two clusters using different resolutions

### decide to use 2 clusters, double-check the number of clusters
res="res.0.5";j=1;resi=1;
   dge=SetAllIdent(dge,id=res)
table(dge@ident)
table(testis1@ident[names(dge@ident)],dge@ident)
                0   1
  Pericyte    116   0
  Subcluster2   1  32

### label cell types for Pericyte subclustering
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident2[which(ident2==0)]<-"m-Pericyte"
ident2[which(ident2==1)]<-"f-Pericyte"
dge=AddMetaData(dge,ident2,"CellType2")
dge=SetAllIdent(dge,id="CellType2")
save(dge,file="MonkeyPericyte_DirectMerged5Subjects-No10pctMT.Robj")

pdf(paste0(dgefile,"PCA_tSNE_UMAP_Pericyte2.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

### Add these two clusters back to All Somatic cells of Monkey
testis=testis14
ident=as.character(testis14@ident)
names(ident)=names(testis14@ident)
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident[names(ident2)]<-ident2
levels=c(levels(testis14@ident)[1:2],"m-Pericyte","f-Pericyte",levels(testis14@ident)[4:6])
#levels=c(levels(testis14@ident)[1:3],"m-Pericyte","f-Pericyte",levels(testis14@ident)[6:8])
table(ident)[levels]
#Macrophage/Tcell      Endothelial       m-Pericyte       f-Pericyte 
#              37              178              117               38 
#           Myoid        ImmLeydig       DiffLeydig 
#             362             1324               42 
ident=factor(ident,levels=levels)
testis=AddMetaData(testis,ident,"CellType2")

testis=SetAllIdent(testis,id="CellType2")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType2))
testis14=testis
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA-No10pctMT.Robj"))
