### Monkey Immune cells subclustering
# 3.11.2019 by Qianyi
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering 
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing doublets
### Removed Cluster 7/7 of SCyte-Somatic doublets from somatic subclustering by CCA
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing C1 and C7 doublets
### Extract Macrophage/Tcell cluster, directly-merged and do subclustering
#-> Decided to use Macrophage/Tcell clusters from 3-species somatic clustering, as that is more accurate
### No cell with >=10% MT genes in the Macrophage/Tcell immune subset 

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



### load object of CCA and somatic subclustering after removing Cluster 1/8 and 7/7 of spermatid-myoid doublets for somatic subclustering 
dgefile=dgename="CCA5Subjects-Somatic-NoC1-NoC7/"
load(file="Monkey5Somatic-NoC1-NoC7_CCA.Robj")
testis # 17469 genes across 2337 samples.
testis1=testis

library(Seurat)
# cowplot enables side-by-side ggplots
library(cowplot)


# setup Seurat objects 
### remove 7/7 from previous CCA and somatic subclustering
table(testis1@ident)
#Macrophage/Tcell      Endothelial         Immune            Myoid 
#              37              170              174              516 
#       ImmLeydig       DiffLeydig 
#            1396               44 
testis=SubsetData(testis1,ident.use="Macrophage/Tcell")
table(testis@ident)
#Macrophage/Tcell 
#              37
table(testis@ident,testis@meta.data$protocol)
#                 Monkey1 Monkey2 Monkey3 Monkey4 Monkey5
#Macrophage/Tcell       4       2      19       6       6
nCellperGene <- rowSums(testis@data>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # [1] 6713

dgedata2=testis@raw.data[genes.use,names(testis@ident)]
dim(dgedata2) #  [1] 6713   37

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project="MonkeyImmune_DirectMerged5Subjects", min.cells=1, min.genes=1)
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("Monkey_Immune_DirectMerged5Subjects_HVG0.2.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # [1] 2669
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 6713 genes across 37 samples     
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
#Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
#   4       2      19       6       6
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")
table(dge@meta.data$indiv)
#Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
#    4       2      19       6       6
table(gsub("_.*","",names(dge@ident)))
#Monkey1.1 Monkey1.2   Monkey2 Monkey3.1 Monkey3.2 Monkey4.1 Monkey4.2 Monkey5.1 
#        1         3         2         8        11         5         1         3 
#Monkey5.2 
#        3 

### PCA
print(Sys.time())  # [1] "2019-03-11 09:57:43 EDT"
#note: given the small sample size (N=37), better to use standard svd
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 20,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)

### PCA Plot
pdf("PCA_Immune_DirectMerged5Subjects.pdf",width=12,height=10)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
### Scree Plot for PCA
numPCs=4;i=1 # HVG
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
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,1.3,by=0.1),force.recalc = TRUE, save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T,perplexity=5)
# if Error: Perplexity is too large, use a smaller perplexity
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])
Sys.time()  # [1] "2019-03-11 10:09:33 EDT"

save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")


### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
# 1 1 1 1 1 1 1 1 1 2 4 2

table(dge@meta.data$res.1)
# 0  1 
#24 13
table(dge@meta.data$res.1.1)
table(dge@meta.data$res.1,dge@meta.data$res.1.1)
#     0  1  2  3
#  0  0  3 21  0
#  1 11  0  0  2


### decide to use 2 clusters, double-check the number of clusters
res="res.1";j=1;resi=1;
   dge=SetAllIdent(dge,id=res)
   print(length(unique(dge@ident)))
   TSNEPlot(dge)

## order cell by cluster ID and randomly shuffle cells within each batch
levels=levels(dge@ident)
ident=factor(dge@ident,levels=levels)

### randomly shuffling cells within each cluster
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

### for each cluster, calculate average normalized expression of each gene
tmpdge=data.frame(t(as.matrix(dge@data[,cells.use])))
# make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dge@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dge@data[,cells.use])))
genecountsall=matrix(,dim(dge@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dge@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}
write.table(genecountsall,paste0(dgefile,"NoC1-NoC7-Immune_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
tmp=genecountsall[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) # 91
da
###### plot with seriation
pdf(file=paste0("Immune_DirectMerged5SubjectsSubjects_Centroid_norm_Seriation_",res,".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO"),col=bluered(10)))
dev.off()

### get order of seriation
do=seriate(da,method="OLO")
print(get_order(do))
levelss=get_order(do)

#levels=rev(levelss)
levels=levelss
print(levels-1)

### Reordered clusters for all cells
cells.use=colnames(dge@data)
# random shuffling cells within ordered clusters
ident=factor(dge@ident,levels=levels-1)

cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

### save ordered cluster ID in dge object
print(which(unique(cells)!=levelss[[j]]-1)) # integer(0)

ordered=paste0(res,"order")

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@meta.data[,ordered]=factor(dge@meta.data[,ordered])
dge=SetAllIdent(dge,ordered)
# save the dge file


### Comparison with separeted Macrophage and T cell for Monkey from CCA of 3-species Somatic
load(file=paste0("../Somatic_3Species_11Subjects_CCA.Robj"))
testisall=testis
testis # 10765 genes across 11158 samples.

table(testisall@ident[names(dge@ident)],dge@ident)
#      1  2
#  9  17  0
#  10  7 13


pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()


## label cell types for Immune subclustering
### 1. label cell types from clustering of 3-species somatic merged by CCA
# note: decided to keep this cluster solution as final separation for monkey Tcell/macrophag
table(testisall@ident[names(dge@ident)])
# 9 10  
#17 20  
ident3=as.character(testisall@ident[names(dge@ident)])
names(ident3)=names(dge@ident)
ident3[which(ident3==9)]<-"Macrophage"
ident3[which(ident3==10)]<-"Tcell"
dge=AddMetaData(dge,ident3,"CellType1")
dge=SetAllIdent(dge,id="CellType1") # used this
save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")


pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune1.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Immune1_3Species_markersall_mindiff0.2_logfc2fold_3.11.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
#Macrophage      Tcell 
#        17         20 
#        41         13 

### 2. label cell types from subclustering by Louvain-Jaccard clustering
table(dge@meta.data$res.1order)
# 1  2 
#24 13
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident2[which(ident2==1)]<-"Macrophage"
ident2[which(ident2==2)]<-"Tcell"
dge=AddMetaData(dge,ident2,"CellType2")
dge=SetAllIdent(dge,id="CellType2")
save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")

pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune2.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Immune2_res.1order_markersall_mindiff0.2_logfc2fold_3.11.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
# 1  2 
#24 13
#32 32 

### 3. label cell types from k-means of UMAP
kclust=kmeans(dge@dr$umap@cell.embeddings,centers=2)
table(kclust$cluster)
# 1  2 
#19 18 

ident3=kclust$cluster
ident3[which(ident3==2)]<-"Macrophage"
ident3[which(ident3==1)]<-"Tcell"
dge=AddMetaData(dge,ident3,"CellType3")
dge=SetAllIdent(dge,id="CellType3")
save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")

pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune3.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Immune3_kmeans2_markersall_mindiff0.2_logfc2fold_3.11.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
#Macrophage      Tcell 
#        18         19 
#        42         11 



### Add these two clusters back to All Somatic cells of Monkey
dge=SetAllIdent(dge,id="CellType1") # used this
testis=testis1
testis=SetAllIdent(testis,id="CellType2")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType2))
ident=as.character(testis@ident)
names(ident)=names(testis@ident)
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident[names(ident2)]<-ident2
levels=c(levels(dge@ident)[2:1],levels(testis@ident)[2:7])
table(ident)[levels]
#  Tcell       Macrophage Endothelial    Pericyte Subcluster2       Myoid 
#         20          17         170         131          43         516 
#  ImmLeydig  DiffLeydig 
#       1396          44 

ident=factor(ident,levels=levels)
testis=AddMetaData(testis,ident,"CellType3")

testis=SetAllIdent(testis,id="CellType3")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType3))
testis1=testis
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA.Robj"))
