### Examining Cross-subjects Batch effect
# by Qianyi on 12.18.2018
### Merged replicates for each single subject, then cross-tabulated cluster centroids for each subject 

### set path and load libraries
home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/monkey"
setwd(home)
library(Seurat)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
redblue100<-rgb(read.table('../redblue100.txt',sep='\t',row.names=1,header=T))
datainfo=read.table("datainfo_monkey")
datainfo
dataset=as.character(datainfo[,2])
subject=unique(gsub("-.*","",dataset))

# 12.18.2018 by Qianyi

# clustering for each individual subject, order clusters, flip coordinates to make consistency
### load object for merged replicates of each individual subject
dgealllist=list()
for(i in 1:length(subject)){
  load(file=paste0(subject[i],".Robj"))
  dgealllist[[i]]=dgeall
  print(dgeall)
}
project Monkey1 19304 genes across 2108 samples.
project Monkey2 19496 genes across 2804 samples.
project Monkey3 21546 genes across 5058 samples.
project Monkey4 20811 genes across 6625 samples.
project Monkey5 20984 genes across 5645 samples.
for(i in 1:length(subject)){
dge=dgealllist[[i]]
print(dge)
}

### Number of highly-variable genes
for(i in 1:length(subject)){
dge=dgealllist[[i]]
print(length(dge@var.genes))
}
[1] 2646
[1] 2666
[1] 2780
[1] 2605
[1] 2315
# check the number of hvg overlapped between datasets
cc=matrix(0,length(subject),length(subject))
for(i in 1:length(subject)){
    for(j in 1:length(subject)){
        cc[i,j]=length(intersect(dgealllist[[i]]@var.genes,dgealllist[[j]]@var.genes))
    }
}
cc
     [,1] [,2] [,3] [,4] [,5]
[1,] 2646 1489 1567 1416 1108
[2,] 1489 2666 1936 1802 1523
[3,] 1567 1936 2780 1768 1433
[4,] 1416 1802 1768 2605 1411
[5,] 1108 1523 1433 1411 2315

### Add a column of subject ID in meta data sheet
for(i in 1:length(subject)){
dge=dgealllist[[i]]
ident<-rep(subject[i],length(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
dgeall=dge
dgealllist[[i]]=dgeall
}
### PCA
print(Sys.time())   # [1] "2018-12-17 11:35:08 EST"
for(i in 1:length(subject)){
dge=dgealllist[[i]]
dge <- PCA(dge, do.print = TRUE, pcs.print = 1, genes.print = 5)
dgeall=dge
dgealllist[[i]]=dgeall
save(dgeall,file=paste0(subject[i],".Robj"))
}
### PCA Plot
plotlist2=list()
for(i in 1:length(subject)){
dge=dgealllist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
}
pdf("Monkey5_PCA.pdf",height=5.5,width=12)
MultiPlotList(plotlist2,cols = 3)
dev.off()
### Scree Plot for PCA
numPCs=c(8,7,9,8,7)
pdf("Monkey5_dge_PCA_Variablel_variation.pdf",width=7.5,height=5)
par(mfrow=c(2,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
for(i in 1:length(subject)){
dge=dgealllist[[i]]
plot(dge@pca.obj[[1]]$sdev[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,5,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=subject[i])
}
dev.off()
### Louvain-Jaccard clustering and tSNE
for(i in 1:length(subject)){
dge=dgealllist[[i]]
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(0.1,2,0.1), print.output = 0, save.SNN = T)
dge=RunTSNE(dge,dims.use = 1:numPCs[i],do.fast=T)   
dgeall=dge
dgealllist[[i]]=dgeall
save(dgeall,file=paste0(subject[i],".Robj"))
}


## find approximately the same number of clusters for each batch
### check the number of clusters
for(i in 1:length(subject)){
   dge=dgealllist[[i]]
   print(c( length(unique(dge@data.info$res.0.1)),length(unique(dge@data.info$res.0.2)),length(unique(dge@data.info$res.0.3)),length(unique(dge@data.info$res.0.4)),length(unique(dge@data.info$res.0.5)),length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)) ))
} 
i=5
   dge=dgealllist[[i]]
   print(c( length(unique(dge@data.info$res.1.3)),length(unique(dge@data.info$res.1.4)),length(unique(dge@data.info$res.1.5)),length(unique(dge@data.info$res.1.5)),length(unique(dge@data.info$res.1.7)),length(unique(dge@data.info$res.1.8)),length(unique(dge@data.info$res.1.9)),length(unique(dge@data.info$res.2))))

### Compare clusters with individual replicate clustering 
###### clustering for individual replicate used all filtered genes, 
###### clustering for each subject of all replicates used highly-variable genes
dgelist=list()
for(i in 10:11){
   load(file=paste0(datainfo[i,3],".Robj"))
   dgelist[[i]]=dge
}
table(dgealllist[[4]]@data.info$res.0.2,dgelist[[10]]@data.info$res.0.2)
table(dgealllist[[5]]@data.info$res.0.5,dgelist[[11]]@data.info$res.0.3)

### decide to use 11 clusters
res=paste0("res.",c(0.6,0.8,0.6,0.7,0.5))
### double-check that I've set the same number of clusters for 5 monkey subjects
for(i in 1:length(subject)){
   dge=dgealllist[[i]]
   dge=SetAllIdent(dge,id=res[i])
   print(length(unique(dge@ident)))
}

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()

for(j in 1:length(subject)){
dge=dgealllist[[j]]
dge=SetAllIdent(dge,id=res[j])
print(c(subject[j],length(unique(dge@ident))))
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
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) expMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

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
pdf(file=paste0("Centroid_norm_Seriation_",subject[j],"_",res[j],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO"),col=bluered(10)))
dev.off()

### get order of seriation
do=seriate(da,method="OLO")
print(get_order(do))
levelss[[j]]=get_order(do)

### levels
levels=levelss[[j]]
if(j==1 | j==5){
  levels=rev(levelss[[j]])
}
if(j==3){
  levels=levelss[[j]][c(4:1,11:5)]
}
if(j==4){
  levels=levelss[[j]][c(10,11,9:1)]
}

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

ordered="ordered"

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@data.info[,ordered]=factor(dge@data.info[,ordered])
dge=SetAllIdent(dge,ordered)
# save the dge file
dgealllist[[j]]=dge
#save(dge,file=paste0(dataset[j],".Robj"))
}


## double-check if I need to reverse the cluster ID orders
### use PRM2 as marker for the last cluster (elongating)
pdf("clusters_ordered0_PRM2_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"PRM2")
}
dev.off()
pdf("clusters_ordered0_ACTA2_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"ACTA2")
}
dev.off()
pdf("clusters_ordered0_ID4_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"ID4")
}
dev.off()
pdf("clusters_ordered0_PRM2_Feature.pdf",height=4,width=4)
for(i in 1:length(subject)){
  FeaturePlot(dgealllist[[i]],"PRM2")
}
dev.off()
#FeaturePlot(dgealllist[[i]],"PRM2")
#TSNEPlot(dgealllist[[i]],do.return=T)
# download the tSNE plot below to check the cluster IDs

## double-check if I need to flip PCs, or flip tSNE
### PCA and tSNE for ordered clusters of each individual replicate
plotlistt=plotlist2=plotlist3=plotlist4=plotlist5=list()
for(i in 1:length(subject)){
dge=dgealllist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist4[[i]]=PCAPlot(dge,1,4,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist5[[i]]=PCAPlot(dge,1,5,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE,label.size=4)
}
pdf("clusters_ordered0.pdf",height=5,width=10)
MultiPlotList(plotlist2,cols = 3)
MultiPlotList(plotlist3,cols = 3)
MultiPlotList(plotlist4,cols = 3)
MultiPlotList(plotlist5,cols = 3)
MultiPlotList(plotlistt,cols = 3)
dev.off()

### switch PC1&2
for(i in c(3,4)){
dge=dgealllist[[i]]
tmp=dge@pca.rot[,1]
dge@pca.rot[,1]=dge@pca.rot[,2]
dge@pca.obj[[1]]$x[,1]=dge@pca.rot[,2]
dge@pca.rot[,2]=tmp
dge@pca.obj[[1]]$x[,2]=tmp
dgealllist[[i]]=dge
}
### flip PC1
for(i in c(2,3,4,5)){
dge=dgealllist[[i]]
dge@pca.rot[,1]=-dge@pca.rot[,1]
dge@pca.obj[[1]]$x[,1]=-dge@pca.obj[[1]]$x[,1]
dgealllist[[i]]=dge
}
### flip PC2
for(i in c(2,5)){
dge=dgealllist[[i]]
dge@pca.rot[,2]=-dge@pca.rot[,2]
dge@pca.obj[[1]]$x[,2]=-dge@pca.obj[[1]]$x[,2]
dgealllist[[i]]=dge
}
### flip PC3
for(i in c(1,3,4,5)){
dge=dgealllist[[i]]
dge@pca.rot[,3]=-dge@pca.rot[,3]
dge@pca.obj[[1]]$x[,3]=-dge@pca.obj[[1]]$x[,3]
dgealllist[[i]]=dge
}
### flip PC4
for(i in c(1,3,5)){
dge=dgealllist[[i]]
dge@pca.rot[,4]=-dge@pca.rot[,4]
dge@pca.obj[[1]]$x[,4]=-dge@pca.obj[[1]]$x[,4]
dgealllist[[i]]=dge
}
### switch tSNE1 and tSNE2
for(i in c(1,2,5)){
dge=dgealllist[[i]]
tmp=dge@tsne.rot[,2]
dge@tsne.rot[,2]=dge@tsne.rot[,1]
dge@tsne.rot[,1]=tmp
dgealllist[[i]]=dge
}
### flip tSNE1
for(i in c(5)){
dge=dgealllist[[i]]
dge@tsne.rot[,1]=-dge@tsne.rot[,1]
dgealllist[[i]]=dge
}
### flip tSNE2
for(i in c(3)){
dge=dgealllist[[i]]
dge@tsne.rot[,2]=-dge@tsne.rot[,2]
dgealllist[[i]]=dge
}



## after reverse cluster order and flip cooridnates, repeat the above plots
# change ordered0 in the file name to order1
pdf("clusters_order1_PRM2_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"PRM2",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_ACTA2_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"ACTA2",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_CLU_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"CLU",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_NR2F2_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"NR2F2",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_STRA8_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"STRA8",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_ID4_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"ID4",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_ZBTB16_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"ZBTB16",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_ACTA2_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"ACTA2",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_CLU_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"CLU",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_NR2F2_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"NR2F2",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_STRA8_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"STRA8",cols.use=myBrewerPalette)
}
dev.off()
pdf("clusters_order1_CYP17A1_Violin.pdf",height=4,width=8)
for(i in 1:length(subject)){
  VlnPlot(dgealllist[[i]],"CYP17A1",cols.use=myBrewerPalette)
}
dev.off()
knownmarkers=c("VIM","CD163","S100A4","TYROBP","LYZ","RGS1","VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN","MYH11","ACTA2","DCN","DLK1","NR2F2","CYP17A1","ID4","CLU","PRM2")
pdf("clusters_order1_Feature.pdf",height=12,width=18)
for(i in 1:length(subject)){
  FeaturePlot(dgealllist[[i]],knownmarkers,col=c("grey80","red"),nCol=6)
}
dev.off()

## double-check if I need to flip PCs, or flip tSNE
### PCA and tSNE for ordered clusters 
### Visualize PC1-4 and tSNE for merged replicates of each individual subject
plotlistt=plotlist2=plotlist3=plotlist4=plotlist5=list()
for(i in 1:length(subject)){
dge=dgealllist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist4[[i]]=PCAPlot(dge,1,4,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist5[[i]]=PCAPlot(dge,1,5,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE,label.size=4)
}
pdf("Monkey5_clusters_order1.pdf",height=5,width=8)
MultiPlotList(plotlist2,cols = 3)
MultiPlotList(plotlist3,cols = 3)
MultiPlotList(plotlist4,cols = 3)
MultiPlotList(plotlist5,cols = 3)
MultiPlotList(plotlistt,cols = 3)
dev.off()

for(i in 1:length(subject)){
dgeall=dgealllist[[i]]
save(dgeall,file=paste0(subject[i],".Robj"))
}


# merge all replicates of all 5 subjects 

### load object for merged replicates of each individual subject
dgealllist=list()
for(i in 1:length(subject)){
  load(file=paste0(subject[i],".Robj"))
  dgealllist[[i]]=dgeall
}

### merge all replicates for each individual
set=list()
for(i in 1:length(subject)){
  set[[i]]=data.frame(GENE=rownames(dgealllist[[i]]@data),dgealllist[[i]]@raw.data)
}
dgedata=Reduce(function(x,y) merge(x,y,all=TRUE), set)
dgedata[is.na(dgedata)] <- 0
row.names(dgedata)=dgedata[,1]
dgedata=dgedata[,-1]
dim(dgedata)                   #[1] 23101 22240
nCellperGene <- rowSums(dgedata>0)
length(which(nCellperGene==0))     #0
nCellperGene[which(nCellperGene==0)] #numeric(0)

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata, project="MonkeyMerged5", min.cells=1, min.genes=1)
MT <- c("ENSMMUG00000028704", "ENSMMUG00000028703", "ENSMMUG00000028702", "ENSMMUG00000028701", "ENSMMUG00000028700", "ENSMMUG00000028699", "ENSMMUG00000028698", "ENSMMUG00000028697", "ENSMMUG00000028696", "ENSMMUG00000028694", "ENSMMUG00000028693", "ENSMMUG00000028692", "ENSMMUG00000028691", "ENSMMUG00000028690", "ENSMMUG00000028688", "ENSMMUG00000028687", "ENSMMUG00000028686", "ENSMMUG00000028685", "ENSMMUG00000028684", "ENSMMUG00000028683", "ENSMMUG00000028681", "ENSMMUG00000028679", "ENSMMUG00000028678", "ENSMMUG00000028676", "ENSMMUG00000028675", "ENSMMUG00000028674", "ENSMMUG00000028671", "ENSMMUG00000028670", "ENSMMUG00000028669", "ENSMMUG00000028668", "ENSMMUG00000028695", "MT-CO1", "MT-CO3", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6")
MT=unique(MT)
mito.genes=MT[which(MT %in% rownames(dge@data))]
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("Monkey_Merged5_HVG0.2.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # 1787
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 46781 genes across 14058 samples           
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 Monkey6 
  4928   4194   2481   3979   2455   4626
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
dgeall=dge
dge5=dge
save(dgeall,file="MonkeyMerged5.Robj")
table(gsub("_.*","",names(dge@ident)))
Monkey1.1 Monkey1.2 Monkey1.3 Monkey1.4 Monkey1.5 Monkey2.1 Monkey2.2 Monkey3.1 
     888     1411      776      890      963     2199     1995      993 
Monkey3.2   Monkey5 
    1488     2455 

### PCA
print(Sys.time())   
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)
dge5=dge
save(dge,file="MonkeyMerged5.Robj")
print(Sys.time())   


# cross-tabulation for all individual subjects (merged replicates)
### load object for merged 5 subjects
load(file="MonkeyMerged5.Robj")
dge6=dgeall

### load object for each individual subject
dgealllist=list()
for(indiv in 1:length(subject)){
  load(file=paste0(subject[indiv],".Robj"))
  dgealllist[[indiv]]=dgeall
}


## Monkey1-5 rank cor
### order cells by batch first, then by clusters of each batch
blockident=NULL
for(i in 1:length(subject)){
  tmp=paste(subject[i],dgealllist[[i]]@ident,sep="_")
  names(tmp)=names(dgealllist[[i]]@ident)
  blockident=c(blockident,tmp)
}

### Clusters ordered first by batches, then by res
batch=subject
nbatch=length(batch)
ncluster=NULL
for(i in 1:length(subject)){
  ncluster=c(ncluster,length(unique(dgealllist[[i]]@ident)))
}
ncluster 
clusters=list()
for(i in 1:nbatch){
  clusters[[i]]=rep(1:ncluster[i])
}
levels2=NULL
for(bb in 1:nbatch){
    cluster=clusters[[bb]]
    levels2=c(levels2,paste(batch[bb],cluster,sep="_"))
}
levels2=levels2[which(levels2 %in% unique(blockident))]
levels=levels2

ident=factor(blockident,levels=levels)

### order cells by batch, then by clusters and randomly shuffle cells within each cluster
cells=sort(ident)
cells.use=NULL
for(sample in 1:length(levels)){
   set.seed(sample)
   tmp=cells[which(cells == levels[sample])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

rm(dgealllist)
rm(dgeall)
gc()

### Calculate correlation for each normalized centroid using HVG
### for each cluster, calculate average normalized expression of each gene
dge=dge6
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
cc=cor(as.matrix(genecountsall)[dge@var.genes,],method="spearman")
dim(cc)
min(cc) # 0.036

data.use=cc[levels,levels]
### save the cluster centroid rank correlation using HVG
write.table(data.use,"Monkey_5subjects_Centroid_rho_HVG.txt",row.names=T,col.names=T,quote=F,sep="\t")
### note: once calculated, next time can directly read table without calculating again

### load cluster centroid rank correlation using HVG
data.use=read.table("Monkey_5subjects_Centroid_rho_HVG.txt",header=T,row.names=1)
colnames(data.use)=rownames(data.use)
levels=rownames(data.use)
batch=subject

### labeling
colsep.use=cumsum(table(gsub("_.*","",levels)))
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels)))-table(gsub("_.*","",levels))/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)

sidecol=do.call(rbind,strsplit(levels,"_"))
batchtmp=batch[which(batch %in% unique(sidecol[,1]))]
for(rep in 1:length(unique(sidecol[,1]))){
a=batchtmp[rep]
sidecol[which(sidecol[,1]==a),1]<-rep
}

library(RColorBrewer)
myBrewerPalette <- brewer.pal(12,"Paired")
rlab=matrix(0,2,length(levels))
rlab[1,]=rep(c("white","black"),6)[as.numeric(sidecol[,1])]
for(i in 1:nrow(sidecol)){
  rlab[2,i]=myBrewerPalette[as.numeric(sidecol[i,2])]
}
clab=cbind(rlab[2,],rlab[1,])
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

col.use=redblue100

pdf(file="Monkey_5subjects_Centroid_RankedCorrelation_HVG.pdf",height=5.5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()


# cross-tabulation for cluster centroids of merged replicates Vs individual replicate for each individual subject
### load individual replicate clustering
dgelist=list()
for(i in 1:length(dataset)){
  load(file=paste0(datainfo[i,1],".Robj"))
  dgelist[[i]]=dge
}

### load object for each individual subject
dgealllist=list()
for(indiv in 1:length(subject)){
  load(file=paste0(subject[indiv],".Robj"))
  dgealllist[[indiv]]=dgeall
}


## Monkey1-5 number of cells
### cross-tabulate the number of cells
for(indiv in 1:length(subject)){
for(i in 1:length(grep(subject[indiv],dataset))){
  rep=grep(subject[indiv],dataset)[i]
  ident1=dgelist[[rep]]@ident
  cells=gsub("-",".",paste(dataset[rep],names(ident1),sep="_"))
  print(table(ident1,dgealllist[[indiv]]@ident[cells]))
}
}

## Monkey1-5 rank cor
datalist=list()
for(indiv in 1:length(subject)){
### order cells by batch first, then by clusters of each batch
blockident=NULL
for(i in 1:length(grep(subject[indiv],dataset))){
  rep=grep(subject[indiv],dataset)[i]
  tmp=paste(dataset[rep],dgelist[[rep]]@ident,sep="_")
  names(tmp)=paste(dataset[rep],names(dgelist[[rep]]@ident),sep="_")
  blockident=c(blockident,tmp)
}
if(indiv==2){
  tmp=paste(paste0(subject[indiv],"HVG"),dgealllist[[indiv]]@ident,sep="_")
  names(tmp)=gsub("_","._",names(dgealllist[[indiv]]@ident))
  batch=c(grep(subject[indiv],dataset,value=T),paste0(subject[indiv],"HVG"))
} else{
  tmp=paste(subject[indiv],dgealllist[[indiv]]@ident,sep="_")
  names(tmp)=names(dgealllist[[indiv]]@ident)  
  batch=c(grep(subject[indiv],dataset,value=T),subject[indiv])
}
  blockident=c(blockident,tmp)

### Clusters ordered first by batches, then by res
nbatch=length(batch)
ncluster=NULL
for(i in 1:length(grep(subject[indiv],dataset))){
  rep=grep(subject[indiv],dataset)[i]
  ncluster=c(ncluster,length(unique(dgelist[[rep]]@ident)))
}
  ncluster=c(ncluster,length(unique(dgealllist[[indiv]]@ident)))
ncluster 
clusters=list()
for(i in 1:nbatch){
  clusters[[i]]=rep(1:ncluster[i])
}
levels2=NULL
for(bb in 1:nbatch){
    cluster=clusters[[bb]]
    levels2=c(levels2,paste(batch[bb],cluster,sep="_"))
}
levels2=levels2[which(levels2 %in% unique(blockident))]
levels=levels2

ident=factor(blockident,levels=levels)

### order cells by batch, then by clusters and randomly shuffle cells within each cluster
cells=sort(ident)
cells.use=NULL
for(sample in 1:length(levels)){
   set.seed(sample)
   tmp=cells[which(cells == levels[sample])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

### Calculate correlation for each normalized centroid using HVG
### for each cluster, calculate average normalized expression of each gene
dge=dgealllist[[indiv]]
data2=dge@data
if(indiv==2){
  colnames(data2)=gsub("_","._",colnames(data2))
  data=cbind(dge@data,data2)
} else{
  colnames(data2)=gsub("\\.","-",colnames(data2))
  data=cbind(data2,dge@data)
}
tmpdge=data.frame(t(as.matrix(data[,cells.use])))
# make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(data)[cells.use])  # integer(0)
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(data[,cells.use])))
genecountsall=matrix(,dim(dge@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dge@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) expMean(as.numeric(x))) # ExpMean # log(mean(exp(x) - 1) + 1)
    print(i)
}
cc=cor(as.matrix(genecountsall)[dge@var.genes,],method="spearman")
dim(cc)
min(cc) # 0.036

data.use=cc[levels,levels]
### save the cluster centroid rank correlation using HVG
write.table(data.use,paste0(subject[indiv],"_IndivVsMergedReplicates_Centroid_rho_HVG.txt"),row.names=T,col.names=T,quote=F,sep="\t")
### note: once calculated, next time can directly read table without calculating again
datalist[[indiv]]=data.use

### load cluster centroid rank correlation using HVG
data.use=read.table(paste0(subject[indiv],"_IndivVsMergedReplicates_Centroid_rho_HVG.txt"),header=T,row.names=1)
colnames(data.use)=rownames(data.use)
levels=rownames(data.use)
batch=c(grep(subject[indiv],dataset,value=T),subject[indiv])

### labeling
colsep.use=cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])-table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))]/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)

sidecol=do.call(rbind,strsplit(levels,"_"))
batchtmp=batch[which(batch %in% unique(sidecol[,1]))]
for(rep in 1:length(unique(sidecol[,1]))){
a=batchtmp[rep]
sidecol[which(sidecol[,1]==a),1]<-rep
}

library(RColorBrewer)
myBrewerPalette <- brewer.pal(12,"Paired")
rlab=matrix(0,2,length(levels))
rlab[1,]=rep(c("white","black"),6)[as.numeric(sidecol[,1])]
for(i in 1:nrow(sidecol)){
  rlab[2,i]=myBrewerPalette[as.numeric(sidecol[i,2])]
}
clab=cbind(rlab[2,],rlab[1,])
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

col.use=redblue100

pdf(file=paste0(subject[indiv],"_IndivVsMergedReplicates_Centroid_RankedCorrelation_HVG.pdf"),height=5.5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()

}
### Concluded Monkey subjects 1-5 are very similar