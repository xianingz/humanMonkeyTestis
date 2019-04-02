### Examining Within-subject Between-Replicates Batch effect
# by Qianyi on 12.14.2018
### Within each single subject, visualized clusters for each replicate in top PCs and t-SNE spaces and cross-tabulated cluster centroids for each replicate 

### set paths and load library
home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/monkey/"
setwd(home)
library(Seurat)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
redblue100<-rgb(read.table('../redblue100.txt',sep='\t',row.names=1,header=T))

# use the same number of clusters for each replicate, and order the clusters for each batch


# 12.14.2018 by Qianyi 
## load all individual replicate batches
datainfo=read.table("datainfo_monkey")
datainfo
      V1        V2                    V3
1  80188 Monkey1-1 mac1_tsne_res0.8.Robj
2  80189 Monkey1-2 mac2_tsne_res1.5.Robj
3 100262   Monkey2 mac3_tsne_res0.7.Robj
4  99907 Monkey3-1 mac4_tsne_res0.6.Robj
5  94670 Monkey3-2 mac5_tsne_res0.6.Robj
6 109004 Monkey4-1 mac6_tsne_res0.6.Robj
7 109005 Monkey4-2 mac7_tsne_res0.6.Robj
8 109006 Monkey5-1 mac8_tsne_res0.6.Robj
9 109007 Monkey5-2 mac9_tsne_res0.6.Robj


dataset=as.character(datainfo[,2])


## load individual replicate Robj
dgelist=list()
for(i in 1:length(dataset)){
   load(file=paste0(datainfo[i,3]))
   dgelist[[i]]=dge
}

## find approximately the same number of clusters for each batch
### check the number of clusters for each replicate using different resolutions
for(i in 1:length(dataset)){
   dge=dgelist[[i]]
   print(length(unique(dge@ident)))
}
for(i in 1:length(dataset)){
   dge=dgelist[[i]]
   print(c( length(unique(dge@data.info$res.0.1)),length(unique(dge@data.info$res.0.2)),length(unique(dge@data.info$res.0.3)),length(unique(dge@data.info$res.0.4)),length(unique(dge@data.info$res.0.5)),length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)) ))
} 
### most of batches have 11 clusters, so I decide to have 11 clusters for each batch
### re-do clustering for 100262 (Monkey2) to check if I can get 11 clusters
j=3
dge=dgelist[[j]]
dgefile=paste0(dataset[j],"_")
numPCs=7;i=1 # need to modify
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
dev.off()
dge <- FindClusters(dge, pc.use = 1:numPCs, resolution = seq(0.1,0.6,0.1), print.output = 0, save.SNN = T)
print(c( length(unique(dge@data.info$res.0.1)),length(unique(dge@data.info$res.0.2)),length(unique(dge@data.info$res.0.3)),length(unique(dge@data.info$res.0.4)),length(unique(dge@data.info$res.0.5)),length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)) ))
# no, I do not find 11 clusters after changing number of top PCs
### double-check that I've set 11 clusters for 8 batches, and set 10 clusters for batch Monkey2-2
for(i in 1:length(dataset)){
   dge=dgelist[[i]]
   print(length(unique(dge@ident)))
}
res=paste0("res.",c(0.8,1.5,0.7,0.6,0.6,0.6,0.6,0.6,0.6))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
for(j in 1:length(dataset)){
dge=dgelist[[j]]
dge=SetAllIdent(dge,id=res[j])
print(c(dataset[j],length(unique(dge@ident))))
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
pdf(file=paste0("Centroid_norm_Seriation_",dataset[j],"_",res[j],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO"),col=bluered(10)))
dev.off()

### get order of seriation
do=seriate(da,method="OLO")
print(get_order(do))
levelss[[j]]=get_order(do)
levelss[[j]]
levels=levelss[[j]]
print(levels)
if(j<4 || j==6){
  levels=rev(levelss[[j]])
}
if(j==4){
  levels=levelss[[j]][c(2:4,1,11:5)]
}
if(j==7){
  levels=levelss[[j]][c(10:11,1:9)]
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
print(which(unique(cells)!=get_order(do)-1)) # integer(0)

ordered="ordered"

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@data.info[,ordered]=factor(dge@data.info[,ordered])
dge=SetAllIdent(dge,ordered)
# save the dge file
dgelist[[j]]=dge
#save(dge,file=paste0(dataset[j],".Robj"))
}



## double-check if I need to reverse the cluster ID orders
### use PRM2 as marker for the last cluster (elongating)
pdf("clusters_ordered0_PRM2_Violin.pdf",height=4,width=8)
for(i in 1:length(dataset)){
  VlnPlot(dgelist[[i]],"PRM2")
}
dev.off()
pdf("clusters_ordered0_PRM2_Feature.pdf",height=4,width=4)
for(i in 1:length(dataset)){
  FeaturePlot(dgelist[[i]],"PRM2")
}
dev.off()
#FeaturePlot(dgelist[[i]],"PRM1")
#TSNEPlot(dgelist[[i]],do.return=T)
# download the tSNE plot below to check the cluster IDs

## double-check if I need to flip PCs, or flip tSNE
### PCA and tSNE for ordered clusters of each individual replicate
plotlistt=plotlist2=plotlist3=plotlist4=plotlist5=list()
for(i in 1:length(dataset)){
dge=dgelist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist4[[i]]=PCAPlot(dge,1,4,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist5[[i]]=PCAPlot(dge,1,5,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE,label.size=4)
}
pdf("clusters_ordered0.pdf",height=8,width=10.5)
MultiPlotList(plotlist2,cols = 3)
MultiPlotList(plotlist3,cols = 3)
MultiPlotList(plotlist4,cols = 3)
MultiPlotList(plotlist5,cols = 3)
MultiPlotList(plotlistt,cols = 3)
dev.off()

### flip PC1
for(i in c(1)){
dge=dgelist[[i]]
dge@pca.rot[,1]=-dge@pca.rot[,1]
dge@pca.obj[[1]]$x[,1]=-dge@pca.obj[[1]]$x[,1]
dgelist[[i]]=dge
}
### flip PC2
for(i in c(2,4,6,7,9)){
dge=dgelist[[i]]
dge@pca.rot[,2]=-dge@pca.rot[,2]
dge@pca.obj[[1]]$x[,2]=-dge@pca.obj[[1]]$x[,2]
dgelist[[i]]=dge
}
### flip PC3
for(i in c(2,4)){
dge=dgelist[[i]]
dge@pca.rot[,3]=-dge@pca.rot[,3]
dge@pca.obj[[1]]$x[,3]=-dge@pca.obj[[1]]$x[,3]
dgelist[[i]]=dge
}
### flip tSNE1
for(i in c(5,6)){
dge=dgelist[[i]]
dge@tsne.rot[,1]=-dge@tsne.rot[,1]
dgelist[[i]]=dge
}
### flip tSNE2
for(i in c(6,8)){
dge=dgelist[[i]]
dge@tsne.rot[,2]=-dge@tsne.rot[,2]
dgelist[[i]]=dge
}
### switch tSNE1 and tSNE2
for(i in c(2,5,6,8)){
dge=dgelist[[i]]
tmp=dge@tsne.rot[,2]
dge@tsne.rot[,2]=dge@tsne.rot[,1]
dge@tsne.rot[,1]=tmp
dgelist[[i]]=dge
}



## after reverse cluster order and flip cooridnates, repeat the above plots
# change ordered0 in the file name to order1
pdf("clusters_order1_PRM2_Violin.pdf",height=4,width=8)
for(i in 1:length(dataset)){
  VlnPlot(dgelist[[i]],"PRM2")
}
dev.off()
pdf("clusters_order1_ACTA2_Violin.pdf",height=4,width=8)
for(i in 1:length(dataset)){
  FeaturePlot(dgelist[[i]],"ACTA2")
}
dev.off()
pdf("clusters_order1_PRM2_Feature.pdf",height=4,width=4)
for(i in 1:length(dataset)){
  FeaturePlot(dgelist[[i]],"PRM2",col=c("grey80","red"))
}
dev.off()
pdf("clusters_order1_ACTA2_Feature.pdf",height=4,width=4)
for(i in 1:length(dataset)){
  FeaturePlot(dgelist[[i]],"ACTA2",col=c("grey80","red"))
}
dev.off()
#FeaturePlot(dgelist[[i]],"PRM1")
#TSNEPlot(dgelist[[i]],do.return=T)
# download the tSNE plot below to check the cluster IDs

## double-check if I need to flip PCs, or flip tSNE
### PCA and tSNE for ordered clusters of each individual replicate
plotlistt=plotlist2=plotlist3=plotlist4=plotlist5=list()
for(i in 1:length(dataset)){
dge=dgelist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist4[[i]]=PCAPlot(dge,1,4,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist5[[i]]=PCAPlot(dge,1,5,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE,label.size=4)
}
pdf("clusters_ordered1.pdf",height=8,width=10.5)
MultiPlotList(plotlist2,cols = 3)
MultiPlotList(plotlist3,cols = 3)
MultiPlotList(plotlist4,cols = 3)
MultiPlotList(plotlist5,cols = 3)
MultiPlotList(plotlistt,cols = 3)
dev.off()


for(i in 1:length(dataset)){
dge=dgelist[[i]]
save(dge,file=paste0(datainfo[i,1],".Robj"))
}

for(i in 1:length(dataset)){
dge=dgelist[[i]]
print(ncol(dge@data))
}
### Visualize PC1-3 and tSNE for each replicate of each subject
ncol=c(2,1,2,2,2)
nrow=c(1,1,1,1,1)
subject=unique(gsub("-.*","",dataset))
for(indiv in 1:length(subject)){
plotlist2=plotlist3=plotlistt=list()
for(i in 1:length(grep(subject[indiv],dataset))){
  rep=grep(subject[indiv],dataset)[i]
dge=dgelist[[rep]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE,label.size=4)
}
pdf(paste0(subject[indiv],"_clusters_order1.pdf"),height=3*nrow[indiv],width=3.8*ncol[indiv])
MultiPlotList(plotlist2,cols = ncol[indiv])
MultiPlotList(plotlist3,cols = ncol[indiv])
MultiPlotList(plotlistt,cols = ncol[indiv])
dev.off()
}


# 12.17.2018
# merge all replicates for each individual subject
dgealllist=list()

## Monkey1-5 individual replicate merged data
for(indiv in 1:length(subject)){

### load whole gene exp for individual replicate
setlist=list()
for(i in 1:length(grep(subject[indiv],dataset))){
rep=grep(subject[indiv],dataset)[i]
### read in data
dgedata=read.table(paste0(datainfo[rep,4],".dge.txt.gz"),header=T,row.names=1)
colnames(dgedata)=paste(dataset[rep],colnames(dgedata),sep="_")
### Filter for cells (>500 genes and <10% MT)
nGeneperCell <- colSums(dgedata>0)
dgedata.tmp=dgedata[,nGeneperCell>500]
MT <- c("ENSMMUG00000028704", "ENSMMUG00000028703", "ENSMMUG00000028702", "ENSMMUG00000028701", "ENSMMUG00000028700", "ENSMMUG00000028699", "ENSMMUG00000028698", "ENSMMUG00000028697", "ENSMMUG00000028696", "ENSMMUG00000028694", "ENSMMUG00000028693", "ENSMMUG00000028692", "ENSMMUG00000028691", "ENSMMUG00000028690", "ENSMMUG00000028688", "ENSMMUG00000028687", "ENSMMUG00000028686", "ENSMMUG00000028685", "ENSMMUG00000028684", "ENSMMUG00000028683", "ENSMMUG00000028681", "ENSMMUG00000028679", "ENSMMUG00000028678", "ENSMMUG00000028676", "ENSMMUG00000028675", "ENSMMUG00000028674", "ENSMMUG00000028671", "ENSMMUG00000028670", "ENSMMUG00000028669", "ENSMMUG00000028668", "ENSMMUG00000028695", "MT-CO1", "MT-CO3", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6")
MT=unique(MT)
mito.genes=MT[which(MT %in% rownames(dgedata.tmp))]
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1] 
print(c(ncol(dgedata),ncol(dgedata[,nGeneperCell>500]),ncol(dgedata.tmp)))
### Filter for genes - all detected genes
nCellperGene <- rowSums(dgedata.tmp>0)
nUMIperGene <- rowSums(dgedata.tmp)
dgedata2=dgedata.tmp[which(nUMIperGene>0),]
print(c(nrow(dgedata),nrow(dgedata2)))
print(summary(rowSums(dgedata2)))
setlist[[i]]=dgedata2
}

### merge all replicates for each individual
set=setlist
for(i in 1:length(grep(subject[indiv],dataset))){
  set[[i]]=data.frame(GENE=rownames(setlist[[i]]),setlist[[i]])
}
dgedataall=Reduce(function(x,y) merge(x,y,all=TRUE), set)
dgedataall[is.na(dgedataall)] <- 0
row.names(dgedataall)=dgedataall[,1]
dgedataall=dgedataall[,-1]
dim(dgedataall)             
nCellperGene <- rowSums(dgedataall>0)
length(which(nCellperGene==0))
nCellperGene[which(nCellperGene==0)]
dge <- new("seurat", raw.data = dgedataall)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = subject[indiv],names.field = 1,names.delim = "_")
dge                 
pdf(paste(subject[indiv],"_VariableGenes0.2.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.2,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"Variable Genes"))
dev.off()
length(dge@var.genes) 
dgeall=dge
dgealllist[[indiv]]=dgeall
save(dgeall,file=paste0(subject[indiv],".Robj"))
print(table(gsub("_.*","",names(dge@ident))))
}

### For subjects with too few HVG, change the criteria of selecting HVG
for(i in 1:length(subject)){
dge=dgealllist[[i]]
if(i==2 | i==5){
pdf(paste(subject[indiv],"_VariableGenes0.1.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.1,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"Variable Genes"))
dev.off()
}
if(i==4){
pdf(paste(subject[indiv],"_VariableGenes0.06.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.06,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"Variable Genes"))
dev.off()
}
dgealllist[[i]]=dge
print(length(dge@var.genes))
}


# cross-tabulation for all replicates for each individual

### load individual replicate clustering
dgelist=list()
for(i in 1:length(dataset)){
  load(file=paste0(datainfo[i,1],".Robj"))
  dgelist[[i]]=dge
}

### load object for each individual 
dgealllist=list()
for(indiv in 1:length(subject)){
  load(file=paste0(subject[indiv],".Robj"))
  dgealllist[[indiv]]=dgeall
}

## Monkey1-5 rank cor
datalist=list()
ExpMean=expMean

for(indiv in 1:length(subject)){
### order cells by batch first, then by clusters of each batch
blockident=NULL
for(i in 1:length(grep(subject[indiv],dataset))){
  rep=grep(subject[indiv],dataset)[i]
  tmp=paste(dataset[rep],dgelist[[rep]]@ident,sep="_")
  names(tmp)=paste(dataset[rep],names(dgelist[[rep]]@ident),sep="_")
  blockident=c(blockident,tmp)
}

### Clusters ordered first by batches, then by res
batch=grep(subject[indiv],dataset,value=T)
nbatch=length(batch)
ncluster=NULL
for(i in 1:length(grep(subject[indiv],dataset))){
  rep=grep(subject[indiv],dataset)[i]
  ncluster=c(ncluster,length(unique(dgelist[[rep]]@ident)))
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

names(blockident)=gsub("-",".",names(blockident))
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
dgeall=dgealllist[[indiv]]
dge=dgeall
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
datalist[[rep]]=data.use
write.table(data.use,paste0(subject[indiv],"_rep_Centroid_rho_HVG.txt"),row.names=T,col.names=T,quote=F,sep="\t")
### note: once calculated, next time can directly read table without calculating again

### load cluster centroid rank correlation using HVG
data.use=read.table(paste0(subject[indiv],"_rep_Centroid_rho_HVG.txt"),header=T,row.names=1)
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
pdf(file=paste(subject[indiv],"_rep_Centroid_RankedCorrelation_HVG.pdf",sep=""),height=5.5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=.8,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
}
# Concluded No between-replicate batch effect within each subject