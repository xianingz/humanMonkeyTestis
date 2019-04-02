### Global Clustering for Directly-merged all 5 monkey subjects
# by Qianyi on 12.18.2018
### visualized each individual batch in global PCs and t-SNE, Examined batch effect and markers

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


# 12.18.2018 by Qianyi 

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

### PCA Plot
pdf("PCA_Merged5.pdf",width=12,height=10)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
### Scree Plot for PCA
numPCs=8;i=1 # HVG
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
Sys.time() # [1] "2018-12-18 16:23:23 EST"
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,3,by=0.1), save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])
Sys.time() # [1] "2018-12-13 16:07:09 EST"
dgeall=dge
save(dgeall,file="MonkeyMerged5.Robj")


# 12.20.2018 
### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
print(c( length(unique(dge@meta.data$res.1.3)),length(unique(dge@meta.data$res.1.4)),length(unique(dge@meta.data$res.1.5)),length(unique(dge@meta.data$res.1.6)),length(unique(dge@meta.data$res.1.7)),length(unique(dge@meta.data$res.1.8)),length(unique(dge@meta.data$res.1.9)),length(unique(dge@meta.data$res.2)),length(unique(dge@meta.data$res.2.1)),length(unique(dge@meta.data$res.2.2)),length(unique(dge@meta.data$res.2.3)),length(unique(dge@meta.data$res.2.4)),length(unique(dge@meta.data$res.2.5)),length(unique(dge@meta.data$res.2.6)),length(unique(dge@meta.data$res.2.7)),length(unique(dge@meta.data$res.2.8)),length(unique(dge@meta.data$res.2.9)),length(unique(dge@meta.data$res.3)) ))
### decide to use 12 clusters, double-check the number of clusters
res="res.0.5";j=1;resi=1;
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
pdf(file=paste0("Merged5Subjects_Centroid_norm_Seriation_",res,".pdf"))
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
dge5=dge
save(dge,file="MonkeyMerged5.Robj")


## double-check if I need to reverse the cluster ID orders
### use PRM2 as marker for the last cluster (elongating)
pdf("clusters_ordered0_PRM2.pdf",height=4,width=8)
VlnPlot(dge,"NR2F2",cols.use=myBrewerPalette)
VlnPlot(dge,"ID4",cols.use=myBrewerPalette)
VlnPlot(dge,"PRM2",cols.use=myBrewerPalette)
VlnPlot(dge,"ACRV1",cols.use=myBrewerPalette)
FeaturePlot(dge,"PRM2")
dev.off()
dge=SetAllIdent(dge,id="res.0.5order")
pdf("clusters_ordered1_PRM2.pdf",height=10,width=10)
VlnPlot(dge,c("ID4","PIWIL4","TCF3","FGFR3","MORC1","ZCWPW1","ACRV1","PRM2"),cols.use=myBrewerPalette,nCol=2,point.size.use=-1)
dev.off()
pdf("clusters_ordered1_BEND2.pdf",height=10,width=5)
VlnPlot(dge,c("BEND2","ZCWPW1","SSX1","PRDM9"),cols.use=myBrewerPalette,nCol=1,point.size.use=-1)
dev.off()


## after reverse cluster order, repeat the above plots
# change ordered0 in the file name to order1
knownmarkers=c("VIM","CD163","S100A4","TYROBP","LYZ","RGS1","VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN","MYH11","ACTA2","DCN","DLK1","NR2F2","CYP17A1","CLU","ID4","FGFR3", "MORC1", "PIWIL4", "TCF3","ZCWPW1","PRM2","ACRV1")
length(knownmarkers) # 26
pdf("knownmarkers.pdf",height=8,width=21)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
knownmarkers=c("ID4","PIWIL4","TCF3","FGFR3","MORC1","BEND2","ZCWPW1","SSX1","PRDM9","ACRV1","PRM2")
length(knownmarkers) # 11
pdf("knownmarkers2.pdf",height=8,width=12)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("grey80","red"),nCol=4)
dev.off()


###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Merged5subjects_res.0.5order_markersall_mindiff0.2_logfc2fold_12.14.2018.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)

1    2    3    4    5    6    7    8    9   10   11   12   13 
496  714 1926  620  178  776  833  551 1083 2099 1827 1747 1208 
253 145 274 145 148 413 238 441 318 215 216 290 258

markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf("markerstop.pdf",height=12,width=12)
FeaturePlot(object = dge, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=4)
dev.off()

### plot for individual plots without Seurat package
redblue100.alpha<-rgb(read.table("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt",sep='\t',row.names=1,header=T),alpha=0.8)
set.ifnull=function(x,y) {
  if(is.null(x)) return(y)
  return(x)
}
translate.dim.code=function(reduction.use) {
  if (reduction.use=="pca") return.code="PC"
  if (reduction.use=="ica") return.code="IC"
  if (reduction.use=="tsne") return.code="tSNE_"
  if (reduction.use=="mds") return.code="MDS"
  if (reduction.use=="umap") return.code="UMAP"
  return(return.code)
}

object=dge
dgename="markers/"

setname="knownmarkers"
genes=knownmarkers

setname="markers"
genes=markers

for(j in 1:length(genes)){
feature=features.plot=genes[j]

features.plot; dim.1 = 1; dim.2 = 2; cells.use = NULL; pt.size = 1;
pch.use = 16; reduction.use = "tsne";
use.imputed = FALSE; no.axes = TRUE; no.legend = FALSE

cells.use <- set.ifnull(cells.use, colnames(object@data))
dim.code <- translate.dim.code(reduction.use)
dim.codes <- paste(dim.code, c(dim.1, dim.2), sep = "")
data.plot <- as.data.frame(object@dr$tsne@cell.embeddings)

x1 <- paste(dim.code, dim.1, sep = "")
x2 <- paste(dim.code, dim.2, sep = "")

data.plot$x <- data.plot[, x1]
data.plot$y <- data.plot[, x2]
data.plot$pt.size <- pt.size
data.use <- data.frame(t(object@data[feature,]))
rownames(data.use)=feature

data.gene0 <- na.omit(data.frame(data.use[feature, ]))
data.plot$gene <- t(data.gene0)

st6<- data.plot
st6<-st6[order(st6[,6]),]
z<-st6[,6]

jpeg(paste0(dgename,setname,"_",feature,"_40-100redblue0.8.jpeg"),height=850,width=800,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- redblue100.alpha[40:100][(z - min(z))/diff(range(z))*60 + 1]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()

jpeg(paste0(dgename,setname,"_",feature,"_redblue0.8.jpeg"),height=850,width=800,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()

}


### PCA and tSNE plot
dge=SetAllIdent(dge,id="orig.ident")
pdf("PCA_tSNE_UMAP_Merged5Monkey_Rep.pdf",width=10.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot3=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot4=PCAPlot(dge,1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot5=PCAPlot(dge,1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
dge=SetAllIdent(dge,id="indiv")
pdf("PCA_tSNE_UMAP_Merged5Monkey_Subject.pdf",width=10,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,6,do.return = TRUE,pt.size = .8,do.label=F)
plot4=PCAPlot(dge,1,7,do.return = TRUE,pt.size = .8,do.label=F)
plot5=PCAPlot(dge,1,8,do.return = TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
dge=SetAllIdent(dge,id="res.0.5order")
pdf("PCA_tSNE_UMAP_Merged5Monkey_res.0.5order.pdf",width=9,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=PCAPlot(dge,1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=PCAPlot(dge,1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

###### Visualize individual batch and subject
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlims=ylims=list()
pos=c("bottomright","topright","bottomleft","topleft","bottomright")
# dge5 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-56,16),c(-56,16),c(-56,16),c(-48,48),c(-13,9))
ylims[[1]]=list(c(-39,25),c(-18,36),c(-23,16),c(-45,48),c(-12,14))

dim=list(dge@dr$pca@cell.embeddings[,1:2],dge@dr$pca@cell.embeddings[,c(1,3)],dge@dr$pca@cell.embeddings[,c(1,4)],dge@dr$tsne@cell.embeddings[,1:2],dge@dr$umap@cell.embeddings[,1:2])

### plot PCs and tSNE for each batch using the other batches as background
dge=SetAllIdent(dge,id="orig.ident")
sets=levels(dge@meta.data$orig.ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette[c(8:7,10,6:5,4:1)]

pdf("OrigSet_PCtSNE.pdf",height=4.6,width=11.5)
par(mfrow=c(2,5),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:length(dim)){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(dge@ident)
names(ident)=names(dge@ident)
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.4,col=c("grey70",alpha(cols[seti],0.8))[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
plot.new()
}
dev.off()

### plot PCs and tSNE for each subject using the other subjects as background
dge=SetAllIdent(dge,id="indiv")
sets=levels(dge@meta.data$indiv)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

pdf("OrigSubject_PCtSNE.pdf",height=4.6,width=6.9)
par(mfrow=c(2,3),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:length(dim)){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(dge@ident)
names(ident)=names(dge@ident)
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.4,col=c("grey70",alpha(cols[seti],0.8))[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
plot.new()
}
dev.off()

### plot PCs and tSNE for each batch using the other batches as background
### visualize 12 clusters
dge=SetAllIdent(dge,id="res.0.5order")
table(dge@meta.data$indiv,dge@meta.data$res.0.5order)

  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@meta.data$orig.ident,dge@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,"ncellspercluster_batch_top8PCs.txt",quote=F,row.names=T,col.names=T,sep="\t")


dge=SetAllIdent(dge,id="res.0.5order")
sets=levels(dge@meta.data$orig.ident)
which(rownames(data)!=names(dge@ident))


cols=myBrewerPalette
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())

j=1
label="PC1-2"
data=data_bg=as.data.frame(dge@dr$pca@cell.embeddings[,1:2])
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(PC1 = median(PC1), PC2 = median(PC2)) -> centers
label.size=6

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(dge@ident)[which(gsub("_.*","",names(dge@ident))==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=PC1,y=PC2,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=cols,limits=levels(data$ident),drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}



j=4
label="tSNE"
data=data_bg=as.data.frame(dge@dr$tsne@cell.embeddings)
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) -> centers
label.size=6

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(dge@ident)[which(gsub("_.*","",names(dge@ident))==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=myBrewerPalette,limits=levels(data$ident),drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}

### plot
source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/Rcode_multiplot.R")

jpeg(file=paste0("OrigSet_12clusters_",label,"_bg.jpeg"),res=300,height=1200,width=3500)
multiplot(plotset,cols = 5)
dev.off()
pdf(file=paste0("OrigSet_12clusters_",label,"_bg.pdf"),height=6,width=15)
multiplot(plotset,cols = 5)
dev.off()

### end



###### Rank correlation for each normalized centroid using HVG
### order cell by cluster ID and randomly shuffle cells within each batch
levels=levels(dge@ident)
ident=factor(dge@ident,levels=levels)
### randomly shuffling cells within each cluster
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
levels(cells.ident)=levels### for each cluster, calculate average normalized expression of each gene
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

write.table(genecountsall,"Merged5subjects_12clusters_Centroid.txt",quote=F,sep="\t",row.names=T,col.names=T)

cc=cor(as.matrix(genecountsall)[dge@var.genes,],method="spearman")
dim(cc)
min(cc) # 0.036

### labeling
data.use=cc[levels,levels]
colsep.use=cumsum(table(gsub("_.*","",levels)))
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels)))-table(gsub("_.*","",levels))/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)

#ncluster=9
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep("white",length(sidecol[1,]))
sidecol[2,]=myBrewerPalette[as.numeric(gsub(".*-","",gsub(".*_","",levels)))]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

col.use=redblue100

pdf(file=paste(dgename,"5subjects_12clusters_Centroid_RankedCorrelation_HVG.pdf",sep=""),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()



### Per-cell attributes
###### per-cell attributes statistics
library(RColorBrewer)
myBrewerPalette <- brewer.pal(12,"Paired")
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
# used this color scheme for 13 clusters

pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()

###### PCA heatmap and tSNE heatmap
##### Plot nGenes, %MT,%x, %y, average cell cycle gene expression in PCA space
redblue100.alpha<-rgb(read.table("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt",sep='\t',row.names=1,header=T),alpha=0.8)

### 1. plot all together
set.ifnull=function(x,y) {
  if(is.null(x)) return(y)
  return(x)
}
translate.dim.code=function(reduction.use) {
  if (reduction.use=="pca") return.code="PC"
  if (reduction.use=="ica") return.code="IC"
  if (reduction.use=="tsne") return.code="tSNE_"
  if (reduction.use=="mds") return.code="MDS"
  if (reduction.use=="umap") return.code="UMAP"
  return(return.code)
}

dims=c("pca","tsne","umap")

### In all cells
dgefile=""
labels=c("nGene","nUMI","%MT")
featurelist=list(dge@meta.data$nGene,dge@meta.data$nUMI,dge@meta.data$percent.mito)

pdf(paste0(dgefile,"PerCellAttributes_heatmap_redblue0.8.pdf"),height=6,width=9)
par(mfrow=c(2,3),mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0),bg="grey")
for(i in 1:length(labels)){
for(dim in dims){
setname=labels[i]
feature=features.plot=featurelist[[i]]
names(feature)=rownames(dge@meta.data)

object=dge; dim.1 = 1; dim.2 = 2; cells.use = NULL; pt.size = 1;
                        pch.use = 16; reduction.use = dim;
                        use.imputed = FALSE; no.axes = TRUE; no.legend = FALSE


            cells.use <- set.ifnull(cells.use, colnames(object@data))
            dim.code <- translate.dim.code(reduction.use)
            dim.codes <- paste(dim.code, c(dim.1, dim.2), sep = "")
            data.plot <- FetchData(object, dim.codes, cells.use = cells.use)

            x1 <- paste(dim.code, dim.1, sep = "")
            x2 <- paste(dim.code, dim.2, sep = "")

            data.plot$x <- data.plot[, x1]
            data.plot$y <- data.plot[, x2]
            data.plot$pt.size <- pt.size
            #data.use <- data.frame(t(FetchData(object, features.plot, cells.use = cells.use,
            #                                   use.imputed = use.imputed)))


   data.plot$gene <- feature

st6<- data.plot
st6<-st6[order(st6[,6]),]
z<-st6[,6]

# redblue100 original
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]

if(labels[1]=="nGene"){
  if(i<3){
    z=log10(st6[,6])
    zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
#    zcolor <- redblue100.alpha[c(1:5,11:15,21:25,31:35,51:101)][(z - min(z))/diff(range(z))*70 + 1]
  } else if(i==3){
    zcolor <- redblue100.alpha[c(1:60,71,81,91,101)][(z - min(z))/diff(range(z))*63 + 1]
  } else if(i==7 | i==8){
    zcolor <- redblue100.alpha[c(1:5,51:101)][(z - min(z))/diff(range(z))*55 + 1]
  }
}


plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=paste0(setname,":",round(min(feature),2),"-",round(max(feature),2)),xlab="",ylab="")
}
}
dev.off()


# 1.22.2019 Qianyi
# label cell types
### label major cell types
1-2: Somatic - need to remove subcluster of somatic-Scyte/spermatid doublets
3-4: SPG 
5-7: spermatocytes
8-9: round spermatids
10-12: elongating spermatids


somatic=1:2
spg=3:4
scyte=5:7
rs=8:9
es=10:12

id=as.numeric(dge@ident)
names(id)=names(dge@ident)
celltype2=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
id[which(id %in% somatic)]<-"Somatic"
id[which(id %in% spg)]<-"Spermatogonia"
id[which(id %in% scyte)]<-"Spermatocyte"
id[which(id %in% rs)]<-"RoundSpermatid"
id[which(id %in% es)]<-"Elongating"

id=factor(id,levels=celltype2,order=T)
table(id)
      Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
          2562           1309           4015           2751          11603 

table(dge@ident)
   1    2    3    4    5    6    7    8    9   10   11   12 
1901  661  552  757 1531 1469 1015 1258 1493 2431 4593 4579 

id=factor(id,levels=celltype2,order=T)
dge=AddMetaData(dge,id,"CellType2")
dge=SetAllIdent(dge,id="CellType2")
dge@ident=id
table(dge@ident)
   Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
          2562           1309           4015           2751          11603 

### plot major cell types (not removing somatic doublet subcluster)

# note: need to flip PC2 and PC1, and tSNE1
# flip PCs


library(RColorBrewer)
myBrewerPalette <- c("gray60",brewer.pal(12,"Paired")[1:4]) # 1 
# used this on 12/20/2017
dge=SetAllIdent(dge,id="CellType2")
dge@ident=factor(dge@ident,levels=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"))

pdf("PCA_tSNE_UMAP_Merged5_5CellTypeswithMergedSomaticArm.pdf",width=9.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=PCAPlot(dge,1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=PCAPlot(dge,1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

save(dge,file="MonkeyMerged5.Robj")
dgeall=dge

## label somatic cell types
### remove the Cluster 1 of somatic-Scyte/spermatid doublet cells
somatic1=read.table("MonkeyMerged5Somatic-NoC1_6celltypes_label.txt",row.names=1)
somatic=as.character(somatic1$V2)
names(somatic)=rownames(somatic1)
tmp=names(dge@ident)[which(dge@ident=="Somatic")]
C1cells=tmp[which(!(tmp %in% names(somatic)))]
length(C1cells) # 207
NoC1=names(dge@ident)[which(!(names(dge@ident) %in% C1cells))]

dge=SubsetData(dge,cells.use=NoC1)
dge # 23101 genes across 22033 samples.
table(dge@ident)
Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
          2355           1309           4015           2751          11603 

### label somatic cell types
dge=SetAllIdent(dge,id="CellType2")
dge@ident=dge@meta.data$CellType2
names(dge@ident)=rownames(dge@meta.data)
germ=as.character(dge@ident)[which(dge@ident!="Somatic")]
names(germ)=names(dge@ident)[which(dge@ident!="Somatic")]
celltype=c(somatic,germ)
levels=c("Macrophage","Endothelial","Leydig","Pericyte","Myoid","IntProgenitor","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
celltype=factor(celltype,levels=levels,order=T)
table(celltype)
  Macrophage    Endothelial         Leydig       Pericyte          Myoid 
            37            191             48            186            604 
 IntProgenitor  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
          1289           1309           4015           2751          11603 

dge=AddMetaData(dge,celltype,"CellType")
dge=SetAllIdent(dge,id="CellType")
dge@ident=celltype
table(dge@ident)
Macrophage    Endothelial         Leydig       Pericyte          Myoid 
            37            191             48            186            604 
 IntProgenitor  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
          1289           1309           4015           2751          11603 

save(dge,file="MonkeyMerged5-NoSomaticC1.Robj")
write.table(dge@ident,"MonkeyMerged5-NoSomaticC1_10celltypes_label.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(dge@dr$pca@cell.embeddings[,1:3],"MonkeyMerged5-NoSomaticC1_PC1-3.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(dge@dr$tsne@cell.embeddings,"MonkeyMerged5-NoSomaticC1_tSNE.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(dge@dr$umap@cell.embeddings,"MonkeyMerged5-NoSomaticC1_UMAP.txt",row.names=T,col.names=F,quote=F,sep="\t")


### plot major cell types (removed somatic subcluster C1 of doublets)
library(RColorBrewer)
myBrewerPalette <- c("gray60",brewer.pal(12,"Paired")[1:4]) # 1 
# used this on 12/20/2017

dge=SetAllIdent(dge,id="CellType2")
dge@ident=dge@meta.data$CellType2
names(dge@ident)=rownames(dge@meta.data)

pdf("PCA_tSNE_UMAP_Merged5-NoSomaticC1_5CellTypeswithMergedSomaticArm.pdf",width=11.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=PCAPlot(dge,1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=PCAPlot(dge,1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()


### plot 10 cell types (removed C1 doublets)
library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(11,10,9,8,6,7,1:4)]  # 1 

dge=SetAllIdent(dge,id="CellType")
dge@ident=dge@meta.data$CellType
names(dge@ident)=rownames(dge@meta.data)

pdf("PCA_tSNE_UMAP_Merged5-NoSomaticC1_10CellTypes.pdf",width=11.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=PCAPlot(dge,1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=PCAPlot(dge,1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()



# 2.25.2019 Qianyi
## Visualize C1/8 and C7/7 Doublets on Global PCA and t-SNE
dge=dgeall
### C1/8 cells
somatic1=read.table("MonkeyMerged5Somatic-NoC1_6celltypes_label.txt",row.names=1)
somatic=as.character(somatic1$V2)
names(somatic)=rownames(somatic1)
tmp=names(dge@ident)[which(dge@ident=="Somatic")]
C1cells=tmp[which(!(tmp %in% names(somatic)))]
length(C1cells) # 207

### C7/7 cells
somatic2=read.table("MonkeyMerged5Somatic-NoC1_CCA_7ClustersID.txt",row.names=1)
C7cells=rownames(somatic2)[which(somatic2$V2==7)]
length(C7cells) # 18

### label C1/8 cells and C7/7 cells
ident=as.character(dge@ident)
names(ident)=names(dge@ident)
levels=levels(dge@ident)
ident[C1cells]<-"C1/8"
ident[C7cells]<-"C7/7"
ident=factor(ident,levels=c("C1/8","C7/7",levels))
table(dge@ident)
table(ident)
       Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
          2562           1309           4015           2751          11603 

          C1/8           C7/7        Somatic  Spermatogonia   Spermatocyte 
           207             18           2337           1309           4015 
RoundSpermatid     Elongating 
          2751          11603 

dge=AddMetaData(dge,ident,"celltypeDoublets")
dge@ident=ident
table(dge@ident,dge@meta.data$indiv)
                 Monkey1 Monkey2 Monkey3 Monkey4 Monkey5
  C1/8                19      21      83      62      22
  C7/7                 4       1       9       3       1
  
library(RColorBrewer)
myBrewerPalette <- c("red","magenta","gray60",brewer.pal(12,"Paired")[1:4]) # 1 
pdf("PCA_tSNE_UMAP_Merged5_Doublets.pdf",width=11,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

save(dge,file="MonkeyMerged5.Robj")
dgeall=dge


