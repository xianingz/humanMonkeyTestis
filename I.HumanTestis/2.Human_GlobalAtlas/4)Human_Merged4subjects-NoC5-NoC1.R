### Global Clustering for Directly-merged 4 subjects (Huamn1,2,3,5) After Removing Doublets Clusters 5 and 1
# by Qianyi on 3.4.2019
### merged 4 human subjects (1,2,3,5); did global clustering -> 13 global clusters
### removed cluster 5 of doublets from global clustering; re-do global clustering -> 13 global clusters
#Cluster 1-4: somatic cells
#Cluster 5-7: SPG cells
#Cluster 8-9: Spermatocyte 
#Cluster 10-11: round spermatids
#Cluster 12-13: elongating spermatids
### furthur removed cluster 1 of doublets from somatic subclustering; Re-did global clustering -> 15 global clusters
#Cluster 1-6: somatic cells
#Cluster 7-9: SPG cells
#Cluster 10-11: Spermatocyte 
#Cluster 12-13: round spermatids
#Cluster 14-15: elongating spermatids
### Incorporated 7 somatic cell types from somatic subclustering
#-> Final 11 major cell types for human testis

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


# 3.4.2019 by Qianyi 

### load object for merged 4 subjects after removing Cluster 5 and Cluster 1 of doublets
# PCA and clustering performed after removing Cluster 5 and before removing Cluster 1
# labeled cell types after removing cluster 1
load(file="HumanMerged4-1235-NoC5-NoSomaticC1celltype.Robj")
dge1=dge
dge # project HumanMerged4-1235-NoC5-NoC1: 46560 genes across 13837 samples

dge=SetAllIdent(dge,id="CellType2")
dge@ident=factor(dge@ident,levels=levels(dge@meta.data$CellType2))
table(dge@ident)
table(dge@ident,dge@meta.data$indiv)
#       Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
#          3722           2158           2200           3215           2542 
#                 Human1 Human2 Human3 Human5
#  Somatic           211    600    638   2273
#  Spermatogonia    1260    323    572      3
#  Spermatocyte      966    716    515      3
#  RoundSpermatid   1249   1432    534      0
#  Elongating       1234   1099    209      0

dim(dge@raw.data) # [1] 46560 13880
cells.use=names(dge@ident)
dgedata=dge@raw.data[,cells.use]
dim(dgedata)      # [1] 46560 13837
nCellperGene <- rowSums(dgedata>0)
length(which(nCellperGene==0))     # 19
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgedata2=dgedata[genes.use,]
dim(dgedata2)     # [1] 46541 13837

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project="HumanMerged4-1235-NoC5-NoC1", min.cells=1, min.genes=1)
mito.genes <- grep(pattern = "^MT-|^mt-", x = rownames(x = dge@data), value = TRUE)
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("Human_Merged4-1235-NoC5-NoC1_HVG0.2.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # 3469
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 46541 genes across 13837 samples.          
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
#Human1 Human2 Human3 Human5 
#  4920   4170   2468   2279 
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
dgeall=dge
save(dgeall,file="HumanMerged4-1235-NoC5-NoC1.Robj")
table(gsub("_.*","",names(dge@ident)))
#Human1.1 Human1.2 Human1.3 Human1.4 Human1.5 Human2.1 Human2.2 Human3.1 
#     888     1409      776      890      963     2198     1994      992 
#Human3.2   Human5 
#    1486     2284

### PCA
print(Sys.time())    # [1] "2019-03-04 12:00:54 EST"
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)
dgeall=dge
save(dge,file="HumanMerged4-1235-NoC5-NoC1.Robj")

### PCA Plot
pdf("PCA_Merged4-1235-NoC5-NoC1.pdf",width=12,height=10)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
### Scree Plot for PCA
numPCs=10;i=1 # HVG
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
Sys.time() # [1] "2019-03-04 12:19:21 EST"
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,3,by=0.1), save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])
dgeall=dge
save(dgeall,file="HumanMerged4-1235-NoC5-NoC1.Robj")
Sys.time() # [1] "2019-03-04 12:58:13 EST"


### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.6)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
print(c( length(unique(dge@meta.data$res.1.3)),length(unique(dge@meta.data$res.1.4)),length(unique(dge@meta.data$res.1.5)),length(unique(dge@meta.data$res.1.6)),length(unique(dge@meta.data$res.1.7)),length(unique(dge@meta.data$res.1.8)),length(unique(dge@meta.data$res.1.9)),length(unique(dge@meta.data$res.2)),length(unique(dge@meta.data$res.2.1)),length(unique(dge@meta.data$res.2.2)),length(unique(dge@meta.data$res.2.3)),length(unique(dge@meta.data$res.2.4)),length(unique(dge@meta.data$res.2.5)),length(unique(dge@meta.data$res.2.6)),length(unique(dge@meta.data$res.2.7)),length(unique(dge@meta.data$res.2.8)),length(unique(dge@meta.data$res.2.9)),length(unique(dge@meta.data$res.3)) ))

table(dge@meta.data$res.0.4,dge1@meta.data[rownames(dge@meta.data),]$CellType2)
     Somatic Spermatogonia Spermatocyte RoundSpermatid Elongating
  0        0             0         1717            566          0
  1        0             0            0              3       2074
  10     507             0            0              0          0
  11       0             1          477              0          0
  2        0             0            0           1386        468
  3        0          1605            4              0          0
  4     1464             0            0              0          0
  5        0             0            1           1260          0
  6      624             0            0              0          0
  7      612             0            0              0          0
  8        0           552            1              0          0
  9      515             0            0              0          0
table(dge@meta.data$res.0.5,dge1@meta.data$CellType2)
     Somatic Spermatogonia Spermatocyte RoundSpermatid Elongating
  0        0             0            0              2       2148
  1        0             0         1562            586          0
  10     515             0            0              0          0
  11     507             0            0              0          0
  12     380             0            0              0          0
  13     232             0            0              0          0
  2        0             0            0           1402        394
  3     1464             0            0              0          0
  4        0             0            1           1225          0
  5        0           805            4              0          0
  6        0           806            0              0          0
  7        0             1          632              0          0
  8      624             0            0              0          0
  9        0           546            1              0          0
table(dge@meta.data$res.0.6,dge1@meta.data$CellType2)
     Somatic Spermatogonia Spermatocyte RoundSpermatid Elongating
  0        0             0         1704             62          0
  1        0             0            6           1731          0
  10     514             0            0              0          0
  11     507             0            0              0          0
  12       0             0          481              0          0
  13     381             0            0              0          0
  14     232             0            0              0          0
  2     1451             0            0              0          0
  3        0             0            0           1417         31
  4        0             0            0              5       1325
  5        0             0            0              0       1186
  6        0          1023            0              0          0
  7      637             0            0              0          0
  8        0           588            4              0          0
  9        0           547            5              0          0


### decide to use 15 clusters, double-check the number of clusters
res="res.0.6";j=1;resi=1;
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
pdf(file=paste0("Merged4-1235-NoC5-NoC1_Centroid_norm_Seriation_",res,".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO"),col=bluered(10)))
dev.off()

### get order of seriation
do=seriate(da,method="OLO")
print(get_order(do))
levelss=get_order(do)
levelss=rev(levelss)
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
print(which(levels(cells)!=levelss[[j]]-1)) # integer(0)

ordered=paste0(res,"order")

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@meta.data[,ordered]=factor(dge@meta.data[,ordered])
dge=SetAllIdent(dge,ordered)
# save the dge file
dgeall=dge
save(dgeall,file="HumanMerged4-1235-NoC5-NoC1.Robj")
Sys.time() # [1] "2019-03-04 16:14:29 EST"




## double-check if I need to reverse the cluster ID orders
### use PRM1 as marker for the last cluster (elongating)
pdf("clusters_ordered0_PRM2.pdf",height=4,width=8)
VlnPlot(dge,"NR2F2",cols.use=myBrewerPalette)
VlnPlot(dge,"ID4",cols.use=myBrewerPalette)
VlnPlot(dge,"PRM2",cols.use=myBrewerPalette)
VlnPlot(dge,"ACRV1",cols.use=myBrewerPalette)
FeaturePlot(dge,"PRM2")
dev.off()
pdf("clusters_ordered1_PRM1.pdf",height=10,width=10)
VlnPlot(dge,c("ID4","PIWIL4","TCF3","FGFR3","MORC1","ZCWPW1","ACRV1","PRM1"),cols.use=myBrewerPalette,nCol=2,point.size.use=-1)
dev.off()
pdf("clusters_ordered1_SSX.pdf",height=10,width=5)
VlnPlot(dge,c("SSX1","SSX2","SSX3","ACRV1"),cols.use=myBrewerPalette,nCol=1,point.size.use=-1)
dev.off()
pdf("clusters_ordered1_BEND2.pdf",height=10,width=5)
VlnPlot(dge,c("BEND2","ZCWPW1","SSX3","PRDM9"),cols.use=myBrewerPalette,nCol=1,point.size.use=-1)
dev.off()

knownmarkers=c("ID4","PIWIL4","TCF3","FGFR3","MORC1","BEND2","ZCWPW1","SSX1","SSX2","SSX3","PRDM9","ACRV1","PRM1")

## after reverse cluster order, repeat the above plots
### compare with previous 13 clusters of merged 4 subjects after removing cluster 5, before removing cluster 1
table(dge@ident,dge1@meta.data[rownames(dge@meta.data),]$CellType2)
     Somatic Spermatogonia Spermatocyte RoundSpermatid Elongating
  1      232             0            0              0          0
  2      381             0            0              0          0
  3      514             0            0              0          0
  4     1451             0            0              0          0
  5      637             0            0              0          0
  6      507             0            0              0          0
  7        0          1023            0              0          0
  8        0           588            4              0          0
  9        0           547            5              0          0
  10       0             0          481              0          0
  11       0             0         1704             62          0
  12       0             0            6           1731          0
  13       0             0            0           1417         31
  14       0             0            0              5       1325
  15       0             0            0              0       1186
table(dge@ident,dge1@meta.data[rownames(dge@meta.data),]$res.0.6order)


# change ordered0 in the file name to order1
knownmarkers=c("VIM","CD163","S100A4","TYROBP","LYZ","RGS1","VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN","MYH11","ACTA2","PDGFRA","DCN","DLK1","NR2F2","CYP17A1","CLU","ID4","FGFR3", "MORC1", "PIWIL4", "TCF3","ZCWPW1","TSPAN33","PRM1","PRM2","ACRV1")
length(knownmarkers) # 30
pdf("knownmarkers.pdf",height=15,width=18)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("grey80","red"),nCol=6)
dev.off()


###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Merged4-1235-NoC5-NoC1_res.0.6order_markersall_mindiff0.2_logfc2fold_3.4.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)

1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
232  381  514 1451  637  507 1023  592  552  481 1766 1737 1448 1330 1186 
264 134 195 284 170 253 442 126 451 387 270 265 217 275 249 

markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Merged4-1235-NoC5-NoC1_11celltypes_markersall_mindiff0.2_logfc2fold_3.20.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
         Tcell     Macrophage    Endothelial     m-Pericyte     f-Pericyte 
            65            231            297            410             79 
           251            266            160            263            216 
         Myoid      ImmLeydig  Spermatogonia   Spermatocyte RoundSpermatid 
          2121            519           2167           2247           3185 
           251            202            521            274            280 
    Elongating 
          2516 
           351 



markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf("markerstop.pdf",height=15,width=18)
FeaturePlot(object = dge, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=6)
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

# edit gene name
#if(features.plot=="Trdmt1"){features.plot="Dnmt2"}
#if(features.plot=="Gm1564"){features.plot="Meioc"}

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
pdf("PCA_tSNE_UMAP_Merged4-1235-NoC5-NoC1_Rep.pdf",width=10.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot3=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot4=PCAPlot(dge,1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot5=PCAPlot(dge,1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
dge=SetAllIdent(dge,id="indiv")
pdf("PCA_tSNE_UMAP_Merged4-1235-NoC5-NoC1_Subject.pdf",width=10,height=8)
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
dge=SetAllIdent(dge,id="res.0.6order")
pdf("PCA_tSNE_UMAP_Merged4-1235-NoC5-NoC1_res.0.6order.pdf",width=9.5,height=8)
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

###### Visualize individual batch and subject
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlims=ylims=list()
pos=c("bottomleft","bottomright","bottomleft","bottomright")
# dge4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-19,39),c(-19,39),c(-19,11),c(-51,58))
ylims[[1]]=list(c(-34,17),c(-24,11),c(-23,12),c(-57,54))


### plot PCs and tSNE for each batch using the other batches as background
dge=SetAllIdent(dge,id="orig.ident")
sets=levels(dge@meta.data$orig.ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette[c(5:9,2:1,4:3,10)]

dim=list(dge@dr$pca@cell.embeddings[,1:2],dge@dr$pca@cell.embeddings[,c(1,3)],dge@dr$umap@cell.embeddings,dge@dr$tsne@cell.embeddings)
xlim=xlims[[i]]
ylim=ylims[[i]]

pdf("OrigSet_PCtSNE.pdf",height=4.6,width=11.5)
par(mfrow=c(2,5),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
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
}
dev.off()

### Save As Figure S1?

### plot PCs and tSNE for each subject using the other subjects as background
dge=SetAllIdent(dge,id="indiv")
sets=levels(dge@meta.data$indiv)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

dim=list(dge@dr$pca@cell.embeddings[,1:2],dge@dr$pca@cell.embeddings[,c(1,3)],dge@dr$pca@cell.embeddings[,c(1,4)],dge@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

pdf("OrigSubject_PCtSNE.pdf",height=4.6,width=4.6)
par(mfrow=c(2,2),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
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
}
dev.off()

### plot PCs and tSNE for each batch using the other batches as background
### visualize 12 clusters
dge=SetAllIdent(dge,id="res.0.6order")
table(dge@meta.data$indiv,dge@meta.data$res.0.6order)

  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@meta.data$orig.ident,dge@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,"ncellspercluster_batch_top10PCs.txt",quote=F,row.names=T,col.names=T,sep="\t")


dge=SetAllIdent(dge,id="res.0.6order")
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

write.table(genecountsall,"Merged4-1235-NoC5-NoC1_15clusters_Centroid.txt",quote=F,sep="\t",row.names=T,col.names=T)
write.table(genecountsall,"Merged4-1235-NoC5-NoC1_11celltypes_Centroid.txt",quote=F,sep="\t",row.names=T,col.names=T)

cc=cor(as.matrix(genecountsall)[dge@var.genes,],method="spearman")
dim(cc)
min(cc) # 0.036

### labeling
data.use=cc[levels,levels]
colsep.use=cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])
row.lab=col.lab=gsub(".*_","",levels)

#ncluster=9
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep("white",length(sidecol[1,]))
sidecol[2,]=myBrewerPalette[as.numeric(gsub(".*-","",gsub(".*_","",levels)))]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")
col.use=redblue100

pdf(file="4subjects_15clusters_Centroid_RankedCorrelation_HVG.pdf",height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()

### visualize for cell types
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep("white",length(sidecol[1,]))
sidecol[2,]=myBrewerPalette
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","CellType")
colnames(clab)=c("CellType","")

col.use=redblue100
pdf(file=paste0(dgefile,"4subjectsNoC5NoC1_11celltypes_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,7))
dev.off()
# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)-1):101]
pdf(file=paste0(dgefile,"4subjectsNoC5NoC1_11celltypes_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,7))
dev.off()



### Per-cell attributes
######### Analyze per-cell attributes 
###### add %ChrX and %ChrY genes to the dge datainfo
x=read.table("../humanChrXgenes",stringsAsFactors=FALSE)[,1]
y=read.table("../humanChrYgenes",stringsAsFactors=FALSE)[,1]
autosome=read.table("../humanautosomegenes",stringsAsFactors=FALSE)
length(x) # 2321
length(y) # 494
dim(autosome) # 50105
table(autosome[,2])
 HUMAN_1 HUMAN_10 HUMAN_11 HUMAN_12 HUMAN_13 HUMAN_14 HUMAN_15 HUMAN_16 
    5202     2185     3140     2738     1201     2177     2018     2274 
HUMAN_17 HUMAN_18 HUMAN_19  HUMAN_2  HUMAN_3  HUMAN_4  HUMAN_5  HUMAN_6 
    2830     1089     2876     3925     2990     2488     2786     2818 
 HUMAN_7  HUMAN_8  HUMAN_9 
    2773     2331     2264

chrx.genes <- x[which(x %in% rownames(dge@data))]
chry.genes <- y[which(y %in% rownames(dge@data))]
autosome.genes <- autosome[which(autosome[,1] %in% rownames(dge@data)),]
length(chrx.genes)  # 1707
length(chry.genes)  # 189
dim(autosome.genes) # 41889  2
table(autosome.genes[,2])

# %x expression
percent.x <- colSums(expm1(dge@data[chrx.genes, ]))/colSums(expm1(dge@data))
percent.y <- colSums(expm1(dge@data[chry.genes, ]))/colSums(expm1(dge@data))
percent.autosome <- colSums(expm1(dge@data[autosome.genes[,1], ]))/colSums(expm1(dge@data))
# named integer(0)

# add to dge
dge <- AddMetaData(dge, percent.x, "percent.x")
dge <- AddMetaData(dge, percent.y, "percent.y")
dge <- AddMetaData(dge, percent.autosome, "percent.autosome")

###### calculate Gini index
library(ineq)
#for(i in 1:(length(dataset)+1)){
#  dge=dgelist[[i]]
GiniAll=apply( dge@raw.data,2,function(x) ineq(x,type=c("Gini")) )
GiniNon0=apply( dge@raw.data,2,function(x) ineq(x[which(x!=0)],type=c("Gini")) )

dge=AddMetaData(dge,GiniAll,"GiniAll")
dge=AddMetaData(dge,GiniNon0,"GiniNon0")

###### Cell cycle index 
G1S=c("ACD","ACYP1","ADAMTS1","ANKRD10","APEX2","ARGLU1","ATAD2","BARD1","BRD7","C1orf63","C7orf41","C14orf142","CAPN7","CASP2","CASP8AP2","CCNE1","CCNE2","CDC6","CDC25A","CDCA7","CDCA7L","CEP57","CHAF1A","CHAF1B","CLSPN","CREBZF","CTSD","DIS3","DNAJC3","DONSON","DSCC1","DTL","E2F1","EIF2A","ESD","FAM105B","FAM122A","FLAD1","GINS2","GINS3","GMNN","HELLS","HOXB4","HRAS","HSF2","INSR","INTS8","IVNS1ABP","KIAA1147","KIAA1586","LNPEP","LUC7L3","MCM2","MCM4","MCM5","MCM6","MDM1","MED31","MRI1","MSH2","NASP","NEAT1","NKTR","NPAT","NUP43","ORC1","OSBPL6","PANK2","PCDH7","PCNA","PLCXD1","PMS1","PNN","POLD3","RAB23","RECQL4","RMI2","RNF113A","RNPC3","SEC62","SKP2","SLBP","SLC25A36","SNHG10","SRSF7","SSR3","TAF15","TIPIN","TOPBP1","TRA2A","TTC14","UBR7","UHRF1","UNG","USP53","VPS72","WDR76","ZMYND19","ZNF367","ZRANB2")                                                                                                          
S=c("ABCC5","ABHD10","ANKRD18A","ASF1B","ATAD2","BBS2","BIVM","BLM","BMI1","BRCA1","BRIP1","C5orf42","C11orf82","CALD1","CALM2","CASP2","CCDC14","CCDC84","CCDC150","CDC7","CDC45","CDCA5","CDKN2AIP","CENPM","CENPQ","CERS6","CHML","COQ9","CPNE8","CREBZF","CRLS1","DCAF16","DEPDC7","DHFR","DNA2","DNAJB4","DONSON","DSCC1","DYNC1LI2","E2F8","EIF4EBP2","ENOSF1","ESCO2","EXO1","EZH2","FAM178A","FANCA","FANCI","FEN1","GCLM","GOLGA8A","GOLGA8B","H1F0","HELLS","HIST1H2AC","HIST1H4C","INTS7","KAT2A","KAT2B","KDELC1","KIAA1598","LMO4","LYRM7","MAN1A2","MAP3K2","MASTL","MBD4","MCM8","MLF1IP","MYCBP2","NAB1","NEAT1","NFE2L2","NRD1","NSUN3","NT5DC1","NUP160","OGT","ORC3","OSGIN2","PHIP","PHTF1","PHTF2","PKMYT1","POLA1","PRIM1","PTAR1","RAD18","RAD51","RAD51AP1","RBBP8","REEP1","RFC2","RHOBTB3","RMI1","RPA2","RRM1","RRM2","RSRC2","SAP30BP","SLC38A2","SP1","SRSF5","SVIP","TOP2A","TTC31","TTLL7","TYMS","UBE2T","UBL3","USP1","ZBED5","ZWINT")                                                                             
G2M=c("ANLN","AP3D1","ARHGAP19","ARL4A","ARMC1","ASXL1","ATL2","AURKB","BCLAF1","BORA","BRD8","BUB3","C2orf69","C14orf80","CASP3","CBX5","CCDC107","CCNA2","CCNF","CDC16","CDC25C","CDCA2","CDCA3","CDCA8","CDK1","CDKN1B","CDKN2C","CDR2","CENPL","CEP350","CFD","CFLAR","CHEK2","CKAP2","CKAP2L","CYTH2","DCAF7","DHX8","DNAJB1","ENTPD5","ESPL1","FADD","FAM83D","FAN1","FANCD2","G2E3","GABPB1","GAS1","GAS2L3","H2AFX","HAUS8","HINT3","HIPK2","HJURP","HMGB2","HN1","HP1BP3","HRSP12","IFNAR1","IQGAP3","KATNA1","KCTD9","KDM4A","KIAA1524","KIF5B","KIF11","KIF20B","KIF22","KIF23","KIFC1","KLF6","KPNA2","LBR","LIX1L","LMNB1","MAD2L1","MALAT1","MELK","MGAT2","MID1","MIS18BP1","MND1","NCAPD3","NCAPH","NCOA5","NDC80","NEIL3","NFIC","NIPBL","NMB","NR3C1","NUCKS1","NUMA1","NUSAP1","PIF1","PKNOX1","POLQ","PPP1R2","PSMD11","PSRC1","RANGAP1","RCCD1","RDH11","RNF141","SAP30","SKA3","SMC4","STAT1","STIL","STK17B","SUCLG2","TFAP2A","TIMP1","TMEM99","TMPO","TNPO2","TOP2A","TRAIP","TRIM59","TRMT2A","TTF2","TUBA1A","TUBB","TUBB2A","TUBB4B","TUBD1","UACA","UBE2C","VPS25","VTA1","WSB1","ZNF587","ZNHIT2")                                        

M=c("AHI1","AKIRIN2","ANKRD40","ANLN","ANP32B","ANP32E","ARHGAP19","ARL6IP1","ASXL1","ATF7IP","AURKA","BIRC2","BIRC5","BUB1","CADM1","CCDC88A","CCDC90B","CCNA2","CCNB2","CDC20","CDC25B","CDC27","CDC42EP1","CDCA3","CENPA","CENPE","CENPF","CEP55","CFLAR","CIT","CKAP2","CKAP5","CKS1B","CKS2","CNOT10","CNTROB","CTCF","CTNNA1","CTNND1","DEPDC1","DEPDC1B","DIAPH3","DLGAP5","DNAJA1","DNAJB1","DR1","DZIP3","E2F5","ECT2","FAM64A","FOXM1","FYN","G2E3","GADD45A","GAS2L3","GOT1","GRK6","GTSE1","HCFC1","HMG20B","HMGB3","HMMR","HN1","HP1BP3","HPS4","HS2ST1","HSPA8","HSPA13","INADL","KIF2C","KIF5B","KIF14","KIF20B","KLF9","LBR","LMNA","MCM4","MDC1","MIS18BP1","MKI67","MLLT4","MZT1","NCAPD2","NCOA5","NEK2","NUF2","NUP35","NUP98","NUSAP1","ODF2","ORAOV1","PBK","PCF11","PLK1","POC1A","POM121","PPP1R10","PRPSAP1","PRR11","PSMG3","PTP4A1","PTPN9","PWP1","QRICH1","RAD51C","RANGAP1","RBM8A","RCAN1","RERE","RNF126","RNF141","RNPS1","RRP1","SEPHS1","SETD8","SFPQ","SGOL2","SHCBP1","SMARCB1","SMARCD1","SPAG5","SPTBN1","SRF","SRSF3","SS18","SUV420H1","TACC3","THRAP3","TLE3","TMEM138","TNPO1","TOMM34","TPX2","TRIP13","TSG101","TSN","TTK","TUBB4B","TXNDC9","TXNRD1","UBE2D3","USP13","USP16","VANGL1","WIBG","WSB1","YWHAH","ZC3HC1","ZFX","ZMYM1","ZNF207")   
MG1=c("AGFG1","AGPAT3","AKAP13","AMD1","ANP32E","ANTXR1","BAG3","BTBD3","CBX3","CDC42","CDK7","CDKN3","CEP70","CNIH4","CTR9","CWC15","DCP1A","DCTN6","DEXI","DKC1","DNAJB6","DSP","DYNLL1","EIF4E","ELP3","FAM60A","FAM189B","FOPNL","FOXK2","FXR1","G3BP1","GATA2","GNB1","GRPEL1","GSPT1","GTF3C4","HIF1A","HMG20B","HMGCR","HSD17B11","HSPA8","ILF2","JMJD1C","KDM5B","KIAA0586","KIF5B","KPNB1","KRAS","LARP1","LARP7","LRIF1","LYAR","MORF4L2","MRPL19","MRPS2","MRPS18B","MSL1","MTPN","NCOA3","NFIA","NFIC","NUCKS1","NUFIP2","NUP37","ODF2","OPN3","PAK1IP1","PBK","PCF11","PLIN3","PPP2CA","PPP2R2A","PPP6R3","PRC1","PSEN1","PTMS","PTTG1","RAD21","RAN","RHEB","RPL13A","SLC39A10","SNUPN","SRSF3","STAG1","SYNCRIP","TAF9","TCERG1","TLE3","TMEM138","TOB2","TOP1","TROAP","TSC22D1","TULP4","UBE2D3","VANGL1","VCL","WIPF2","WWC1","YY1","ZBTB7A","ZCCHC10","ZNF24","ZNF281","ZNF593")                                                                                             

print(c(length(G1S),length(S),length(G2M),length(M),length(MG1)))
#[1]  100 113 133 151 106
G1S=G1S[which(G1S %in% rownames(dge@data))]
S=S[which(S %in% rownames(dge@data))]
G2M=G2M[which(G2M %in% rownames(dge@data))]
M=M[which(M %in% rownames(dge@data))]
MG1=MG1[which(MG1 %in% rownames(dge@data))]
print(c(length(G1S),length(S),length(G2M),length(M),length(MG1)))
# [1]  99 112 133 150 106
labels=c("G1-S","S","G2-M","M","M-G1")
cellcycle=list(G1S,S,G2M,M,MG1)


#### Plot expression of cell cycle genes passing cor>0.3 with average expression pattern
cellcyclegene=cellcycle

### %expression of cell cycle genes in all cells 
dge@data=as.matrix(dge@data)
exp0=colSums(expm1(dge@data[unlist(cellcyclegene), ]))/colSums(expm1(dge@data))
exp1=colSums(expm1(dge@data[cellcyclegene[[1]], ]))/colSums(expm1(dge@data))
exp2=colSums(expm1(dge@data[cellcyclegene[[2]], ]))/colSums(expm1(dge@data))
exp3=colSums(expm1(dge@data[cellcyclegene[[3]], ]))/colSums(expm1(dge@data))
exp4=colSums(expm1(dge@data[cellcyclegene[[4]], ]))/colSums(expm1(dge@data))
exp5=colSums(expm1(dge@data[cellcyclegene[[5]], ]))/colSums(expm1(dge@data))

ccinall=cbind(exp0,exp1,exp2,exp3,exp4,exp5)
write.table(ccinall,"cellcycle_Merged4-1235-NoC5-NoC1.txt",col.names=T,row.names=T,sep="\t",quote=F)
dge=AddMetaData(dge,exp0,"FracAllCellCycle")
which(exp0!=dge@meta.data$FracAllCellCycle) # named integer(0)
dge3=dge


###### per-cell attributes statistics
pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.x", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.y", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.autosome", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="GiniAll", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="GiniNon0", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="FracAllCellCycle", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()



# 3.19.2019
# label cell types
### label major cell types
Cluster 1-6: somatic cells
Cluster 7-9: SPG cells
Cluster 10-11: Spermatocyte # 10 is pachytene 
Cluster 12-13: round spermatids
Cluster 14-15: elongating spermatids

## label major cell types
somatic=1:6
spg=7:9
scyte=10:11
rs=12:13
es=14:15

id=as.numeric(dge@ident)
names(id)=names(dge@ident)
celltype2=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
id[which(id %in% somatic)]<-"Somatic"
id[which(id %in% spg)]<-"Spermatogonia"
id[which(id %in% scyte)]<-"Spermatocyte"
id[which(id %in% rs)]<-"RoundSpermatid"
id[which(id %in% es)]<-"Elongating"

id=factor(id,levels=celltype2,order=T)
dge=AddMetaData(dge,id,"CellType2")
dge=SetAllIdent(dge,id="CellType2")
dge@ident=id
table(dge@ident)
Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
3722           2167           2247           3185           2516 

### plot major cell types (re-did clustering after removing C5 and C1 doublets)

# note: do not need to flip PC2 and PC1, and tSNE1
# flip PCs

library(RColorBrewer)
myBrewerPalette <- c("gray60",brewer.pal(12,"Paired")[1:4]) # 1 
# used this on 12/20/2017
dge=SetAllIdent(dge,id="CellType2")
dge@ident=factor(dge@ident,levels=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"))

pdf("PCA_tSNE_UMAP_Merged4-1235-NoC5-NoC1_5CellTypeswithMergedSomaticArm.pdf",width=11.3,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
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

save(dge,file="HumanMerged4-1235-NoC5-NoC1.Robj")

### visualize in different color schemes:

library(RColorBrewer)
col1 <- c("gray60",brewer.pal(4,"Paired"))  
col2 <- c("gray60",brewer.pal(4, "Set2"))  
col3 <- c("gray60",brewer.pal(4, "Set3")) 
library(viridis)
col4 <- c("gray60",viridis(4)) 
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
col5 <- c("gray60",gg_color_hue(4)) #  used this
library("colorspace") 
col6 <- c("gray60",rainbow_hcl(4)) 
col7 <- c("gray60",topo.colors(4)) 
myBrewerPalette=col5

library(RColorBrewer)
col11 <- c(brewer.pal(12,"Paired")[c(12,11,10,6,8,7,9)])  #used this 1 
col22 <- brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3)]
col33 <- c(brewer.pal(8,"Set1")[c(7,6,4)],gg_color_hue(4))
myBrewerPalette=col11
myBrewerPalette=col22


pdf("UMAP_col_Merged4-1235-NoC5-NoC1_5CellTypeswithMergedSomaticArm.pdf",width=11.3,height=8)
plot1=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col1)
plot2=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col2)
plot3=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col3)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col4)
plot_grid(plot1,plot2,plot3,plot4,ncol = 2)
plot1=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col5)
plot2=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col6)
plot3=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col7)
plot_grid(plot1,plot2,plot3,plot4,ncol = 2)
dev.off()


## label somatic cell types
somatic1=read.table("Human4-1235-NoC5SomaticNoC1_CCA_7celltypes_label.txt",row.names=1)
somatic=as.character(somatic1$V2)
names(somatic)=rownames(somatic1)

dge=SetAllIdent(dge,id="CellType2")
dge@ident=dge@meta.data$CellType2
names(dge@ident)=rownames(dge@meta.data)
germ=as.character(dge@ident)[which(dge@ident!="Somatic")]
names(germ)=names(dge@ident)[which(dge@ident!="Somatic")]
celltype=c(somatic,germ)
levels=c("Tcell","Macrophage","Endothelial","m-Pericyte","f-Pericyte","Myoid","ImmLeydig","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
celltype=factor(celltype,levels=levels,order=T)
table(celltype)
         Tcell     Macrophage    Endothelial     m-Pericyte     f-Pericyte 
            65            231            297            410             79 
         Myoid      ImmLeydig  Spermatogonia   Spermatocyte RoundSpermatid 
          2121            519           2167           2247           3185 
    Elongating 
          2516 


dge=AddMetaData(dge,celltype,"CellType")
dge=SetAllIdent(dge,id="CellType")
dge@ident=celltype
table(dge@ident)
Tcell     Macrophage    Endothelial     m-Pericyte     f-Pericyte 
65            231            297            410             79 
         Myoid      ImmLeydig  Spermatogonia   Spermatocyte RoundSpermatid 
          2121            519           2167           2247           3185 
    Elongating 
          2516 



### plot 11 major cell types (clustering after removing C5 and C1 doublets)

# note: do not need to flip UMAP2
#dge@dr$umap@cell.embeddings[,2]=-dge@dr$umap@cell.embeddings[,2]
# flip PCs
library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,8,7,9,1:4)]  #used this 1 

dge=SetAllIdent(dge,id="CellType")
dge@ident=dge@meta.data$CellType
names(dge@ident)=rownames(dge@meta.data)

pdf("PCA_tSNE_UMAP_Merged4-1235-NoC5-NoC1_11CellTypes.pdf",width=11.3,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
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

col11 <- c(brewer.pal(12,"Paired")[c(12,11,10,6,8,7,9,1:4)])  #used this 1 
col22 <- c(brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3)],gg_color_hue(4))
col33 <- c(brewer.pal(8,"Set1")[c(7,6,4)],gg_color_hue(4),brewer.pal(12,"Paired")[1:4])
col44 <- c(brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3)],rainbow_hcl(4))
myBrewerPalette=col22
pdf("UMAP_col_Merged4-1235-NoC5-NoC1_11CellTypes.pdf",width=11.3,height=8)
plot1=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot2=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot3=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col33)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col44)
plot_grid(plot1,plot2,plot3,plot4,ncol = 2)
dev.off()

save(dge,file="HumanMerged4-1235-NoC5-NoC1.Robj")
write.table(dge@ident,"HumanMerged4-1235-NoC5-NoC1_11celltypes_label.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(dge@dr$pca@cell.embeddings[,1:3],"HumanMerged4-1235-NoC5-NoC1_PC1-3.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(dge@dr$tsne@cell.embeddings,"HumanMerged4-1235-NoC5-NoC1_tSNE.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(dge@dr$umap@cell.embeddings,"HumanMerged4-1235-NoC5-NoC1_UMAP.txt",row.names=T,col.names=F,quote=F,sep="\t")

myBrewerPalette=col22
pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
VlnPlot(dge, features.plot="percent.x", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
VlnPlot(dge, features.plot="percent.y", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
VlnPlot(dge, features.plot="percent.autosome", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
VlnPlot(dge, features.plot="GiniAll", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
VlnPlot(dge, features.plot="GiniNon0", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))
dev.off()

pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot_ChrY.pdf"),height=4,width=8)
VlnPlot(dge, features.plot="percent.y", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))+coord_cartesian(ylim=c(0, 0.01))
VlnPlot(dge, features.plot="percent.autosome", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))+coord_cartesian(ylim=c(0.75, 1))
VlnPlot(dge, features.plot="percent.y", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)+theme(axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=45, size=12,hjust=1))+coord_cartesian(ylim=c(0, 0.0075))
dev.off()


###### Correlation between Gini and nUMI
datainfo=dge@meta.data
datainfo_bg=datainfo[,! names(datainfo) %in% "CellType"]

### Gini Vs nUMI
pdf(file=paste(dgefile,"GiniVsnUMI.pdf",sep=""),height=3.5,width=5.5)
ggplot(datainfo,aes(x=nUMI,y=GiniAll,color=CellType))+
geom_point(size=1) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
dev.off()

# Gini for all genes of each cell
pdf(file=paste(dgefile,"GiniVsnUMI2.pdf",sep=""),height=4.5,width=7)
ggplot(datainfo,aes(x=nUMI,y=GiniAll,color=CellType))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~CellType,ncol=4) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
# Gini for non-0 genes of each cell
ggplot(datainfo,aes(x=nUMI,y=GiniNon0,color=CellType))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~CellType,ncol=4) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
# Gini for highly-variable genes of each cell
ggplot(datainfo,aes(x=nUMI,y=GiniHVG,color=CellType))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~CellType,ncol=4) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
dev.off()

cor.test(datainfo$GiniAll,datainfo$nUMI)
for(i in 1:length(levels(datainfo$CellType))){
print(cor(datainfo$GiniAll[which(datainfo$CellType==levels(datainfo$CellType)[i])],datainfo$nUMI[which(datainfo$CellType==levels(datainfo$CellType)[i])]))
}
for(i in 1:length(levels(datainfo$CellType))){
print(cor.test(datainfo$GiniAll[which(datainfo$CellType==levels(datainfo$CellType)[i])],datainfo$nUMI[which(datainfo$CellType==levels(datainfo$CellType)[i])])$p.val)
}



######## Heatmap for all markers
dge=dgeall
celltype=dge@meta.data$CellType
names(celltype)=rownames(dge@meta.data)

dge=SetAllIdent(dge,id="CellType")
dge@ident=dge@meta.data$CellType
names(dge@ident)=rownames(dge@meta.data)

### centroids
centroid=matrix(,dim(dge@data)[1],length(levels(dge@ident)))
rownames(centroid)=rownames(dge@data)
colnames(centroid)=levels(dge@ident)
for(i in levels(dge@ident)){
    centroid[,i]=apply(dge@data[,which(dge@ident==i)],1,function(x) ExpMean(as.numeric(x))) 
    print(i)
}
centroid2=AverageExpression(dge)
which(round(centroid,2) != round(log(centroid2+1),2)) # integer(0)
write.table(centroid,"HumanMerged4-1235-NoC5-NoC1_11celltypes_centroids_allgenes.txt",row.names=T,col.names=T,quote=F,sep="\t")

### Genes Standardized Across Cell Types
# note: used this
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
write.table(centroid.std,"HumanMerged4-1235-NoC5-NoC1_11celltypes_centroid_allgenes_std.txt",row.names=T,col.names=T,quote=F,sep="\t")

centroid=read.table(paste0("HumanMerged4-1235-NoC5-NoC1_11celltypes_centroids_allgenes.txt"),header=T,row.names=1)
centroid.std=read.table(paste0("HumanMerged4-1235-NoC5-NoC1_11celltypes_centroid_allgenes_std.txt"),header=T,row.names=1)
markers=read.table("Merged4-1235-NoC5-NoC1_11celltypes_markersall_mindiff0.2_logfc2fold_3.20.2019.txt",header=T,row.names=1,stringsAsFactors=F)

### Visualize markers in heatmap across all cell types
genes=markers$gene
data.use=centroid.std

levels=colnames(centroid.std)

colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
col.lab=rep("",length(levels))
col.lab=gsub(".*_","",levels)

ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cell Type")
colnames(clab)=c("Cell Type","")

col.use=redblue100

data.use=centroid.std[markers$gene,]
row.lab=rownames(data.use)
jpeg(file=paste0(dgename,"centroid_std_markersall.jpeg"),res=300,height=2600,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
jpeg(file=paste0(dgename,"centroid_std_markersall2.jpeg"),res=300,height=1800,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()





