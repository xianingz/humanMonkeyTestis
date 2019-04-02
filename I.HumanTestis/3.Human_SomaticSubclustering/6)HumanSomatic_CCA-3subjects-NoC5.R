### Somatic SubClustering for Subject 1,2,3 only
# 1.8.2019 by Qianyi
### Corrected for batch effect using Multi-CCA
### notes: concerned that there are too many myoid cells in Human5, with progenitor cells spread evenly in the myoid cluster.
### to test whether the clustering was biased by Human subject 5

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
### Use Human Subject 1, 2, 3, and 5 
dataset=dataset[c(1:9,11)]
subject=unique(gsub("-.*","",dataset))
length(dataset) # 10
length(subject) # 4


dgefile=dgename="CCA-3subjects123-Somatic/"

### load object for CCA of 4 subjects after removing Cluster 5
load(file="Human4-1235-NoC5Somatic_CCA.Robj")
pbmc
load(file=paste0("Human4-1235-NoC5Somatic_CCAcelltypeNoC1.Robj"))
# labeled with 6 somatic cell types


# setup Seurat objects 
### since both count matrices have already filtered cells, we do no additional filtering here
cells.use=names(pbmc@ident)[which(pbmc@meta.data$protocol %in% subject)]
protocol=pbmc@meta.data$protocol[which(pbmc@meta.data$protocol %in% subject)]
dgedata=pbmc@raw.data[,cells.use]
dim(dgedata)  # [1] 35063  1487
nCellperGene <- rowSums(dgedata>0)
length(which(nCellperGene==0))     #11497
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgedata2=dgedata[genes.use,]
dim(dgedata2) # [1] 28991  1487
datalist=list()
for(i in 1:length(subject)){
### Extract Somatic cells for each human subject
cells.use=colnames(dgedata2)[which(protocol == subject[i])]
dgedata=dgedata2[,cells.use]
# nCellperGene <- rowSums(dgedata>0)
# genes.use=names(nCellperGene[which(nCellperGene!=0)])
# dgedata2=dgedata[genes.use,]
# skip this because I do not want to select hvg that are present in every dataset
  dge <- CreateSeuratObject(raw.data = dgedata)
  dge <- NormalizeData(object = dge)
  dge <- ScaleData(object = dge)
  dge <- FindVariableGenes(object = dge, do.plot = FALSE)
  dge@meta.data[,"protocol"] <- subject[i]
  datalist[[i]]=dge
}
datalist
[[1]]
An object of class seurat in project SeuratProject 
 28991 genes across 217 samples.

[[2]]
An object of class seurat in project SeuratProject 
 28991 genes across 622 samples.

[[3]]
An object of class seurat in project SeuratProject 
 28991 genes across 648 samples.


# we will take the union of the top 2k variable genes in each dataset for
# alignment, note that we use 1k genes in the manuscript examples, you can
# try different number of top HVG here with negligible changes to the overall results
hvglist=list()
for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(dge)
  hvglist[[i]] <- rownames(x = head(x = dge@hvg.info, n = 1000))
  print(length(hvglist[[i]])) # 1000
}
hvg.union=unique(unlist(hvglist))
length(hvg.union) # 2415

# extract hvg that are present in every dataset
# it may be unfair to require hvg to be present in every dataset
# because Human5 may have hvg specific to that subject
# so skip this
for(i in 1:length(subject)){
  dge=datalist[[i]]
  hvg.union=hvg.union[which(hvg.union %in% rownames(dge@data))]
  print(length(hvg.union))
}

# check the number of selected hvg that fall within top 2k hvg for each dataset
hvg.union=unique(unlist(hvglist))
length(hvg.union) # 2415

for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(length(which(hvg.union %in% rownames(x = head(x = dge@hvg.info, n = 2000)))))
}
[1] 1107
[1] 1136
[1] 1135

# check the number of hvg overlapped between datasets
cc=matrix(0,length(subject),length(subject))
for(i in 1:length(subject)){
    for(j in 1:length(subject)){
        cc[i,j]=length(intersect(rownames(x = head(x = datalist[[i]]@hvg.info, n = 1000)),rownames(x = head(x = datalist[[j]]@hvg.info, n = 1000))))
    }
}
cc
     [,1] [,2] [,3]
[1,] 1000  252  230
[2,]  252 1000  278
[3,]  230  278 1000

### run multi CCA
pbmc <- RunMultiCCA(object.list = datalist, genes.use = hvg.union, num.ccs=20)
save(pbmc,file=paste0("Human3-123-NoC5Somatic_CCA.Robj"))


dgefile=dgename="CCA-3subjects123-Somatic/"

p1 <- DimPlot(object = pbmc, reduction.use = "cca",cols.use=myBrewerPalette, group.by = "protocol", pt.size = 1.5, 
    do.return = TRUE)
p2 <- VlnPlot(object = pbmc, features.plot = "CC1",cols.use=myBrewerPalette, group.by = "protocol", do.return = TRUE)
pdf(paste(dgefile,"dge_CCA_orig.pdf",sep=""),height=4.5,width=11)
plot_grid(p1, p2)
dev.off()


pdf(paste(dgefile,"dge_CCs_top.pdf",sep=""),height=5)
#DimElbowPlot(object = pbmc,reduction.type="cca",dims.plot=20)
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 1:6, 
    do.balanced = TRUE)
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 7:12, 
    do.balanced = TRUE)
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 13:18, 
    do.balanced = TRUE)
dev.off()

numCCs=10

###### Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
pbmc <- AlignSubspace(object = pbmc, reduction.type = "cca", grouping.var = "protocol", 
    dims.align = 1:numCCs)

### Visualize the aligned CCA and perform integrated analysis

p1 <- VlnPlot(object = pbmc, features.plot = "CC1", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p2 <- VlnPlot(object = pbmc, features.plot = "CC2", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p3 <- VlnPlot(object = pbmc, features.plot = "CC3", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p4 <- VlnPlot(object = pbmc, features.plot = "CC4", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p11 <- VlnPlot(object = pbmc, features.plot = "ACC1", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p22 <- VlnPlot(object = pbmc, features.plot = "ACC2", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p33 <- VlnPlot(object = pbmc, features.plot = "ACC3", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p44 <- VlnPlot(object = pbmc, features.plot = "ACC4", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
pdf(paste(dgefile,"dge_ACCs_top.pdf",sep=""),height=6,width=15)
plot_grid(p1, p2,p3,p4,p11,p22,p33,p44,ncol=4)
dev.off()

###### Now we can run a single integrated analysis on all cells!

pbmc <- RunTSNE(object = pbmc, reduction.use = "cca.aligned", dims.use = 1:numCCs, 
    do.fast = TRUE)
pbmc <- RunUMAP(pbmc, reduction.use = "cca.aligned", dims.use = 1:numCCs)

###### Visualize the original batch identity in ACC and tSNE spaces
plotlist=list()
plotlist[[1]]=DimPlot(pbmc,group.by = "protocol",reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE)
plotlist[[2]]=DimPlot(pbmc,group.by = "protocol",reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE)
plotlist[[3]]=DimPlot(pbmc,group.by = "protocol",reduction.use = "umap",cols.use=myBrewerPalette,do.return=TRUE)
plotlist[[4]]=TSNEPlot(pbmc,group.by = "protocol",colors.use=myBrewerPalette,do.return=TRUE)
pdf(paste(dgefile,"ACCs_tSNE_orig.pdf",sep=""),height=8,width=10.5)
multiplot(plotlist,cols = 2)
dev.off()
### plot ACCs and tSNE for each batch
pbmc=SetAllIdent(pbmc,id="orig.ident")
myBrewerPalette <- brewer.pal(12,"Paired")[c(1:6,9:10)] # used this for 8 batches
cols=myBrewerPalette
# note: fix xlim and ylim in DimPlot
# ACC2 +coord_cartesian(ylim=c(-5,2),xlim=c(-1,4.5))
# tSNE +coord_cartesian(ylim=c(-40,40),xlim=c(-38,35))
sets=unique(pbmc@meta.data$protocol)
plot2set=plot3set=plot4set=plottset=NULL
for(i in 1:length(sets)){
set=sets[i]
plot2set[[i]]=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,cols.use=cols[i],do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(pbmc@meta.data)[which(gsub("_.*","",rownames(pbmc@meta.data))==set)])
}
for(i in 1:length(sets)){
set=sets[i]
plot3set[[i]]=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,cols.use=cols[i],do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(pbmc@meta.data)[which(gsub("_.*","",rownames(pbmc@meta.data))==set)])
}
for(i in 1:length(sets)){
set=sets[i]
plot4set[[i]]=DimPlot(pbmc,reduction.use = "cca.aligned",1,4,cols.use=cols[i],do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(pbmc@meta.data)[which(gsub("_.*","",rownames(pbmc@meta.data))==set)])
}
for(i in 1:length(sets)){
set=sets[i]
plottset[[i]]=TSNEPlot(pbmc,colors.use=cols[i],do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(pbmc@meta.data)[which(gsub("_.*","",rownames(pbmc@meta.data))==set)])
}
pdf(paste(dgefile,"origSet-AAC_tSNE.pdf",sep=""),height=6,width=12)
par(mar=c(1,1,3,1))
multiplot(plot2set,cols = 4)
#multiplot(plot3set,cols = 4)
#multiplot(plot4set,cols = 4)
multiplot(plottset,cols = 4)
dev.off()
pdf(paste(dgefile,"origSet-AAC.pdf",sep=""),height=6,width=12)
par(mar=c(1,1,3,1))
multiplot(plot2set,cols = 4)
dev.off()
pdf(paste(dgefile,"origSet-tSNE.pdf",sep=""),height=6,width=12)
par(mar=c(1,1,3,1))
multiplot(plottset,cols = 4)
dev.off()


###### Louvain-jaccard clustering
pbmc <- FindClusters(object = pbmc, reduction.type = "cca.aligned", dims.use = 1:numCCs, 
    resolution=seq(0.1,3,by=0.1),save.SNN = TRUE)
save(pbmc,file=paste0("Human3-123-NoC5Somatic_CCA.Robj"))

print(c( length(unique(pbmc@meta.data$res.0.1)),length(unique(pbmc@meta.data$res.0.2)),length(unique(pbmc@meta.data$res.0.3)),length(unique(pbmc@meta.data$res.0.4)),length(unique(pbmc@meta.data$res.0.5)),length(unique(pbmc@meta.data$res.0.6)),length(unique(pbmc@meta.data$res.0.7)),length(unique(pbmc@meta.data$res.0.8)),length(unique(pbmc@meta.data$res.0.9)),length(unique(pbmc@meta.data$res.1)) ))
print(c( length(unique(pbmc@meta.data$res.1.1)),length(unique(pbmc@meta.data$res.1.2)),length(unique(pbmc@meta.data$res.1.3)),length(unique(pbmc@meta.data$res.1.4)),length(unique(pbmc@meta.data$res.1.5)),length(unique(pbmc@meta.data$res.1.6)),length(unique(pbmc@meta.data$res.1.7)),length(unique(pbmc@meta.data$res.1.8)),length(unique(pbmc@meta.data$res.1.9)),length(unique(pbmc@meta.data$res.2)) ))



ACCPlot <- function(object, ...) {
  return(DimPlot(object = object, reduction.use = "cca.aligned", label.size = 4, ...))
}

res=paste0("res.",c(paste0("0.",1:9),1)) # res.0.8 for 10 clusters
plotlist=list()
pdf(paste(dgefile,"clusters.pdf",sep=""),height=7.5,width=9)
for(resi in 1:length(res)){
pbmc=SetAllIdent(pbmc,id=res[resi])
plotlist[[1]]=DimPlot(pbmc,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(pbmc,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=TSNEPlot(pbmc,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)
}
dev.off()


###### order clusters for each dataset
res=c("res.0.2","res.0.5")
for(resi in 1:length(res)){

pbmc=SetAllIdent(pbmc,id=res[resi])
pbmc <- BuildClusterTree(pbmc, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs)
pbmc@meta.data[,res[resi]]=pbmc@ident
levels=levels(pbmc@ident)
ident=factor(pbmc@ident,levels=levels)

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
tmppbmc=data.frame(t(as.matrix(pbmc@data[,cells.use])))
# make sure same order for cells.ident and pbmc before combining
which(names(cells.ident)!=colnames(pbmc@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(pbmc@data[,cells.use])))
genecountsall=matrix(,dim(pbmc@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(pbmc@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(subject)
bb=1
pdf(file=paste0(dgename,"Centroid_norm_Seriation_Aligned_",res[resi],".pdf"))
tmp=genecountsall[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) # 91
da
### methods of seriation
# methods <- c("ARSA","HC_single", "HC_complete", "OLO", "GW", "R2E", "VAT","TSP", "Spectral", "SPIN", "MDS", "Identity", "Random")
# note: the default method fo dissplot is "ARSA"; the default method for hmap is "OLO"
# order by seriation
do=seriate(da,method="ARSA")
# do=seriate(da,method="OLO")
# get order of seriation
print(get_order(seriate(da,method="ARSA")))
print(get_order(seriate(da,method="Spectral")))
print(get_order(seriate(da,method="OLO")))
 # plot original matrix
 dissplot(da, method = NA,options = list(main = paste("Dissimilarity")))
 pimage(da) # original matrix
 # plot with seriation
 dissplot(da, options = list(main = paste("Dissimilarity with seriation"))) # default method="ARSA"
 dissplot(da, method="ARSA", options = list(main = paste("Dissimilarity with seriation ARSA")))
 dissplot(da, method="Spectral",options = list(main = paste("Dissimilarity with seriation Spectral"))) 
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO")))
 hmap(da) # default method="OLO"

# dataset[bb]
dev.off()

levelss=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss
levels=levelss

### Calculate correlation for each normalized centroid
cc=cor(as.matrix(genecountsall),method="spearman")
dim(cc)
min(cc)

### Reordered cluster centroid 
data.use=cc[levels,levels]
colsep.use=cumsum(table(gsub("_.*","",levels)))
row.lab=col.lab=gsub(".*_","",levels)

ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

col.use=redblue100

pdf(file=paste(dgename,"Centroid_RankedCorrelation_Aligned_",res[resi],".pdf",sep=""),height=5.5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()

### Reordered clusters for all cells
cells.use=colnames(pbmc@data)
# random shuffling cells within ordered clusters
ident=factor(pbmc@ident,levels=levels-1)

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

### save ordered cluster ID in pbmc object
which(unique(cells.ident)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"order")

pbmc=AddMetaData(pbmc,cells.ident.ordered,ordered)
pbmc@meta.data[,ordered]=factor(pbmc@meta.data[,ordered])
pbmc=SetAllIdent(pbmc,ordered)
}
# save the pbmc file
save(pbmc,file=paste0("Human3-123-NoC5Somatic_CCA.Robj"))
pbmc3=pbmc
# go above to re-plot PCA and tSNE





###### re-plot PCA and tSNE for ordered clusters
pbmc=SetAllIdent(pbmc,"orderedclusters")
myBrewerPalette <- brewer.pal(12,"Paired")[c(4:1,5:6,9:12)] 
# used this color scheme for 9 clusters
pdf(paste(dgefile,"Aligned_clusters_ordered.pdf",sep=""),height=4.5,width=5)
DimPlot(pbmc,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
TSNEPlot(object = pbmc, do.label = TRUE, colors.use=myBrewerPalette,label.size=4,pt.size=1.2)
dev.off()
#pdf(paste(dgefile,"Aligned_orig.pdf",sep=""),height=4,width=5)
#TSNEPlot(object = pbmc, group.by = "protocol", do.return = TRUE,pt.size=1.2)
#dev.off()
table(pbmc@meta.data$protocol,pbmc@meta.data$orderedclusters)
table(pbmc@ident,dge4@ident)
table(pbmc@ident[names(dgeold@ident)],dgeold@ident)

load(file="Human4-1235-NoC5Somatic_CCAcelltypeNoC1.Robj")
pbmc # 35063 genes across 3722 samples.
pbmc4=pbmc
table(pbmc@ident[names(pbmc3@ident)],pbmc3@ident)

load(file="Human4-1235-NoC5Somatic_CCA.Robj")
pbmc # 35063 genes across 3722 samples.
pbmc4=pbmc
table(pbmc4@ident[names(pbmc3@ident)],pbmc3@ident)

sets=unique(pbmc@meta.data$protocol)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette
for(i in 1:length(sets)){
set=sets[i]
plot2set[[i]]=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,cols.use=myBrewerPalette,do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(pbmc@meta.data)[which(gsub("_.*","",rownames(pbmc@meta.data))==set)])
}
for(i in 1:length(sets)){
set=sets[i]
plottset[[i]]=TSNEPlot(pbmc,colors.use=myBrewerPalette,do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(pbmc@meta.data)[which(gsub("_.*","",rownames(pbmc@meta.data))==set)])
}
pdf(paste(dgefile,"origSet_clusters.pdf",sep=""),height=6,width=12)
par(mar=c(1,1,3,1))
multiplot(plot2set,cols = 4)
multiplot(plottset,cols = 4)
dev.off()


### visualize known markers
knownmarkers=c("VIM",
  "CD163","S100A4","TYROBP","LYZ","RGS1",
  "VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN",
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","PDGFRA","PDGFRB","PDGFB",
  "DLK1","HSD17B3","STAR","CYP17A1","NR5A1","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "CLU","SOX9","AMH"  )
length(knownmarkers) # 34
knownmarkers[which(!(knownmarkers %in% rownames(pbmc@data)))] 
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(pbmc@data))]
length(knownmarkers) # 34
pdf(paste0(dgefile,"knownmarkers_Feature.pdf"),height=15,width=21)
FeaturePlot(object = pbmc, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf(paste0(dgefile,"knownmarkers_Violin.pdf"),height=10,width=22)
VlnPlot(pbmc,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()

# 1.26.2019 Visualize IntProg/Myoid/Leydig markers 
### all known IntProg/Myoid/Leydig markers
knownmarkers=c(
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","TCF21","PDGFRA","PDGFRB","PDGFB",
  "DLK1","HSD17B3","STAR","CYP17A1","GATA4","NR5A1","EGR2","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "THY1","ENG","VCAM1","ALCAM","CD44","NT5E","PDGFA","CASP3","CLU" )
length(knownmarkers) # 31
knownmarkers[which(!(knownmarkers %in% rownames(pbmc@data)))] #[1] "PDGFRA" "CASP3"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(pbmc@data))]
length(knownmarkers) # 31

pdf("knownmarkers_IntProg_Feature.pdf",height=15,width=21)
FeaturePlot(object = pbmc, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf("knownmarkers_IntProg_Violin.pdf",height=10,width=22)
VlnPlot(pbmc,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()


###### assign cell types based on markers for each ordered clusters in each dataset
pbmc=SetAllIdent(pbmc,id="res.0.2")
markers=FindAllMarkers(pbmc,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,"Human3-123-NoC5Somatic_res.0.2order_markersall_mindiff0.2_logfc2fold_1.30.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(pbmc@ident)
table(markers$cluster)
  1   2   3   4   5   6 
 19 129  74 368 790 107 
194 203  92 119  62  66


pbmc=SetAllIdent(pbmc,id="res.0.5order")
markers=FindAllMarkers(pbmc,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,"Human3-123-NoC5Somatic_res.0.5order_markersall_mindiff0.2_logfc2fold_1.30.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(pbmc@ident)
table(markers$cluster)
  1   2   3   4   5   6   7   8 
 20 133  74 369 367 384 107  33 
164 203  92 119  31  11  66 337 

markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf(paste0(dgefile,"markerstop_res.0.5order.pdf"),height=12,width=15)
FeaturePlot(object = pbmc, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=5)
dev.off()

### PCA and tSNE plot
pbmc=SetAllIdent(pbmc,id="orig.ident")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-123-NoC5Somatic_Rep.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
pbmc=SetAllIdent(pbmc,id="protocol")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-123-NoC5Somatic_Subject.pdf"),width=10,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
pbmc=SetAllIdent(pbmc,id="res.0.2order")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-123-NoC5Somatic_orderedclusters.pdf"),width=9.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
pbmc=SetAllIdent(pbmc,id="res.0.5order")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-123-NoC5Somatic_res.0.5orderedclusters.pdf"),width=9.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
