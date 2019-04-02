### Somatic Subclustering for 5 monkey subjects corrected for batch effect by CCA after removing doublets cluster 1 and 7
# 2.22.2019 by Qianyi
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering 
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1 cells after removing C1 doublets
### Removed Cluster 7/7 of SCyte-Somatic doublets from somatic subclustering by CCA
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing C1 and C7 doublets

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



### load object for CCA of 5 subjects after NoC1 before removing C7
# CCA and somatic subclustering after removing Cluster 1/8 of spermatid-myoid doublets for somatic subclustering 
load(file="Monkey5Somatic-NoC1_CCA.Robj")
testis1=pbmc
testis1 # 18370 genes across 2355 samples.

### directly-merged 5 subjects labeled NoC1
# did somatic subclustering before removing cluster 1/8 doublets
# labeled 6 cell types after removing Cluster 1/8 of spermatid-myoid doublets from somatic subclustering 
load(file="MonkeyMerged5SomaticcelltypeNoC1.Robj")
dge5celltype=dge
dge5=dge


library(Seurat)
# cowplot enables side-by-side ggplots
library(cowplot)


# setup Seurat objects 
### remove 7/7 from previous CCA and somatic subclustering
table(testis1@ident)
   1    2    3    4    5    6    7 
  37  175   42  180  521 1382   18 
testis=SubsetData(testis1,ident.remove=7)
table(testis@ident)
   1    2    3    4    5    6 
  37  175   42  180  521 1382 
nCellperGene <- rowSums(testis@data>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # 17469

datalist=list()
for(i in 1:length(subject)){
### Extract Somatic-NoC1-NoC7 cells for each monkey subject
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
[[1]]
An object of class seurat in project SeuratProject 
 17469 genes across 386 samples.

[[2]]
An object of class seurat in project SeuratProject 
 17469 genes across 160 samples.

[[3]]
An object of class seurat in project SeuratProject 
 17469 genes across 1436 samples.

[[4]]
An object of class seurat in project SeuratProject 
 17469 genes across 208 samples.

[[5]]
An object of class seurat in project SeuratProject 
 17469 genes across 147 samples.


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
length(hvg.union) # 3531

# extract hvg that are present in every dataset
# it may be unfair to require hvg to be present in every dataset
# because Monkey5 may have hvg specific to that subject
# so skip this
for(i in 1:length(subject)){
  dge=datalist[[i]]
  hvg.union=hvg.union[which(hvg.union %in% rownames(dge@data))]
  print(length(hvg.union))
}

# check the number of selected hvg that fall within top 2k hvg for each dataset
hvg.union=unique(unlist(hvglist))
length(hvg.union) # 3531

for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(length(which(hvg.union %in% rownames(x = head(x = dge@hvg.info, n = 2000)))))
}
[1] 1338
[1] 1255
[1] 1313
[1] 1289
[1] 1267


# check the number of hvg overlapped between datasets
cc=matrix(0,length(subject),length(subject))
for(i in 1:length(subject)){
    for(j in 1:length(subject)){
        cc[i,j]=length(intersect(rownames(x = head(x = datalist[[i]]@hvg.info, n = 1000)),rownames(x = head(x = datalist[[j]]@hvg.info, n = 1000))))
    }
}
cc
     [,1] [,2] [,3] [,4] [,5]
[1,] 1000  220  286  290  265
[2,]  220 1000  219  231  218
[3,]  286  219 1000  278  257
[4,]  290  231  278 1000  263
[5,]  265  218  257  263 1000



### run multi CCA
testis <- RunMultiCCA(object.list = datalist, genes.use = hvg.union, num.ccs=20)
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA.Robj"))


dgefile=dgename="CCA5Subjects-Somatic-NoC1-NoC7/"

p1 <- DimPlot(object = testis, reduction.use = "cca",cols.use=myBrewerPalette, group.by = "protocol", pt.size = 1.5, 
    do.return = TRUE)
p2 <- VlnPlot(object = testis, features.plot = "CC1",cols.use=myBrewerPalette, group.by = "protocol", do.return = TRUE)
pdf(paste(dgefile,"dge_CCA_orig.pdf",sep=""),height=4.5,width=11)
plot_grid(p1, p2)
dev.off()

pdf(paste(dgefile,"dge_CCs_top.pdf",sep=""),height=5)
MetageneBicorPlot(testis, grouping.var = "protocol", dims.eval = 1:20, 
    display.progress = FALSE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 1:6, 
    do.balanced = TRUE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 7:12, 
    do.balanced = TRUE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 13:18, 
    do.balanced = TRUE)
dev.off()

#Rescaling group 1
#Rescaling group 2
#Rescaling group 3
#Rescaling group 4
#Rescaling group 5
#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
#`geom_smooth()` using method = 'loess' and formula 'y ~ x'

numCCs=12
###### Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
testis <- AlignSubspace(object = testis, reduction.type = "cca", grouping.var = "protocol", 
    dims.align = 1:numCCs)

### Visualize the aligned CCA and perform integrated analysis

p1 <- VlnPlot(object = testis, features.plot = "CC1", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p2 <- VlnPlot(object = testis, features.plot = "CC2", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p3 <- VlnPlot(object = testis, features.plot = "CC3", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p4 <- VlnPlot(object = testis, features.plot = "CC4", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p11 <- VlnPlot(object = testis, features.plot = "ACC1", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p22 <- VlnPlot(object = testis, features.plot = "ACC2", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p33 <- VlnPlot(object = testis, features.plot = "ACC3", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p44 <- VlnPlot(object = testis, features.plot = "ACC4", group.by = "protocol",cols.use=myBrewerPalette, 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
pdf(paste(dgefile,"dge_ACCs_top.pdf",sep=""),height=6,width=15)
plot_grid(p1, p2,p3,p4,p11,p22,p33,p44,ncol=4)
dev.off()

###### Now we can run a single integrated analysis on all cells!

testis <- RunTSNE(object = testis, reduction.use = "cca.aligned", dims.use = 1:numCCs, 
    do.fast = TRUE)
testis <- RunUMAP(testis, reduction.use = "cca.aligned", dims.use = 1:numCCs)
testis <- RunPCA(testis, pc.genes = testis@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
testis <- ProjectPCA(object = testis, do.print = FALSE)

###### Louvain-jaccard clustering
testis <- FindClusters(object = testis, reduction.type = "cca.aligned", dims.use = 1:numCCs, 
    resolution=seq(0.1,3,by=0.1),save.SNN = TRUE)
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA.Robj"))

print(c( length(unique(testis@meta.data$res.0.1)),length(unique(testis@meta.data$res.0.2)),length(unique(testis@meta.data$res.0.3)),length(unique(testis@meta.data$res.0.4)),length(unique(testis@meta.data$res.0.5)),length(unique(testis@meta.data$res.0.6)),length(unique(testis@meta.data$res.0.7)),length(unique(testis@meta.data$res.0.8)),length(unique(testis@meta.data$res.0.9)),length(unique(testis@meta.data$res.1)) ))
print(c( length(unique(testis@meta.data$res.1.1)),length(unique(testis@meta.data$res.1.2)),length(unique(testis@meta.data$res.1.3)),length(unique(testis@meta.data$res.1.4)),length(unique(testis@meta.data$res.1.5)),length(unique(testis@meta.data$res.1.6)),length(unique(testis@meta.data$res.1.7)),length(unique(testis@meta.data$res.1.8)),length(unique(testis@meta.data$res.1.9)),length(unique(testis@meta.data$res.2)) ))


res=paste0("res.",c(paste0("0.",1:9),1,paste0("1.",1:9))) # res.0.8 for 10 clusters
plotlist=list()
pdf(paste(dgefile,"clusters.pdf",sep=""),height=7.5,width=9)
for(resi in 1:length(res)){
testis=SetAllIdent(testis,id=res[resi])
plotlist[[1]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=DimPlot(testis,reduction.use = "umap",cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[4]]=TSNEPlot(testis,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)

}
dev.off()

table(dge5@ident[names(testis@ident)],testis@meta.data$res.0.1)
table(dge5@ident[names(testis@ident)],testis@meta.data$res.0.2)
table(dge5@ident[names(testis@ident)],testis@meta.data$res.0.3)
table(dge5@ident[names(testis@ident)],testis@meta.data$res.0.4)
                   0    1    2    3    4    5
  Macrophage       0    0    0    0    0   37
  Endothelial     18    0  170    2    0    0
  Leydig           5    1    0    0   42    0
  Pericyte        34    0    0  149    0    0
  Myoid          129  473    0    1    1    0
  IntProgenitor 1246   28    0    1    0    0

table(dge5@ident[names(testis@ident)],testis@meta.data$res.0.5)
                   0    1    2    3    4    5
  Macrophage       0    0    0    0    0   37
  Endothelial     18    0    2  170    0    0
  Leydig           4    1    0    0   43    0
  Pericyte        13    0  170    0    0    0
  Myoid          122  480    1    0    1    0
  IntProgenitor 1239   35    1    0    0    0


table(testis1@ident[names(testis@ident)],testis@meta.data$res.0.4)
       0    1    2    3    4    5
  1    0    0    0    0    0   37
  2   14    0  161    0    0    0
  3    0    1    0    0   41    0
  4   29    0    0  151    0    0
  5   45  475    0    1    0    0
  6 1344   26    9    1    2    0
  7    0    0    0    0    0    0

table(testis1@ident[names(testis@ident)],testis@meta.data$res.0.5)
       0    1    2    3    4    5
  1    0    0    0    0    0   37
  2   14    0    0  161    0    0
  3    0    1    0    0   41    0
  4   12    0  168    0    0    0
  5   35  485    1    0    0    0
  6 1335   30    5    9    3    0
  7    0    0    0    0    0    0
table(testis@meta.data$res.0.3,testis@meta.data$res.0.5)
table(testis@meta.data$res.0.4,testis@meta.data$res.0.5)


TSNEPlot(testis,group.by="res.0.5")



###### order clusters for each dataset
res="res.0.5"
resi=1
testis=SetAllIdent(testis,id=res[resi])
testis <- BuildClusterTree(testis, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs)
testis@meta.data[,res[resi]]=testis@ident
levels=levels(testis@ident)
ident=factor(testis@ident,levels=levels)

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
tmptestis=data.frame(t(as.matrix(testis@data[,cells.use])))
# make sure same order for cells.ident and testis before combining
which(names(cells.ident)!=colnames(testis@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(testis@data[,cells.use])))
genecountsall=matrix(,dim(testis@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(testis@data)
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
pdf(file=paste0(dgename,"Centroid_norm_Seriation_",res[resi],".pdf"))
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

pdf(file=paste(dgename,"Centroid_RankedCorrelation_",res[resi],".pdf",sep=""),height=5.5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()

### Reordered clusters for all cells
cells.use=colnames(testis@data)
# random shuffling cells within ordered clusters
ident=factor(testis@ident,levels=levels-1)

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

### save ordered cluster ID in testis object
which(unique(cells.ident)!=get_order(do)) # integer(0)

ordered="orderedclusters"

testis=AddMetaData(testis,cells.ident.ordered,ordered)
testis@meta.data[,ordered]=factor(testis@meta.data[,ordered])
testis=SetAllIdent(testis,ordered)
# save the testis file
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA.Robj"))

# go above to re-plot PCA and tSNE



### compare with previous clusters
table(dge5celltype@ident[names(testis@ident)],testis@ident)
                   1    2    3    4    5    6
  Macrophage      37    0    0    0    0    0
  Endothelial      0  170    2    0   18    0
  Leydig           0    0    0    1    4   43
  Pericyte         0    0  170    0   13    0
  Myoid            0    0    1  480  122    1
  IntProgenitor    0    0    1   35 1239    0

table(testis1@ident[names(testis@ident)],testis@ident)
       1    2    3    4    5    6
  1   37    0    0    0    0    0
  2    0  161    0    0   14    0
  3    0    0    0    1    0   41
  4    0    0  168    0   12    0
  5    0    0    1  485   35    0
  6    0    9    5   30 1335    3
  7    0    0    0    0    0    0


###### re-plot PCA and tSNE for ordered clusters
testis=SetAllIdent(testis,"orderedclusters")
myBrewerPalette <- brewer.pal(12,"Paired")[c(4:1,5:6,9:12)] 
# used this color scheme for 9 clusters
pdf(paste(dgefile,"Aligned_clusters_ordered.pdf",sep=""),height=4.5,width=5)
DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
TSNEPlot(object = testis, do.label = TRUE, colors.use=myBrewerPalette,label.size=4,pt.size=1.2)
dev.off()
#pdf(paste(dgefile,"Aligned_orig.pdf",sep=""),height=4,width=5)
#TSNEPlot(object = testis, group.by = "protocol", do.return = TRUE,pt.size=1.2)
#dev.off()
table(testis@meta.data$protocol,testis@meta.data$orderedclusters)
table(testis@ident,dge4@ident)
table(testis@ident[names(dgeold@ident)],dgeold@ident)

sets=unique(testis@meta.data$protocol)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette
for(i in 1:length(sets)){
set=sets[i]
plot2set[[i]]=DimPlot(testis,reduction.use = "cca.aligned",1,2,cols.use=myBrewerPalette,do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(testis@meta.data)[which(gsub("_.*","",rownames(testis@meta.data))==set)])
}
for(i in 1:length(sets)){
set=sets[i]
plottset[[i]]=TSNEPlot(testis,colors.use=myBrewerPalette,do.return = TRUE,do.label=T,no.legend=TRUE,pt.size = 1,cells.use=rownames(testis@meta.data)[which(gsub("_.*","",rownames(testis@meta.data))==set)])
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
  "DCN","PDGFRA","PDGFRB","PDGFB","ALDH1",
  "DLK1","HSD17B3","STAR","CYP17A1","NR5A1","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "CLU","SOX9","AMH"  )
length(knownmarkers) # 35
knownmarkers[which(!(knownmarkers %in% rownames(testis@data)))] #[1] "ALDH1"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(testis@data))]
length(knownmarkers) # 32
pdf(paste0(dgefile,"knownmarkers_Feature.pdf"),height=15,width=21)
FeaturePlot(object = testis, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf(paste0(dgefile,"knownmarkers_Violin.pdf"),height=10,width=22)
VlnPlot(testis,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()

knownmarkers=c("NR2F2","TCF21","GATA4","NR5A1","EGR2","PDGFRA","PDGFA","PDGFRB","PDGFB",
    "THY1","ENG","VCAM1","ALCAM","CD44","NT5E","CASP3" )
# NGF1B: EGR2

length(knownmarkers) # 16
knownmarkers[which(!(knownmarkers %in% rownames(testis@data)))] 
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(testis@data))]
length(knownmarkers) # 14

pdf(paste0(dgefile,"knownmarkers_MSC_Feature.pdf"),height=12,width=12)
FeaturePlot(object = testis, features.plot = knownmarkers,col=c("grey80","red"),nCol=4)
dev.off()
pdf(paste0(dgefile,"knownmarkers_MSC_Violin.pdf"),height=8,width=13)
VlnPlot(testis,knownmarkers,cols.use=myBrewerPalette,nCol=4,point.size.use=-1)
dev.off()

### all known IntProg/Myoid/Leydig markers
knownmarkers=c(
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","TCF21","PDGFRA","PDGFRB","PDGFB",
  "DLK1","HSD17B3","STAR","CYP17A1","GATA4","NR5A1","EGR2","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "THY1","ENG","VCAM1","ALCAM","CD44","NT5E","PDGFA","CASP3","CLU" )
length(knownmarkers) # 31
knownmarkers[which(!(knownmarkers %in% rownames(testis@data)))] #[1] "PDGFRA" "CASP3"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(testis@data))]
length(knownmarkers) # 31

pdf("knownmarkers_IntProg_Feature.pdf",height=15,width=21)
FeaturePlot(object = testis, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf("knownmarkers_IntProg_Violin.pdf",height=10,width=22)
VlnPlot(testis,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()

### TCF21 paralogs from Adrienne
knownmarkers=c(
  "TCF24","TCF23",
  "TCF15","MSC","TWIST1","TWIST2",
  "TCF4","TCF19","TCF7","TCF12",
  "TCF7L1","TCF25","TCF19","TCF7L2","TCF20",
  "TCF3" ,
  "HAND2","BHLHA9","PTF1A","HAND1","TCF7L1-IT1",
  "SCX",
  "FIGLA","FERD3L")
length(knownmarkers) # 24
knownmarkers[which(!(knownmarkers %in% rownames(testis@data)))] 
#[1] "FIGLA"  "FERD3L" "SCX"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(testis@data))]
length(knownmarkers) # 21

pdf("knownmarkers_TCF21paralogs_Feature.pdf",height=9,width=21)
FeaturePlot(object = testis, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf("knownmarkers_TCF21paralogs_Violin.pdf",height=6,width=22)
VlnPlot(testis,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()




###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(testis,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,"res.0.5order_markersall_mindiff0.2_logfc2fold_2.22.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(testis@ident)
table(markers$cluster)
   1    2    3    4    5    6 
  37  170  174  516 1396   44 
  110 128 114  25  29  90 

markers=FindAllMarkers(testis,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0("monkey/NoC1NoC7CCA5Somatic_8celltypes_markersall_mindiff0.2_logfc2fold_3.13.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(testis@ident)
table(markers$cluster)
Tcell  Macrophage Endothelial    Pericyte Subcluster2       Myoid 
         20          17         170         131          43         516 
        123         173         128         150          57          25 
  ImmLeydig  DiffLeydig 
       1396          44 
         29          90 

markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf(paste0(dgefile,"markerstop.pdf"),height=9,width=12)
FeaturePlot(object = testis, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=4)
dev.off()

### PCA and tSNE plot
testis=SetAllIdent(testis,id="orig.ident")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1-NoC7_Rep.pdf"),width=10.5,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
testis=SetAllIdent(testis,id="protocol")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1-NoC7_Subject.pdf"),width=10,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
testis=SetAllIdent(testis,id="orderedclusters")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1-NoC7_orderedclusters.pdf"),width=9.5,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

###### Visualize individual batch and subject
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlims=ylims=list()
pos=c("bottomleft","topleft","topleft","bottomleft")
# testis4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-1.4,3.8),c(-1.4,3.8),c(-8,11),c(-28,35))
ylims[[1]]=list(c(-5.5,1.9),c(-3.5,4.5),c(-4.9,9.5),c(-45,35))

dim=list(testis@dr$cca.aligned@cell.embeddings[,1:2],testis@dr$cca.aligned@cell.embeddings[,c(1,3)],testis@dr$umap@cell.embeddings[,c(1,2)],testis@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

### plot PCs and tSNE for each batch using the other batches as background
testis=SetAllIdent(testis,id="orig.ident")
sets=levels(testis@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette[c(8:7,10,6:5,4:1)]

pdf(paste0(dgefile,"OrigSet_PCtSNE.pdf"),height=4.6,width=11.5)
par(mfrow=c(2,5),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:length(dim)){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(testis@ident)
names(ident)=names(testis@ident)
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
testis=SetAllIdent(testis,id="protocol")
sets=levels(testis@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

pdf(paste0(dgefile,"OrigSubject_PCtSNE.pdf"),height=4.6,width=6.9)
par(mfrow=c(2,3),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:length(dim)){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(testis@ident)
names(ident)=names(testis@ident)
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
testis=SetAllIdent(testis,id="orderedclusters")
table(testis@meta.data$protocol,testis@meta.data$orderedclusters)

  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(testis@meta.data$orig.ident,testis@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0(dgefile,"ncellspercluster_batch_top15CCs.txt"),quote=F,row.names=T,col.names=T,sep="\t")


testis=SetAllIdent(testis,id="orderedclusters")
sets=subject
which(rownames(data)!=names(testis@ident))


cols=myBrewerPalette1
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())


j=3
label="UMAP"
data=data_bg=as.data.frame(testis@dr$umap@cell.embeddings)
data$ident=testis@ident
data %>% dplyr::group_by(ident) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)) -> centers
label.size=6

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(testis@ident)[which(gsub("\\..*","",gsub("_.*","",names(testis@ident)))==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=UMAP1,y=UMAP2,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=myBrewerPalette,limits=levels(data$ident),drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}


j=4
label="tSNE"
data=data_bg=as.data.frame(testis@dr$tsne@cell.embeddings)
data$ident=testis@ident
data %>% dplyr::group_by(ident) %>% summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) -> centers
label.size=6

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(testis@ident)[which(gsub("\\..*","",gsub("_.*","",names(testis@ident)))==set)]
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
source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_testis/Rcode_multiplot.R")

jpeg(file=paste0(dgefile,"OrigSubject_7clusters_",label,"_bg.jpeg"),res=300,height=1200,width=1900)
multiplot(plotset,cols = 3)
dev.off()
pdf(file=paste0(dgefile,"OrigSubject_7clusters_",label,"_bg.pdf"),height=6,width=9)
multiplot(plotset,cols = 3)
dev.off()



# 2.27.2019 Qianyi
###### label cell types 
testis=SetAllIdent(testis,id="orderedclusters")
table(testis@ident)
   1    2    3    4    5    6
  37  170  174  516 1396   44

### label cell types
id=as.numeric(testis@ident)
names(id)=names(testis@ident)
celltype=c("Macrophage/Tcell","Endothelial","Pericyte","Myoid","ImmLeydig","DiffLeydig")
for(i in 1:length(celltype)){
id[which(id==i)]<-celltype[i]
}
id=factor(id,levels=celltype,order=T)
table(id)
Macrophage/Tcell      Endothelial         Pericyte            Myoid
              37              170              174              516
       ImmLeydig       DiffLeydig
            1396               44

testis=AddMetaData(testis,id,"CellType")
testis=SetAllIdent(testis,id="CellType")
testis@ident=id
table(testis@ident)
Macrophage/Tcell      Endothelial         Pericyte            Myoid
              37              170              174              516
       ImmLeydig       DiffLeydig
            1396               44

save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA.Robj"))
write.table(testis@ident,"Monkey5Somatic-NoC1-NoC7_CCA_6celltypes_label.txt",row.names=T,col.names=F,quote=F,sep="\t")

# try to match color scheme for corresponding mouse cell types
#levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown")
#myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7)] # 1

celltype=c("Macrophage/Tcell","Endothelial","Pericyte","Myoid","ImmLeydig","DiffLeydig")
testis=SetAllIdent(testis,id="CellType")
testis@ident=factor(testis@ident,levels=celltype)
table(testis@ident)
library(RColorBrewer)
myBrewerPalette1 <- c(brewer.pal(6,"Set1")) 
myBrewerPalette2 <- c(brewer.pal(12,"Paired"))[c(11,10,6,7,9,5)] 
# used this color scheme

pdf(paste0(dgefile,"PCA_tSNE_UMAP_6celltypes.pdf"),width=11,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

pdf(paste0(dgefile,"PCA_tSNE_UMAP_6celltypes1.pdf"),width=11,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette1)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette1)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot4=DimPlot(testis,reduction.use = "cca.aligned",1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot5=DimPlot(testis,reduction.use = "cca.aligned",1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(testis,1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot3=PCAPlot(testis,1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot4=PCAPlot(testis,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot5=PCAPlot(testis,1,5,do.return=TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=PCAPlot(testis,1,6,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot3=PCAPlot(testis,1,7,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot4=PCAPlot(testis,1,8,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot5=PCAPlot(testis,1,9,do.return=TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

pdf(paste0(dgefile,"PCA_tSNE_UMAP_6celltypes2.pdf"),width=11,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette2)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette2)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette2)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette2)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()


###### Rank correlation for all cells
nrho=cor(as.matrix(testis@data)[testis@var.genes,],method="spearman")
length(testis@var.genes) # 2494
testcor=nrho

###### Jaccard distance
testcor=as.matrix(testis@snn)


### order cells by clusters 
res="orderedclusters";j=1;resi=1;
   testis=SetAllIdent(testis,id=res)
   print(length(unique(testis@ident)))
   TSNEPlot(testis)

## order cell by cluster ID and randomly shuffle cells within each batch
levels=levels(testis@ident)
ident=factor(testis@ident,levels=levels)

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

            data.use2=testcor[cells.use,cells.use]
            #data.use2=minmax(data.use2,min=disp.min,max=disp.max)

            lab2=rep("",length(cells.use))
            lab2[round(cumsum(table(cells.ident)[levels(cells.ident)])-table(cells.ident)[levels(cells.ident)]/2)+1]=levels(cells.ident)

            row.lab2=gsub(".*_","",lab2)
            orig.ident=factor(gsub("_.*","",cells.ident),levels=unique(gsub("_.*","",cells.ident)))
            col.lab2=rep("",length(cells.use))
            col.lab2[round(cumsum(table(orig.ident)[levels(orig.ident)])-table(orig.ident)[levels(orig.ident)]/2)+1]=levels(orig.ident)
            colsep.use2=cumsum(table(cells.ident)[levels(cells.ident)])
            colsep.use2=cumsum(table(orig.ident)[levels(orig.ident)]) # draw a line between datasets
sidecol2=do.call(rbind,strsplit(as.character(cells.ident),"_"))
sidecol2=cbind(sidecol2,sidecol2)

library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])

rlab2=rbind(rep("white",length(sidecol2[,1])),myBrewerPalette[as.numeric(gsub(".*-","",sidecol2[,2]))])
clab2=cbind(rlab2[2,],rlab2[1,])
colnames(clab2)=c("Cluster","")
rownames(rlab2)=c("","Cluster")

midrange=median(testcor[which(testcor!=1)])
maxrange=max(testcor[which(testcor!=1)])
maxrange2=maxrange
midrange2=maxrange/2
col.use2=redblue100[c(rep(1:50,each=round(midrange2*100)),rep(50:100,each=round((maxrange2-midrange2)*100)),rep(100,50*round((1-maxrange2)*100)))]
length(col.use2)

jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7_6clusters_RankCor_HVG_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7_6clusters_RankCor_HVG.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()

col.use2=redblue100[c(rep(c(41:100),each=10),rep(100,100000))]
length(col.use2)
jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7_6clusters_snn_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7_6clusters_snn.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()


###### Rank correlation and Dissimilarity matrix for each normalized centroid using HVG
### for each cluster, calculate average normalized expression of each gene
tmptestis=data.frame(t(as.matrix(testis@data[,cells.use])))
# make sure same order for cells.ident and testis before combining
which(names(cells.ident)!=colnames(testis@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(testis@data[,cells.use])))
genecountsall=matrix(,dim(testis@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(testis@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

write.table(genecountsall,paste0(dgefile,"Merged5Somatic-NoC1-NoC7_CCA_6clusters_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)
write.table(genecountsall,paste0(dgefile,"Merged5Somatic-NoC1-NoC7_CCA_8celltypes_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)

cc=cor(as.matrix(genecountsall)[testis@var.genes,],method="spearman")
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

pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7_6clusters_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()

# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)):101]
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7_6clusters_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
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
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7_6celltypes_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()
# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)-1):101]
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7_6celltypes_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()


### Per-cell attributes
###### calculate Gini index
library(ineq)
#for(i in 1:(length(dataset)+1)){
#  testis=testislist[[i]]
GiniAll=apply( testis@raw.data,2,function(x) ineq(x,type=c("Gini")) )
GiniNon0=apply( testis@raw.data,2,function(x) ineq(x[which(x!=0)],type=c("Gini")) )

testis=AddMetaData(testis,GiniAll,"GiniAll")
testis=AddMetaData(testis,GiniNon0,"GiniNon0")


###### Correlation between Gini and nUMI
datainfo=testis@meta.data
datainfo_bg=datainfo[,! names(datainfo) %in% "CellType"]
myBrewerPalette=brewer.pal(7,"Set1")[c(1:5,7)]

### Gini Vs nUMI
pdf(file=paste(dgefile,"GiniVsnUMI.pdf",sep=""),height=4,width=5.5)
ggplot(datainfo,aes(x=nUMI,y=GiniAll,color=CellType))+
geom_point(size=1) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
dev.off()

# Gini for all genes of each cell
pdf(file=paste(dgefile,"GiniVsnUMI2.pdf",sep=""),height=3.5,width=6)
ggplot(datainfo,aes(x=nUMI,y=GiniAll,color=CellType))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~CellType,ncol=3) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
# Gini for non-0 genes of each cell
ggplot(datainfo,aes(x=nUMI,y=GiniNon0,color=CellType))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~CellType,ncol=3) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
# Gini for highly-variable genes of each cell
ggplot(datainfo,aes(x=nUMI,y=GiniHVG,color=CellType))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~CellType,ncol=3) +
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

### Visualize three outlier cells
rownames(datainfo)[which(datainfo$GiniAll>0.982)]
"Monkey2_GCAGTTCGTATC"   "Monkey4.1_ATGTCTAGCGTG" "Monkey4.2_CTAACCGAGATA"
which(datainfo$GiniAll>0.982)
[1]  545 2187 2190
which(datainfo$GiniNon0>0.52)
[1]  545 2187 2190

ident=rep(1,length(testis@ident))
names(ident)=names(testis@ident)
ident[which(testis@meta.data$GiniAll>0.982)] <- 2
ident=sort(ident)
par(mar=c(4,4,1,1),mgp=c(2,0.5,0))
plot(testis@dr$tsne@cell.embeddings[names(ident),1:2],cex=0.8,col=c("grey70","red")[ident],pch=16,xlab="tSNE_1",ylab="tSNE_2")

datainfo$nUMI[which(datainfo$GiniAll>0.982)]

pdf(file=paste0("PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(testis, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="GiniAll", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="GiniNon0", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()
