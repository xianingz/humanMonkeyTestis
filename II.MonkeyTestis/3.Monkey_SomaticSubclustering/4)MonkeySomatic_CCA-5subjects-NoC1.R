### Somatic Subclustering for 5 monkey subjects corrected for batch effect by CCA after removing doublets cluster 1
# 2.4.2019 by Qianyi
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering of directly-merged 5 subjects
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1 cells after removing doublets cluster 1

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



### directly-merged 5 subjects 
# did somatic subclustering before removing cluster 1/8 doublets
# labeled 6 cell types after removing Cluster 1/8 of spermatid-myoid doublets from somatic subclustering 
load(file="MonkeyMerged5SomaticcelltypeNoC1.Robj")
dge5celltype=dge
dge5=dge


library(Seurat)
# cowplot enables side-by-side ggplots
library(cowplot)


# setup Seurat objects 
### since both count matrices have already filtered cells, we do no additional filtering here
datalist=list()
for(i in 1:length(subject)){
### Extract Somatic-NoC1 cells for each monkey subject
cells.use=names(dge5@ident)[which(dge5@meta.data$indiv == subject[i])]
dgedata=dge5@raw.data[,cells.use]
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
 18370 genes across 390 samples.

[[2]]
An object of class seurat in project SeuratProject 
 18370 genes across 161 samples.

[[3]]
An object of class seurat in project SeuratProject 
 18370 genes across 1445 samples.

[[4]]
An object of class seurat in project SeuratProject 
 18370 genes across 211 samples.

[[5]]
An object of class seurat in project SeuratProject 
 18370 genes across 148 samples.

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
length(hvg.union) # 3520

# extract hvg that are present in every dataset
for(i in 1:length(subject)){
  dge=datalist[[i]]
  hvg.union=hvg.union[which(hvg.union %in% rownames(dge@data))]
  print(length(hvg.union))
}

# check the number of selected hvg that fall within top 2k hvg for each dataset
hvg.union=unique(unlist(hvglist))
length(hvg.union) # 3520

for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(length(which(hvg.union %in% rownames(x = head(x = dge@hvg.info, n = 2000)))))
}
[1] 1336
[1] 1246
[1] 1310
[1] 1294
[1] 1264


# check the number of hvg overlapped between datasets
cc=matrix(0,length(subject),length(subject))
for(i in 1:length(subject)){
    for(j in 1:length(subject)){
        cc[i,j]=length(intersect(rownames(x = head(x = datalist[[i]]@hvg.info, n = 1000)),rownames(x = head(x = datalist[[j]]@hvg.info, n = 1000))))
    }
}
cc
     [,1] [,2] [,3] [,4] [,5]
[1,] 1000  222  288  295  261
[2,]  222 1000  217  230  221
[3,]  288  217 1000  277  260
[4,]  295  230  277 1000  265
[5,]  261  221  260  265 1000


### run multi CCA
pbmc <- RunMultiCCA(object.list = datalist, genes.use = hvg.union, num.ccs=20)
save(pbmc,file=paste0("Monkey5Somatic-NoC1_CCA.Robj"))


dgefile=dgename="CCA5Subjects-Somatic-NoC1/"

p1 <- DimPlot(object = pbmc, reduction.use = "cca",cols.use=myBrewerPalette, group.by = "protocol", pt.size = 1.5, 
    do.return = TRUE)
p2 <- VlnPlot(object = pbmc, features.plot = "CC1",cols.use=myBrewerPalette, group.by = "protocol", do.return = TRUE)
pdf(paste(dgefile,"dge_CCA_orig.pdf",sep=""),height=4.5,width=11)
plot_grid(p1, p2)
dev.off()
pdf(paste(dgefile,"dge_CCs_top.pdf",sep=""),height=5)
MetageneBicorPlot(pbmc, grouping.var = "protocol", dims.eval = 1:20, 
    display.progress = FALSE)
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 1:6, 
    do.balanced = TRUE)
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 7:12, 
    do.balanced = TRUE)
DimHeatmap(object = pbmc, reduction.type = "cca", cells.use = 500, dim.use = 13:18, 
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


###### Louvain-jaccard clustering
pbmc <- FindClusters(object = pbmc, reduction.type = "cca.aligned", dims.use = 1:numCCs, 
    resolution=seq(0.1,3,by=0.1),save.SNN = TRUE)
save(pbmc,file=paste0("Monkey5Somatic-NoC1_CCA.Robj"))

print(c( length(unique(pbmc@meta.data$res.0.1)),length(unique(pbmc@meta.data$res.0.2)),length(unique(pbmc@meta.data$res.0.3)),length(unique(pbmc@meta.data$res.0.4)),length(unique(pbmc@meta.data$res.0.5)),length(unique(pbmc@meta.data$res.0.6)),length(unique(pbmc@meta.data$res.0.7)),length(unique(pbmc@meta.data$res.0.8)),length(unique(pbmc@meta.data$res.0.9)),length(unique(pbmc@meta.data$res.1)) ))
print(c( length(unique(pbmc@meta.data$res.1.1)),length(unique(pbmc@meta.data$res.1.2)),length(unique(pbmc@meta.data$res.1.3)),length(unique(pbmc@meta.data$res.1.4)),length(unique(pbmc@meta.data$res.1.5)),length(unique(pbmc@meta.data$res.1.6)),length(unique(pbmc@meta.data$res.1.7)),length(unique(pbmc@meta.data$res.1.8)),length(unique(pbmc@meta.data$res.1.9)),length(unique(pbmc@meta.data$res.2)) ))



res=paste0("res.",c(paste0("0.",1:9),1)) # res.0.8 for 10 clusters
plotlist=list()
pdf(paste(dgefile,"clusters.pdf",sep=""),height=7.5,width=9)
for(resi in 1:length(res)){
pbmc=SetAllIdent(pbmc,id=res[resi])
plotlist[[1]]=DimPlot(pbmc,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(pbmc,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=DimPlot(pbmc,reduction.use = "umap",cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[4]]=TSNEPlot(pbmc,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)

}
dev.off()

table(dge5@ident[names(pbmc@ident)],pbmc@meta.data$res.0.1)
table(dge5@ident[names(pbmc@ident)],pbmc@meta.data$res.0.2)
table(dge5@ident[names(pbmc@ident)],pbmc@meta.data$res.0.3)
table(dge5@ident[names(pbmc@ident)],pbmc@meta.data$res.0.4)
table(dge5@ident[names(pbmc@ident)],pbmc@meta.data$res.0.5)
                   0    1    2    3    4    5    6
  Macrophage       0    0    0    0    0   37    0
  Endothelial     14    0    2  174    0    0    1
  Leydig           6    1    0    0   41    0    0
  Pericyte         8    1  174    0    0    0    3
  Myoid          123  479    1    0    1    0    0
  IntProgenitor 1231   40    3    1    0    0   14



###### order clusters for each dataset
res="res.0.5"
resi=1
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

ordered="orderedclusters"

pbmc=AddMetaData(pbmc,cells.ident.ordered,ordered)
pbmc@meta.data[,ordered]=factor(pbmc@meta.data[,ordered])
pbmc=SetAllIdent(pbmc,ordered)
# save the pbmc file
save(pbmc,file=paste0("Monkey5Somatic-NoC1_CCA.Robj"))

# go above to re-plot PCA and tSNE



### compare with previous clusters
table(dge5celltype@ident,pbmc@ident[names(dge5celltype@ident)])
                   1    2    3    4    5    6    7
  Macrophage      37    0    0    0    0    0    0
  Endothelial      0  174    0    2    0   14    1
  Leydig           0    0   41    0    1    6    0
  Pericyte         0    0    0  174    1    8    3
  Myoid            0    0    1    1  479  123    0
  IntProgenitor    0    1    0    3   40 1231   14

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


### visualize known markers
knownmarkers=c("VIM",
  "CD163","S100A4","TYROBP","LYZ","RGS1",
  "VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN",
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","PDGFRA","PDGFRB","PDGFB","ALDH1",
  "DLK1","HSD17B3","STAR","CYP17A1","NR5A1","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "CLU","SOX9","AMH"  )
length(knownmarkers) # 35
knownmarkers[which(!(knownmarkers %in% rownames(pbmc@data)))] #[1] "ALDH1"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(pbmc@data))]
length(knownmarkers) # 32
pdf(paste0(dgefile,"knownmarkers_Feature.pdf"),height=15,width=21)
FeaturePlot(object = pbmc, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf(paste0(dgefile,"knownmarkers_Violin.pdf"),height=10,width=22)
VlnPlot(pbmc,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()

knownmarkers=c("NR2F2","TCF21","GATA4","NR5A1","EGR2","PDGFRA","PDGFA","PDGFRB","PDGFB",
    "THY1","ENG","VCAM1","ALCAM","CD44","NT5E","CASP3" )
# NGF1B: EGR2

length(knownmarkers) # 16
knownmarkers[which(!(knownmarkers %in% rownames(pbmc@data)))] 
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(pbmc@data))]
length(knownmarkers) # 14

pdf(paste0(dgefile,"knownmarkers_MSC_Feature.pdf"),height=12,width=12)
FeaturePlot(object = pbmc, features.plot = knownmarkers,col=c("grey80","red"),nCol=4)
dev.off()
pdf(paste0(dgefile,"knownmarkers_MSC_Violin.pdf"),height=8,width=13)
VlnPlot(pbmc,knownmarkers,cols.use=myBrewerPalette,nCol=4,point.size.use=-1)
dev.off()

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
knownmarkers[which(!(knownmarkers %in% rownames(pbmc@data)))] 
#[1] "FIGLA"  "FERD3L" "SCX"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(pbmc@data))]
length(knownmarkers) # 21

pdf("knownmarkers_TCF21paralogs_Feature.pdf",height=9,width=21)
FeaturePlot(object = pbmc, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf("knownmarkers_TCF21paralogs_Violin.pdf",height=6,width=22)
VlnPlot(pbmc,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()


###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(pbmc,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,"res.0.5order_markersall_mindiff0.2_logfc2fold_1.31.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
write.table(pbmc@ident,"MonkeyMerged5Somatic-NoC1_CCA_7ClustersID.txt",row.names=T,col.names=F,sep="\t",quote=F)

table(pbmc@ident)
table(markers$cluster)

   1    2    3    4    5    6    7 
  37  175   42  180  521 1382   18 
  110 132  97 107  21  27 126

markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf(paste0(dgefile,"markerstop.pdf"),height=9,width=15)
FeaturePlot(object = pbmc, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=5)
dev.off()

### PCA and tSNE plot
pbmc=SetAllIdent(pbmc,id="orig.ident")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1_Rep.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
pbmc=SetAllIdent(pbmc,id="protocol")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1_Subject.pdf"),width=10,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
pbmc=SetAllIdent(pbmc,id="orderedclusters")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1_orderedclusters.pdf"),width=9.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

###### Visualize individual batch and subject
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlims=ylims=list()
pos=c("bottomleft","topleft","topleft","bottomleft")
# pbmc4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-1.4,3.8),c(-1.4,3.8),c(-8,11),c(-28,35))
ylims[[1]]=list(c(-5.5,1.9),c(-3.5,4.5),c(-4.9,9.5),c(-45,35))

dim=list(pbmc@dr$cca.aligned@cell.embeddings[,1:2],pbmc@dr$cca.aligned@cell.embeddings[,c(1,3)],pbmc@dr$umap@cell.embeddings[,c(1,2)],pbmc@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

### plot PCs and tSNE for each batch using the other batches as background
pbmc=SetAllIdent(pbmc,id="orig.ident")
sets=levels(pbmc@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette[c(8:7,10,6:5,4:1)]

pdf(paste0(dgefile,"OrigSet_PCtSNE.pdf"),height=4.6,width=11.5)
par(mfrow=c(2,5),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:length(dim)){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(pbmc@ident)
names(ident)=names(pbmc@ident)
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
pbmc=SetAllIdent(pbmc,id="protocol")
sets=levels(pbmc@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

pdf(paste0(dgefile,"OrigSubject_PCtSNE.pdf"),height=4.6,width=6.9)
par(mfrow=c(2,3),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:length(dim)){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(pbmc@ident)
names(ident)=names(pbmc@ident)
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
pbmc=SetAllIdent(pbmc,id="orderedclusters")
table(pbmc@meta.data$protocol,pbmc@meta.data$orderedclusters)

  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(pbmc@meta.data$orig.ident,pbmc@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0(dgefile,"ncellspercluster_batch_top15CCs.txt"),quote=F,row.names=T,col.names=T,sep="\t")


pbmc=SetAllIdent(pbmc,id="orderedclusters")
sets=subject
which(rownames(data)!=names(pbmc@ident))


cols=myBrewerPalette1
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())


j=3
label="UMAP"
data=data_bg=as.data.frame(pbmc@dr$umap@cell.embeddings)
data$ident=pbmc@ident
data %>% dplyr::group_by(ident) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)) -> centers
label.size=6

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(pbmc@ident)[which(gsub("\\..*","",gsub("_.*","",names(pbmc@ident)))==set)]
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
data=data_bg=as.data.frame(pbmc@dr$tsne@cell.embeddings)
data$ident=pbmc@ident
data %>% dplyr::group_by(ident) %>% summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) -> centers
label.size=6

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(pbmc@ident)[which(gsub("\\..*","",gsub("_.*","",names(pbmc@ident)))==set)]
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
jpeg(file=paste0(dgefile,"OrigSubject_7clusters_",label,"_bg.jpeg"),res=300,height=1200,width=1900)
multiplot(plotset,cols = 3)
dev.off()
pdf(file=paste0(dgefile,"OrigSubject_7clusters_",label,"_bg.pdf"),height=6,width=9)
multiplot(plotset,cols = 3)
dev.off()


###### Rank correlation and Dissimilarity matrix for each normalized centroid using HVG
## order cell by cluster ID and randomly shuffle cells within each batch
levels=levels(pbmc@ident)
ident=factor(pbmc@ident,levels=levels)

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

write.table(genecountsall,paste0(dgefile,"Merged5Somatic-NoC1_CCA_7clusters_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)

cc=cor(as.matrix(genecountsall)[pbmc@var.genes,],method="spearman")
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

pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1_7clusters_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()

# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)):101]
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1_7clusters_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
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
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1NoC1_6celltypes_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()
# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)-1):101]
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1NoC1_6celltypes_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()


pdf(file=paste0("PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(pbmc, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()
