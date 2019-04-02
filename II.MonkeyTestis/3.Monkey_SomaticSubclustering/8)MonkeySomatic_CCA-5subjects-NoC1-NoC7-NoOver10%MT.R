### Somatic Subclustering for 5 monkey subjects corrected for batch effect by CCA after removing doublets cluster 1 and 7 and excluding cells with >=10% MT
# 3.23.2019 by Qianyi
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering 
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1 cells after removing C1 doublets
### Removed Cluster 7/7 of SCyte-Somatic doublets from somatic subclustering by CCA
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing C1 and C7 doublets
### Removed 239 cells with >=10% MT transcripts that we failed to exclude as we missed 7 MT gene symbols using Adrienne MT gene list
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7-No10pctMT cells 


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



library(Seurat)
# cowplot enables side-by-side ggplots
library(cowplot)


# setup Seurat objects 
### load object for CCA of 5 subjects before removing cells with >=10% MT
# CCA and somatic subclustering after removing Cluster 1/8 and 7/7 of doublets for somatic subclustering 
load(file="Monkey5Somatic-NoC1-NoC7_CCA.Robj")
testis1=testis
testis1 # 17469 genes across 2337 samples
table(testis@ident)
table(testis@ident[cells10pctMT])
#      Tcell  Macrophage Endothelial    Pericyte Subcluster2       Myoid 
#         20          17         170         131          43         516 
#          0           0           6          14           1          45 
#  ImmLeydig  DiffLeydig 
#       1396          44 
#        168           5

### remove cells with >=10% MT 
cells.use=names(testis1@ident)[which(!(names(testis1@ident)) %in% cells10pctMT)]
length(cells.use) # 2098
length(testis1@ident)-length(cells.use) # [1] 239
testis=SubsetData(testis1,cells.use=cells.use)
nCellperGene <- rowSums(testis@data>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # 17251

datalist=list()
for(i in 1:length(subject)){
### Extract Somatic-NoC1-NoC7-No10pctMT cells for each monkey subject
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
 17251 genes across 366 samples.

[[2]]
An object of class seurat in project SeuratProject 
 17251 genes across 152 samples.

[[3]]
An object of class seurat in project SeuratProject 
 17251 genes across 1264 samples.

[[4]]
An object of class seurat in project SeuratProject 
 17251 genes across 196 samples.

[[5]]
An object of class seurat in project SeuratProject 
 17251 genes across 120 samples.



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
length(hvg.union) # 3540

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
length(hvg.union) # 3540

for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(length(which(hvg.union %in% rownames(x = head(x = dge@hvg.info, n = 2000)))))
}
[1] 1323
[1] 1261
[1] 1300
[1] 1299
[1] 1263

# check the number of hvg overlapped between datasets
cc=matrix(0,length(subject),length(subject))
for(i in 1:length(subject)){
    for(j in 1:length(subject)){
        cc[i,j]=length(intersect(rownames(x = head(x = datalist[[i]]@hvg.info, n = 1000)),rownames(x = head(x = datalist[[j]]@hvg.info, n = 1000))))
    }
}
cc
     [,1] [,2] [,3] [,4] [,5]
[1,] 1000  221  284  288  246
[2,]  221 1000  216  225  211
[3,]  284  216 1000  277  243
[4,]  288  225  277 1000  255
[5,]  246  211  243  255 1000



### run multi CCA
testis <- RunMultiCCA(object.list = datalist, genes.use = hvg.union, num.ccs=20)
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA.Robj"))


dgefile=dgename="CCA5Subjects-Somatic-NoC1-NoC7/No10pctMT_"

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


###### Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
testislist=list()
for(numCCs in 11:15){
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
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA_",numCCs,"CCs.Robj"))

testislist[[numCCs]]=testis
}
print(c( length(unique(testis@meta.data$res.0.1)),length(unique(testis@meta.data$res.0.2)),length(unique(testis@meta.data$res.0.3)),length(unique(testis@meta.data$res.0.4)),length(unique(testis@meta.data$res.0.5)),length(unique(testis@meta.data$res.0.6)),length(unique(testis@meta.data$res.0.7)),length(unique(testis@meta.data$res.0.8)),length(unique(testis@meta.data$res.0.9)),length(unique(testis@meta.data$res.1)) ))
print(c( length(unique(testis@meta.data$res.1.1)),length(unique(testis@meta.data$res.1.2)),length(unique(testis@meta.data$res.1.3)),length(unique(testis@meta.data$res.1.4)),length(unique(testis@meta.data$res.1.5)),length(unique(testis@meta.data$res.1.6)),length(unique(testis@meta.data$res.1.7)),length(unique(testis@meta.data$res.1.8)),length(unique(testis@meta.data$res.1.9)),length(unique(testis@meta.data$res.2)) ))



for(numCCs in 11:15){
testis=testislist[[numCCs]]

res=paste0("res.",c(paste0("0.",1:9),1,paste0("1.",1:9))) # res.0.8 for 10 clusters
plotlist=list()
pdf(paste(dgefile,"clusters_",numCCs,"CCs.pdf",sep=""),height=7.5,width=9)
for(resi in 1:length(res)){
testis=SetAllIdent(testis,id=res[resi])
plotlist[[1]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=DimPlot(testis,reduction.use = "umap",cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[4]]=TSNEPlot(testis,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)
}
dev.off()

for(numCCs in 11:15){
testis=testislist[[numCCs]]
testis=AddMetaData(testis,testis1@ident,"CellType0")
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA_",numCCs,"CCs.Robj"))
testislist[[numCCs]]=testis
}

pdf(paste0(dgefile,"PCA_tSNE_UMAP_8celltypes0.pdf"),width=10.5,height=8)
for(numCCs in 11:15){
testis=testislist[[numCCs]]
testis=AddMetaData(testis,testis1@ident,"CellType0")
testis=SetAllIdent(testis,id="CellType0")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType0))
plotlist[[1]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette1,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette1,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=DimPlot(testis,reduction.use = "umap",cols.use=myBrewerPalette1,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[4]]=TSNEPlot(testis,colors.use=myBrewerPalette1,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)
}
dev.off()


### decided to use 14 CCs 
numCCs=14

print(table(testis1@ident[names(testis@ident)],testis@meta.data$res.0.5))
table(testis1@ident[names(testis@ident)],testis@meta.data$res.0.4)
table(testis1@ident[names(testis@ident)],testis@meta.data$res.0.5)
table(testis@meta.data$res.0.3,testis@meta.data$res.0.5)
table(testis@meta.data$res.0.4,testis@meta.data$res.0.5)


###### order clusters for each dataset
numCCs=14
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
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA.Robj"))

# go above to re-plot PCA and tSNE



### compare with previous clusters
table(testis1@ident[names(testis@ident)],testis@ident)
                 1    2    3    4    5    6
  Tcell         20    0    0    0    0    0
  Macrophage    17    0    0    0    0    0
  Endothelial    0  162    0    0    2    0
  Pericyte       0    0  116    0    1    0
  Subcluster2    0    3   33    0    6    0
  Myoid          0    0    0  355  113    3
  ImmLeydig      0   13    6    7 1201    1
  DiffLeydig     0    0    0    0    1   38


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
knownmarkers=c("NR2F2","MYH11","ACTA2","PTCH1","PTCH2","DCN","TCF21","PDGFRA","PDGFRB","PDGFB","DLK1","HSD17B3","STAR","CYP17A1","GATA4","SF1","NR5A1","EGR2","IGF1","IGFBP5","IGFBP3","INHBA","VIT","THY1","ENG","VCAM1","ALCAM","CD44","NT5E","PDGFA","CASP3","CLU","ALDH1A1","INSL3","PTGDS","CDC25C","ID4","SFRP1","KLF6","RSPO1","LCN2","ITLN1","RGS5","MCAM","CSPG4","MYL9","FRZB","CD36","NOTCH3","ADIRF","CRIP1","GATA4","TAGLN","VIM","CD163","S100A4","TYROBP","LYZ","RGS1","VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN","SOX9","AMH")
#NG2/CSPG4
#NGF1B/EGR2
#SF1/NR5A1 
#CD90/THY1
#CD105/ENG
#CD106/VCAM1
#CD166/ALCAM
#CD73/NT5E
#TAGLN
#STRO1
length(knownmarkers) # 67
knownmarkers[which(!(knownmarkers %in% rownames(testis@data)))]  # [1] "PDGFRA" "CASP3"  "ADIRF"  "CRIP1"  "AMH"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(testis@data))]
length(knownmarkers) # 62
pdf(paste0(dgefile,"knownmarkers_Feature_2pages.pdf"),height=27,width=21)
FeaturePlot(object = testis,reduction.use="umap", features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
FeaturePlot(object = testis, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf(paste0(dgefile,"knownmarkers_Violin.pdf"),height=20,width=22)
VlnPlot(testis,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
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

object=testis
dgename="markers/"

setname="knownmarkers"
genes=knownmarkers

setname="markers"
genes=markers


for(j in 1:length(genes)){
feature=features.plot=genes[j]

features.plot; dim.1 = 1; dim.2 = 2; cells.use = NULL; pt.size = 1;
pch.use = 16; reduction.use = "umap";
use.imputed = FALSE; no.axes = TRUE; no.legend = FALSE

cells.use <- set.ifnull(cells.use, colnames(object@data))
dim.code <- translate.dim.code(reduction.use)
dim.codes <- paste(dim.code, c(dim.1, dim.2), sep = "")
data.plot <- as.data.frame(object@dr$umap@cell.embeddings)

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



###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(testis,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,"res.0.5order_markersall_mindiff0.2_logfc2fold_2.22.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(testis@ident)
table(markers$cluster)
   1    2    3    4    5    6 
  37  170  174  516 1396   44 
  110 128 114  25  29  90 

markers=FindAllMarkers(testis,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0("monkey/NoC1NoC7No10pctMTCCA5Somatic_8celltypes_markersall_mindiff0.2_logfc2fold_3.27.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(testis@ident)
table(markers$cluster)
Tcell  Macrophage Endothelial  m-Pericyte  f-Pericyte       Myoid 
         20          17         178         116          39         362 
        121         173         123         149          72          31 
  ImmLeydig  DiffLeydig 
       1324          42
         40          89     


markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf(paste0(dgefile,"markerstop.pdf"),height=9,width=12)
FeaturePlot(object = testis, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=4)
dev.off()

### PCA and tSNE plot
testis=SetAllIdent(testis,id="orig.ident")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1-NoC7-No10pctMT_Rep.pdf"),width=10.5,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(8:7,10,6:5,4:1)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
testis=SetAllIdent(testis,id="protocol")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1-NoC7-No10pctMT_Subject.pdf"),width=10,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=T)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=T)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
testis=SetAllIdent(testis,id="orderedclusters")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged5Somatic-NoC1-NoC7-No10pctMT_orderedclusters.pdf"),width=9.5,height=8)
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



# 3.24.2019 Qianyi
###### label cell types 
testis=SetAllIdent(testis,id="orderedclusters")
table(testis@ident)
   1    2    3    4    5    6 
  37  178  155  362 1324   42 
table(testis1@ident[names(testis@ident)],testis@ident)
             
                 1    2    3    4    5    6
  Tcell         20    0    0    0    0    0
  Macrophage    17    0    0    0    0    0
  Endothelial    0  162    0    0    2    0
  Pericyte       0    0  116    0    1    0
  Subcluster2    0    3   33    0    6    0
  Myoid          0    0    0  355  113    3
  ImmLeydig      0   13    6    7 1201    1
  DiffLeydig     0    0    0    0    1   38

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
              37              178              155              362 
       ImmLeydig       DiffLeydig 
            1324               42 

testis=AddMetaData(testis,id,"CellType")
testis=SetAllIdent(testis,id="CellType")
testis@ident=id
table(testis@ident)
Macrophage/Tcell      Endothelial         Pericyte            Myoid 
              37              178              155              362 
       ImmLeydig       DiffLeydig 
            1324               42


### Visualize Somatic cell types
testis=SetAllIdent(testis,id="CellType3")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType3))
table(testis@ident)
      Tcell  Macrophage Endothelial  m-Pericyte  f-Pericyte       Myoid 
         20          17         178         116          39         362 
  ImmLeydig  DiffLeydig 
       1324          42 

save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA.Robj"))
write.table(testis@ident,"Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA_8celltypes_label.txt",row.names=T,col.names=F,quote=F,sep="\t")

# try to match color scheme for corresponding mouse cell types
#levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown")
#myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7)] # 1    
library(RColorBrewer)
col11 <- c(brewer.pal(12,"Paired")[c(12,11,10,6,8,7,9)],brewer.pal(8,"Set1")[8])  #used this 1 
col22 <- brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3,8)]
col33 <- c(brewer.pal(8,"Set1")[c(7,6,4)],gg_color_hue(5))
myBrewerPalette=col11
myBrewerPalette=col22

pdf(paste0(dgefile,"PCA_tSNE_UMAP_8celltypes1.pdf"),width=11,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=col11)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

pdf(paste0(dgefile,"PCA_tSNE_UMAP_8celltypes2.pdf"),width=11,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=col22)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

pdf(paste0(dgefile,"PCA_tSNE_UMAP_8celltypes3.pdf"),width=11,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col33)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col33)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col33)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=col33)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()



###### Rank correlation for all cells
nrho=cor(as.matrix(testis@data)[testis@var.genes,],method="spearman")
length(testis@var.genes) # 2490
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
for(i in 1:length(levels(cells.ident))){
  sidecol2[which(sidecol2[,1] == levels(cells.ident)[i]),2] <- i
}

library(RColorBrewer)
#myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
myBrewerPalette=col22
rlab2=rbind(rep("white",length(sidecol2[,1])),myBrewerPalette[as.numeric(sidecol2[,2])])
clab2=cbind(rlab2[2,],rlab2[1,])
colnames(clab2)=c("CellType","")
rownames(rlab2)=c("","CellType")

midrange=median(testcor[which(testcor!=1)])
maxrange=max(testcor[which(testcor!=1)])
maxrange2=maxrange
midrange2=maxrange/2
col.use2=redblue100[c(rep(1:50,each=round(midrange2*100)),rep(50:100,each=round((maxrange2-midrange2)*100)),rep(100,50*round((1-maxrange2)*100)))]
length(col.use2)

jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_8celltypes_RankCor_HVG_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(15,15))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_8celltypes_RankCor_HVG.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(8,8))                    # symm=F,symkey=F,symbreaks=F,
dev.off()

col.use2=redblue100[c(rep(c(41:100),each=10),rep(100,100000))]
length(col.use2)
jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_8celltypes_snn_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(15,15))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_8celltypes_snn.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(8,8))                    # symm=F,symkey=F,symbreaks=F,
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

write.table(genecountsall,paste0(dgefile,"Merged5Somatic-NoC1-NoC7-No10pctMT_CCA_6clusters_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)
write.table(genecountsall,paste0(dgefile,"Merged5Somatic-NoC1-NoC7-No10pctMT_CCA_8celltypes_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)

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

pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_6clusters_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()

# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)):101]
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_6clusters_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
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
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_8celltypes_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()
# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)-1):101]
pdf(file=paste0(dgefile,"5subjectsSomatic-NoC1-NoC7-No10pctMT_8celltypes_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()


### visualize dissimilarity matrix for cluster centroid
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
do=seriate(da,method="OLO")
print(get_order(do))
###### plot with seriation
pdf(file=paste0(dgefile,"Merged5Somatic-NoC1-NoC7-No10pctMT_Centroid_norm_Seriation_orderedclusters.pdf"))
 dissplot(da,method = NA, options = list(main = paste("Dissimilarity"),col=bluered(10)))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO"),col=bluered(10)))
 hmap(da,col=bluered(10)) # default method="OLO"
dev.off()



### Per-cell attributes
######### Analyze per-cell attributes 
library(ineq)
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


pdf(file=paste0("PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(testis, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="GiniAll", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()
