### Somatic SubClustering for CCA-merged 4 subjects (Huamn1,2,3,5) After Removing Doublets Clusters 5 and 1
# 2.1.2019 by Qianyi
# note: CCA was performed after removing 5/13 and 1/7 doublets
# after removing cells of Cluster 5/13 of SPG-myoid doublets for global clustering
# after removing cells of Cluster 1/7 of spermatid-myoid doublets for somatic subclustering by CCA
# re-do CCA and somatic subclustering -> 7 clusters -> labeled 7 cell types
#-> Final 7 somatic cell types for human testis

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


dgefile=dgename="CCA-4subjects1235NoC5-Somatic-NoC1/NoC5NoC1CCA_"

### load object for CCA of 4 subjects somatic before removing 1/7 doublets
# note: CCA was performed after removing 5/13 but before removing 1/7
# before removing cells of Cluster 5/13 of SPG-myoid doublets for global clustering
load(file="Human4-1235-NoC5Somatic_CCA.Robj")
pbmcC1=pbmc
# label cell types after removing cells of Cluster 1/7 of spermatid-myoid doublets for somatic subclustering by CCA
load(file="Human4-1235-NoC5Somatic_CCAcelltypeNoC1.Robj")
pbmc4=pbmc
pbmc # 35063 genes across 3722 samples.
# labeled with 6 somatic cell types
table(pbmc4@ident,pbmcC1@ident[names(pbmc4@ident)])
                 1    2    3    4    5    6    7
  Macrophage     0    0    0    0    0    0   65
  Monocyte       0    0    0    0    0  233    0
  Endothelial    0    0    0    0  293    0    0
  Leydig         0    0    0  521    0    0    0
  Myoid          0    0 2143    0    0    0    0
  Pericyte       0  467    0    0    0    0    0

library(Seurat)
# cowplot enables side-by-side ggplots
library(cowplot)

# setup Seurat objects 
### remove genes without expression in all cells
dgedata=pbmc4@raw.data[,names(pbmc4@ident)]
dim(dgedata)                       # [1] 35063  3722
nCellperGene <- rowSums(dgedata>0)
length(which(nCellperGene==0))     # 548
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgedata2=dgedata[genes.use,]
dim(dgedata2)                      # [1] 34515  3722
length(genes.use)                  # 34515
### since both count matrices have already filtered cells, we do no additional filtering here
datalist=list()
for(i in 1:length(subject)){
### Extract Somatic cells for each human subject
cells.use=names(pbmc4@ident)[which(pbmc4@meta.data$protocol == subject[i])]
dgedata=pbmc4@raw.data[genes.use,cells.use]
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
 34515 genes across 211 samples.

[[2]]
An object of class seurat in project SeuratProject 
 34515 genes across 600 samples.

[[3]]
An object of class seurat in project SeuratProject 
 34515 genes across 638 samples.

[[4]]
An object of class seurat in project SeuratProject 
 34515 genes across 2273 samples.

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
length(hvg.union) # 3198

# extract hvg that are present in every dataset
# skipped this
for(i in 1:length(subject)){
  dge=datalist[[i]]
  hvg.union=hvg.union[which(hvg.union %in% rownames(dge@data))]
  print(length(hvg.union))
}

# check the number of selected hvg that fall within top 2k hvg for each dataset
hvg.union=unique(unlist(hvglist))
length(hvg.union) # 3198

for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(length(which(hvg.union %in% rownames(x = head(x = dge@hvg.info, n = 2000)))))
}
[1] 1127
[1] 1144
[1] 1164
[1] 1140

# check the number of hvg overlapped between datasets
cc=matrix(0,length(subject),length(subject))
for(i in 1:length(subject)){
    for(j in 1:length(subject)){
        cc[i,j]=length(intersect(rownames(x = head(x = datalist[[i]]@hvg.info, n = 1000)),rownames(x = head(x = datalist[[j]]@hvg.info, n = 1000))))
    }
}
cc
     [,1] [,2] [,3] [,4]
[1,] 1000  247  222  135
[2,]  247 1000  276  178
[3,]  222  276 1000  177
[4,]  135  178  177 1000

### run multi CCA
pbmc <- RunMultiCCA(object.list = datalist, genes.use = hvg.union, num.ccs=20)
save(pbmc,file=paste0("Human4-1235-NoC5SomaticNoC1_CCA.Robj"))


dgefile=dgename="CCA-4subjects1235NoC5-Somatic-NoC1/NoC5NoC1CCA_"

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
#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
#`geom_smooth()` using method = 'loess' and formula 'y ~ x'

for(numCCs in c(7:15)){
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
pdf(paste(dgefile,numCCs,"CCs_dge_ACCs_top.pdf",sep=""),height=6,width=15)
plot_grid(p1, p2,p3,p4,p11,p22,p33,p44,ncol=4)
dev.off()

###### Now we can run a single integrated analysis on all cells!
pbmc <- RunTSNE(object = pbmc, reduction.use = "cca.aligned", dims.use = 1:numCCs, 
    do.fast = TRUE)
pbmc <- RunUMAP(pbmc, reduction.use = "cca.aligned", dims.use = 1:numCCs)

###### Louvain-jaccard clustering
pbmc <- FindClusters(object = pbmc, reduction.type = "cca.aligned", dims.use = 1:numCCs, 
    resolution=seq(0.1,3,by=0.1),save.SNN = TRUE)
save(pbmc,file=paste0("Human4-1235-NoC5SomaticNoC1_CCA_",numCCs,"CCs.Robj"))

res=paste0("res.",c(paste0("0.",1:9),1,paste0("1.",1:9))) # res.0.8 for 10 clusters
plotlist=list()
pdf(paste(dgefile,numCCs,"CCs_clusters.pdf",sep=""),height=7.5,width=9)
for(resi in 1:length(res)){
pbmc=SetAllIdent(pbmc,id=res[resi])
plotlist[[1]]=DimPlot(pbmc,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=DimPlot(pbmc,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(pbmc,reduction.use = "umap",cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[4]]=TSNEPlot(pbmc,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)
}
dev.off()

}


numCCslist=7:15
dgelist=list()
for(resi in 1:length(numCCslist)){
  numCCs=numCCslist[resi]
  load(file=paste0("Human4-1235-NoC5SomaticNoC1_CCA_",numCCs,"CCs.Robj"))
  dgelist[[resi]]=pbmc
}
for(resi in 1:length(numCCslist)){
  numCCs=numCCslist[resi]
  pbmc=dgelist[[resi]]
print(c( length(unique(pbmc@meta.data$res.0.1)),length(unique(pbmc@meta.data$res.0.2)),length(unique(pbmc@meta.data$res.0.3)),length(unique(pbmc@meta.data$res.0.4)),length(unique(pbmc@meta.data$res.0.5)),length(unique(pbmc@meta.data$res.0.6)),length(unique(pbmc@meta.data$res.0.7)),length(unique(pbmc@meta.data$res.0.8)),length(unique(pbmc@meta.data$res.0.9)),length(unique(pbmc@meta.data$res.1)) ))
print(c( length(unique(pbmc@meta.data$res.1.1)),length(unique(pbmc@meta.data$res.1.2)),length(unique(pbmc@meta.data$res.1.3)),length(unique(pbmc@meta.data$res.1.4)),length(unique(pbmc@meta.data$res.1.5)),length(unique(pbmc@meta.data$res.1.6)),length(unique(pbmc@meta.data$res.1.7)),length(unique(pbmc@meta.data$res.1.8)),length(unique(pbmc@meta.data$res.1.9)),length(unique(pbmc@meta.data$res.2)) ))
print(table(pbmc4@ident,pbmc@meta.data$res.1))
}

### Used top 15 CCs
numCCs=15
  load(file=paste0("Human4-1235-NoC5SomaticNoC1_CCA_15CCs.Robj"))
  pbmc15=pbmc

###### order clusters for each dataset
pbmc=pbmc15;
numCCs=15;res="res.0.5";resi=1

pbmc=SetAllIdent(pbmc,id=res[resi])
#pbmc <- BuildClusterTree(pbmc, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs)
#pbmc@meta.data[,res[resi]]=pbmc@ident
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
pdf(file=paste0(dgename,"Centroid_norm_Seriation_",numCCs,"CCs_",res[resi],".pdf"))
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

pdf(file=paste(dgename,"Centroid_RankedCorrelation_",numCCs,"CCs_",res[resi],".pdf",sep=""),height=5.5,width=5)
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
# save the pbmc file
save(pbmc,file=paste0("Human4-1235-NoC5SomaticNoC1_CCA_",numCCs,"CCs.Robj"))
save(pbmc,file=paste0("Human4-1235-NoC5SomaticNoC1_CCA.Robj"))
pbmc15=pbmc

# go above to re-plot PCA and tSNE

###### assign cell types based on markers for each ordered clusters in each dataset
numCCs=15;res="res.0.5order";resi=1
pbmc=pbmc15
pbmc=SetAllIdent(pbmc,id=res[resi])

table(pbmc4@ident[names(pbmc@ident)],pbmc@ident)
                 1    2    3    4    5    6    7
  Macrophage     0    0    0    0    0    0   65
  Monocyte       0    1    3    1    0  228    0
  Endothelial    0    0    0    0  293    0    0
  Leydig         0    1   22  497    1    0    0
  Myoid          0   35 2086   17    2    3    0
  Pericyte      79  373   10    4    1    0    0


markers=FindAllMarkers(pbmc,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,numCCs,"CCs_",res[resi],"order_markersall_mindiff0.2_logfc2fold_2.1.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(pbmc@ident)
table(markers$cluster)
   1    2    3    4    5    6    7 
  79  410 2121  519  297  231   65 
  55 100  32  42  56 161 117 

markers=FindAllMarkers(pbmc,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0("human/NoC5NoC1CCA4Somatic_7celltypes_markersall_mindiff0.2_logfc2fold_3.13.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(pbmc@ident)
table(markers$cluster)
     Tcell  Macrophage Endothelial   ImmLeydig       Myoid  m-Pericyte 
         65         231         297         519        2121         410 
        117         161          56          42          32         100 
 f-Pericyte 
         79      
         55 


markers=FindMarkers(pbmc,3,4,test.use="bimod",logfc.threshold = log(1.6),min.diff.pct=0.1,do.print = TRUE)
write.table(markers,paste0(dgefile,numCCs,"CCs_",res[resi],"order_markers3Vs4_mindiff0.1_logfc1.6fold_2.1.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
dim(markers[which(markers[,2]<0),])
dim(markers[which(markers[,2]>0),])

markers=FindMarkers(pbmc,c(3,4),only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,numCCs,"CCs_",res[resi],"order_markers34VsOthers_mindiff0.2_logfc2fold_2.1.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")

markers=FindMarkers(pbmc,1,2,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,numCCs,"CCs_",res[resi],"order_markers1Vs2_mindiff0.2_logfc2fold_2.1.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
dim(markers[which(markers[,2]<0),]) # 42
dim(markers[which(markers[,2]>0),]) # 58

markers=FindMarkers(pbmc,c(1,2),only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,numCCs,"CCs_",res[resi],"order_markers12VsOthers_mindiff0.2_logfc2fold_2.1.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
dim(markers)


markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf(paste0(dgefile,"markerstop.pdf"),height=9,width=15)
FeaturePlot(object = pbmc, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=5)
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
knownmarkers[which(!(knownmarkers %in% rownames(pbmc@data)))]  # [1] "RSPO1" "LCN2"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(pbmc@data))]
length(knownmarkers) # 65
pdf(paste0(dgefile,"knownmarkers_Feature_2pages.pdf"),height=30,width=21)
FeaturePlot(object = pbmc, reduction.use="umap",features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
FeaturePlot(object = pbmc, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf(paste0(dgefile,"knownmarkers_Violin.pdf"),height=22,width=22)
VlnPlot(pbmc,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
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

object=pbmc
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



### PCA and tSNE plot
pbmc=SetAllIdent(pbmc,id="orig.ident")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-1235-NoC5Somatic_Rep.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
pbmc=SetAllIdent(pbmc,id="protocol")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-1235-NoC5Somatic_Subject.pdf"),width=10,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
pbmc=SetAllIdent(pbmc,id=res[resi])
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-1235-NoC5Somatic_orderedclusters.pdf"),width=9.5,height=8)
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
pos=c("bottomright","bottomleft","topleft","bottomleft")
# pbmc4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-1.5,4.5),c(-1.5,4.5),c(-6,7.5),c(-35,38))
ylims[[1]]=list(c(-4.2,1.9),c(-4.2,2.8),c(-5.2,9.8),c(-38,32))

dim=list(pbmc@dr$cca.aligned@cell.embeddings[,1:2],pbmc@dr$cca.aligned@cell.embeddings[,c(1,3)],pbmc@dr$umap@cell.embeddings[,c(1,2)],pbmc@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

### plot PCs and tSNE for each batch using the other batches as background
pbmc=SetAllIdent(pbmc,id="orig.ident")
sets=levels(pbmc@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette[c(5:9,2:1,4:3,10)]

pdf(paste0(dgefile,"OrigSet_ACC_UMAP_tSNE.pdf"),height=4.6,width=11.5)
par(mfrow=c(2,5),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(pbmc@ident)
names(ident)=names(pbmc@ident)
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.4,col=c("grey70",cols[seti])[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
}
dev.off()

### plot PCs and tSNE for each subject using the other subjects as background
pbmc=SetAllIdent(pbmc,id="protocol")
sets=levels(pbmc@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

pdf(paste0(dgefile,"OrigSubject_ACC_UMAP_tSNE.pdf"),height=4.6,width=4.6)
par(mfrow=c(2,2),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
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
}
dev.off()

### plot PCs and tSNE for each batch using the other batches as background
### visualize 12 clusters
pbmc=SetAllIdent(pbmc,id=res)
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


pbmc=SetAllIdent(pbmc,id=res)
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
source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_pbmc/Rcode_multiplot.R")
jpeg(file=paste0(dgefile,"OrigSubject_7clusters_",label,"_bg.jpeg"),res=300,height=1500,width=1500)
multiplot(plotset,cols = 2)
dev.off()
pdf(file=paste0(dgefile,"OrigSubject_7clusters_",label,"_bg.pdf"),height=6,width=6)
multiplot(plotset,cols = 2)
dev.off()


###### Pathway analysis by Gorilla
write.table(rownames(pbmc@data),paste0(dgefile,"allgenes_background.txt"),row.names=F,col.names=F,quote=F,sep="\t")

## Generate 1 target lists of markers for each cluster
  markers=read.table(paste0(dgefile,"15CCs_res.0.5orderorder_markersall_mindiff0.2_logfc2fold_2.1.2019.txt"),row.names=1,header=T)
for(i in 1:length(unique(markers$cluster))){
  write.table(rownames(markers)[which(markers$cluster==i,)],paste0(dgefile,"markers_",i,".txt"),row.names=F,col.names=F,quote=F,sep="\t")
}

cbl-gorilla.cs.technion.ac.il
Step1: species: Homo sapiens
Step2: mode: choose two unranked lists of genes (target and background)
Step3: list of genes; each line is a gene name
    Target set: paste genes or upload choose file: GOinput_celltype1.txt
    background set: upload choose file: background.txt
Step4: choose ontology: all
Advanced: analysis name
      email: qzm@umich.edu
          check Output results in Microsoft Excel format
          check Show output also in REViGO
Click: serach enriched GO terms

save webpage in bookmark for a month
or save webpage locally
copy and paste into Excel for all three tabs: Process, Function, Component
-> GOrilla10celltypes_8.14.2017.xlsx
paste total number of genes as Info sheet
add a column of GO to indicate Process, Function, Component copied from each tab

Sort: by FDR p-value from smallest to largest
Add a column of -log10(FDR)

Filter: Remove pathways that are too abundant or too rare
    Unmerged p-value color scale
    Select Column E: Text to Column: Delimited/Other: , Finish
    Select all table
    Data/Filter: B (the total number of Genes associated with a GO term): Number Filters: <=2000
    Number Filter: Between 50-500
    Select All, copy and paste into
-> GOrilla10celltypes_8.14.2017.xlsx

Plot name and -log10(FDR), color by corresponding RGB colors of RColorBrewer
http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
http://www.color-hex.com/color/cab2d6


### end

# 2.26.2019 Qianyi
###### label cell types 
pbmc=SetAllIdent(pbmc,id="res.0.5order")
table(pbmc@ident)
   1    2    3    4    5    6    7
  79  410 2121  519  297  231   65

### label cell types
id=as.numeric(pbmc@ident)
names(id)=names(pbmc@ident)
celltype=c("f-Pericyte","m-Pericyte","Myoid","ImmLeydig","Endothelial","Macrophage","Tcell")
for(i in 1:length(celltype)){
id[which(id==i)]<-celltype[i]
}
celltype=rev(celltype)
id=factor(id,levels=celltype,order=T)
table(id)
      Tcell  Macrophage Endothelial   ImmLeydig       Myoid  m-Pericyte
         65         231         297         519        2121         410
 f-Pericyte
         79

pbmc=AddMetaData(pbmc,id,"CellType")
pbmc=SetAllIdent(pbmc,id="CellType")
pbmc@ident=id
table(pbmc@ident)
      Tcell  Macrophage Endothelial   ImmLeydig       Myoid  m-Pericyte
         65         231         297         519        2121         410
 f-Pericyte
         79

save(pbmc,file=paste0("Human4-1235-NoC5SomaticNoC1_CCA.Robj"))
write.table(pbmc@ident,"Human4-1235-NoC5SomaticNoC1_CCA_7celltypes_label.txt",row.names=T,col.names=F,quote=F,sep="\t")

# try to match color scheme for corresponding mouse cell types
#levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown")
#myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7)] # 1

celltype=c("Tcell","Macrophage","Endothelial","ImmLeydig","Myoid","m-Pericyte","f-Pericyte")
library(RColorBrewer)
myBrewerPalette1 <- c(brewer.pal(7,"Set1"))
myBrewerPalette2 <- c(brewer.pal(12,"Paired"))[c(12,11,10,9,7,6,8)] 

pdf(paste0(dgefile,"PCA_tSNE_UMAP_7celltypes.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

pdf(paste0(dgefile,"PCA_tSNE_UMAP_7celltypes1.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette1)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

pdf(paste0(dgefile,"PCA_tSNE_UMAP_7celltypes2.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette2)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette2)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette2)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette2)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()


# 3.20.2019 Qianyi change levels order
levels=c("Tcell","Macrophage","Endothelial","m-Pericyte","f-Pericyte","Myoid","ImmLeydig")
pbmc@ident=factor(pbmc@ident,levels=levels)
library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,8,7,9)]  
library(RColorBrewer)
col11 <- c(brewer.pal(12,"Paired")[c(12,11,10,6,8,7,9)])   
col22 <- brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3)]  #used this 1
col33 <- c(brewer.pal(8,"Set1")[c(7,6,4)],gg_color_hue(4))

pdf(paste0(dgefile,"PCA_tSNE_UMAP_7celltypes11.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

pdf(paste0(dgefile,"PCA_tSNE_UMAP_7celltypes22.pdf"),width=10.5,height=8)
plot2=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot3=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=col11)
plot4=DimPlot(pbmc,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot5=TSNEPlot(pbmc,do.return=TRUE,pt.size = .8,do.label=F,colors.use=col22)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
plot2=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot3=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot4=DimPlot(pbmc,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot5=DimPlot(pbmc,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

myBrewerPalette=col22


###### Rank correlation for all cells
nrho=cor(as.matrix(pbmc@data)[pbmc@var.genes,],method="spearman")
length(pbmc@var.genes) # 3518
testcor=nrho

###### Jaccard distance
testcor=as.matrix(pbmc@snn)


### order cells by clusters 
res=res;j=1;resi=1;
   pbmc=SetAllIdent(pbmc,id=res)
   print(length(unique(pbmc@ident)))
   TSNEPlot(pbmc)

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

rlab2=rbind(rep("white",length(sidecol2[,1])),myBrewerPalette[as.numeric(gsub(".*-","",sidecol2[,2]))])
clab2=cbind(rlab2[2,],rlab2[1,])
colnames(clab2)=c("CellType","")
rownames(rlab2)=c("","CellType")

midrange=median(testcor[which(testcor!=1)])
maxrange=max(testcor[which(testcor!=1)])
maxrange2=maxrange
midrange2=maxrange/2
col.use2=redblue100[c(rep(1:50,each=round(midrange2*100)),rep(50:100,each=round((maxrange2-midrange2)*100)),rep(100,50*round((1-maxrange2)*100)))]
length(col.use2)

jpeg(file=paste(dgefile,"4subjectsSomatic-NoC5-NoC1_7celltypes_RankCor_HVG_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(15,15))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgefile,"4subjectsSomatic-NoC5-NoC1_7celltypes_RankCor_HVG.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(8,8))                    # symm=F,symkey=F,symbreaks=F,
dev.off()

col.use2=redblue100[c(rep(c(41:100),each=10),rep(100,100000))]
length(col.use2)
jpeg(file=paste(dgefile,"4subjectsSomatic-NoC5-NoC1_7celltypes_snn_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(15,15))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgefile,"4subjectsSomatic-NoC5-NoC1_7celltypes_snn.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(8,8))                    # symm=F,symkey=F,symbreaks=F,
dev.off()


###### Rank correlation and Dissimilarity matrix for each normalized centroid using HVG
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

write.table(genecountsall,paste0(dgefile,"CCA4-1235-NoC5SomaticNoC1_",numCCs,"CCs_7clusters_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)
write.table(genecountsall,paste0(dgefile,"CCA4-1235-NoC5SomaticNoC1_7celltypes_Centroid.txt"),quote=F,sep="\t",row.names=T,col.names=T)

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

pdf(file=paste0(dgefile,"4subjectsSomatic_7clusters_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()

# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)):101]
pdf(file=paste0(dgefile,"4subjectsSomatic_7clusters_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
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
pdf(file=paste0(dgefile,"4subjectsSomaticNoC1_7celltypes_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()
# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)-1):101]
pdf(file=paste0(dgefile,"4subjectsSomaticNoC1_7celltypes_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()


### Per-cell attributes
######### Analyze per-cell attributes 
###### add %ChrX and %ChrY genes to the pbmc datainfo
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

chrx.genes <- x[which(x %in% rownames(pbmc@data))]
chry.genes <- y[which(y %in% rownames(pbmc@data))]
autosome.genes <- autosome[which(autosome[,1] %in% rownames(pbmc@data)),]
length(chrx.genes)  # 1562
length(chry.genes)  # 133
dim(autosome.genes) # 35769  2
table(autosome.genes[,2])

# %x expression
percent.x <- colSums(expm1(pbmc@data[chrx.genes, ]))/colSums(expm1(pbmc@data))
percent.y <- colSums(expm1(pbmc@data[chry.genes, ]))/colSums(expm1(pbmc@data))
percent.autosome <- colSums(expm1(pbmc@data[autosome.genes[,1], ]))/colSums(expm1(pbmc@data))
mito.genes <- grep(pattern = "^MT-|^mt-", x = rownames(x = dge@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
# named integer(0)

# add to pbmc
pbmc <- AddMetaData(pbmc, percent.x, "percent.x")
pbmc <- AddMetaData(pbmc, percent.y, "percent.y")
pbmc <- AddMetaData(pbmc, percent.autosome, "percent.autosome")
pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")

###### calculate Gini index
library(ineq)
#for(i in 1:(length(dataset)+1)){
#  pbmc=pbmclist[[i]]
GiniAll=apply( pbmc@raw.data,2,function(x) ineq(x,type=c("Gini")) )
GiniNon0=apply( pbmc@raw.data,2,function(x) ineq(x[which(x!=0)],type=c("Gini")) )

pbmc=AddMetaData(pbmc,GiniAll,"GiniAll")
pbmc=AddMetaData(pbmc,GiniNon0,"GiniNon0")

save(pbmc,file=paste0("Human4-1235-NoC5SomaticNoC1_CCA_",numCCs,"CCs.Robj"))
save(pbmc,file=paste0("Human4-1235-NoC5SomaticNoC1_CCA.Robj"))
pbmc15=pbmc


###### Correlation between Gini and nUMI
datainfo=pbmc@meta.data
datainfo_bg=datainfo[,! names(datainfo) %in% "CellType"]
myBrewerPalette=brewer.pal(7,"Set1")

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
pdf(file=paste(dgefile,"GiniVsnUMI2.pdf",sep=""),height=3.5,width=7)
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


pdf(file=paste0("PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(pbmc, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="percent.x", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="percent.y", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="percent.autosome", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="GiniAll", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(pbmc, features.plot="GiniNon0", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()

