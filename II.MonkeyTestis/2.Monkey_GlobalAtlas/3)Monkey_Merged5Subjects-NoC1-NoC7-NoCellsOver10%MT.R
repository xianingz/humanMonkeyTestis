### Global Clustering for Directly-merged all 5 monkey subjects after removing two doublets clusters and excluding cells with >=10% MT using correct 37 MT gene list
# by Qianyi on 3.26.2019
### directly-merged 5 monkey subjects on 12/18/2018
#-> 12 clusters
#Cluster 1-2: somatic cells
#Cluster 3-4: SPG
#Cluster 5-8: Spermatocytes
#Cluster 9-10: Round spermatids
#Cluster 11-12: Elongating spermatids
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering 
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing doublets
### Removed Cluster 7/7 of SCyte-Somatic doublets from somatic subclustering by CCA
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing C1 and C7 doublets
### clustering of 3-species somatic cells merged by CCA -> Separated T-cells and Macrophages 
### Pericyte subclustering -> separated two pericyte subclusters
### directly-merge 5 monkey subjects after removing two doublets clusters (1/8 and 7/7) and identified cells failed to be removed with >=10% MT genes
### Re-did CCA and Subclustering for Somatic-NoC1-NoC7-NoCellsOver10%MT cells after removing cells with >=10% MT genes
### directly-merge 5 monkey subjects after removing two doublets clusters (1/8 and 7/7) and cells with >=10% MT genes


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



### 3.26.2019 by Qianyi 

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

### load object before removing doublets
load(file="MonkeyMerged5.Robj")
dge5=dge

### load object for merged 5 Monkey subjects after removing doublets and before removing cells with >=10% MT genes
load(file="MonkeyMerged5-NoC1-NoC7.Robj")
dge=dgeall
dge # 23089 genes across 22015 samples.

### remove cells with >=10% MT
dge=dgeall
cells10pctMT=rownames(dge@meta.data)[which(dge@meta.data$percent.mito>=0.1)]
length(cells10pctMT)  # 441
cells.use=rownames(dge@meta.data)[which(dge@meta.data$percent.mito<0.1)]

dgedata2=dgedata[,cells.use]
dim(dgedata2)                   #[1] 23101 21574
nCellperGene <- rowSums(dgedata2>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use)               #[1] 23029
dgedata2=dgedata[genes.use,cells.use]
dim(dgedata2)                   #[1] 23029 21574

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project="MonkeyDirectlyMerged5-NoC1-NoC7-No10pctMT", min.cells=1, min.genes=1)
dge # 23029 genes across 21574 samples. 

### add %ChrX and %ChrY genes to the dge datainfo
#note: convert Monkey Ensembl IDs to Human gene symbols after calculating %ChrX,Y,MT,autosome transcripts
#not right; must have already converted some Ensembl ID to gene symbols????
x=read.table("../monkeyChrXgenes",stringsAsFactors=FALSE)[,1]
y=read.table("../monkeyChrYgenes",stringsAsFactors=FALSE)[,1]
MT=read.table("../monkeyChrMTgenes",stringsAsFactors=FALSE)[,1]
MT2 <- c("ENSMMUG00000028704", "ENSMMUG00000028703", "ENSMMUG00000028702", "ENSMMUG00000028701", "ENSMMUG00000028700", "ENSMMUG00000028699", "ENSMMUG00000028698", "ENSMMUG00000028697", "ENSMMUG00000028696", "ENSMMUG00000028694", "ENSMMUG00000028693", "ENSMMUG00000028692", "ENSMMUG00000028691", "ENSMMUG00000028690", "ENSMMUG00000028688", "ENSMMUG00000028687", "ENSMMUG00000028686", "ENSMMUG00000028685", "ENSMMUG00000028684", "ENSMMUG00000028683", "ENSMMUG00000028681", "ENSMMUG00000028679", "ENSMMUG00000028678", "ENSMMUG00000028676", "ENSMMUG00000028675", "ENSMMUG00000028674", "ENSMMUG00000028671", "ENSMMUG00000028670", "ENSMMUG00000028669", "ENSMMUG00000028668", "ENSMMUG00000028695", "MT-CO1", "MT-CO3", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6")
MT[which(!(MT %in% MT2))] # [1] "ATP6" "ATP8" "COX2" "CYTB" "ND1"  "ND2"  "ND4L"
MT2[which(!(MT2 %in% MT))] # [1] "ENSMMUG00000028699" "ENSMMUG00000028686" "ENSMMUG00000028684" [4] "ENSMMUG00000028683" "ENSMMUG00000028678" "ENSMMUG00000028670" [7] "ENSMMUG00000028695"
length(MT[which(!(MT %in% MT2))]) # 7
dge@data[MT[which(!(MT %in% MT2))],1:2] # yes, found, use MT
dge@data[MT2[which(!(MT2 %in% MT))],1:2] # nothing
MT[which(!(MT %in% rownames(dge@data)))] # [1] "ENSMMUG00000028679" "ENSMMUG00000028681"
#note: the gene names in gene expression matrix were already converted to Monkey gene symbols from Monkey Ensembl IDs
autosome=read.table("../monkeyautosomegenes",stringsAsFactors=FALSE)
length(x) # 1197
length(y) # 60
length(MT) # 37
dim(autosome) # 25352
table(autosome[,2])
1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2566 1542 1560 1434 1141 1293 1794 1021 1039 1267 1416  905  956 1523 1020 1380 
  17   18   19   20 
 551  424 1522  998 

chrx.genes <- x[which(x %in% rownames(dge@data))]
chry.genes <- y[which(y %in% rownames(dge@data))]
mito.genes=MT[which(MT %in% rownames(dge@data))]
mito.genes2=MT2[which(MT2 %in% rownames(dge@data))]
autosome.genes <- autosome[which(autosome[,1] %in% rownames(dge@data)),]
length(chrx.genes)  # 956
length(chry.genes)  # 48
length(mito.genes)  # 35
length(mito.genes2) # 28
dim(autosome.genes) # 21990  2
table(autosome.genes[,2])
   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2263 1338 1356 1266 1001 1109 1517  857  902 1114 1262  761  812 1293  874 1197 
  17   18   19   20 
 463  355 1365  885 

# %x expression
percent.x <- colSums(expm1(dge@data[chrx.genes, ]))/colSums(expm1(dge@data))
percent.y <- colSums(expm1(dge@data[chry.genes, ]))/colSums(expm1(dge@data))
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
percent.mito2 <- Matrix::colSums(dge@raw.data[mito.genes2, ])/Matrix::colSums(dge@raw.data)
percent.autosome <- colSums(expm1(dge@data[autosome.genes[,1], ]))/colSums(expm1(dge@data))

which(percent.mito>=0.1) # named integer(0)

# add to dge
dge <- AddMetaData(dge, percent.x, "percent.x")
dge <- AddMetaData(dge, percent.y, "percent.y")
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
dge <- AddMetaData(dge, percent.autosome, "percent.autosome")

### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("Monkey_Merged5_HVG0.1.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.1, x.high.cutoff = 10, y.cutoff = 0.1, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # 3284
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 23029 genes across 21574 samples.          
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
   2013    2742    4725    6533    5561
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
dgeall=dge
save(dgeall,file="MonkeyMerged5-NoC1-NoC7-No10pctMT.Robj")
table(gsub("_.*","",names(dge@ident)))
Monkey1.1 Monkey1.2   Monkey2 Monkey3.1 Monkey3.2 Monkey4.1 Monkey4.2 Monkey5.1 
     1217       796      2742      2314      2411      2915      3618      2925 
Monkey5.2 
     2636 

### PCA
print(Sys.time())   #[1] "2019-03-26 15:13:57 EDT" 
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)
dgeall=dge

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
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,3,by=0.1), save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])
Sys.time() #[1] "2019-03-27 09:02:08 EDT"
dgeall=dge
save(dgeall,file="MonkeyMerged5-NoC1-NoC7-No10pctMT.Robj")


### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.6)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
print(c( length(unique(dge@meta.data$res.1.3)),length(unique(dge@meta.data$res.1.4)),length(unique(dge@meta.data$res.1.5)),length(unique(dge@meta.data$res.1.6)),length(unique(dge@meta.data$res.1.7)),length(unique(dge@meta.data$res.1.8)),length(unique(dge@meta.data$res.1.9)),length(unique(dge@meta.data$res.2)),length(unique(dge@meta.data$res.2.1)),length(unique(dge@meta.data$res.2.2)),length(unique(dge@meta.data$res.2.3)),length(unique(dge@meta.data$res.2.4)),length(unique(dge@meta.data$res.2.5)),length(unique(dge@meta.data$res.2.6)),length(unique(dge@meta.data$res.2.7)),length(unique(dge@meta.data$res.2.8)),length(unique(dge@meta.data$res.2.9)),length(unique(dge@meta.data$res.3)) ))
# 5  9 11 12 13 13 16 17 19 19 19 20

table(dge@meta.data$res.0.4,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.5,dge@meta.data$res.0.6)
table(dge@meta.data$res.0.5,dge@meta.data$res.0.7)
table(dge@meta.data$res.0.8,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.9,dge@meta.data$res.0.5)
table(dge@meta.data$res.1,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.4,dge5@meta.data[rownames(dge@meta.data),]$CellType2)
table(dge@meta.data$res.0.5,dge5@meta.data[rownames(dge@meta.data),]$CellType2)
table(dge@meta.data$res.0.6,dge5@meta.data[rownames(dge@meta.data),]$CellType2)
table(dge@meta.data$res.0.7,dge5@meta.data[rownames(dge@meta.data),]$CellType2)
table(dge@meta.data$res.0.8,dge5@meta.data[rownames(dge@meta.data),]$CellType2)
table(dge@meta.data$res.0.9,dge5@meta.data[rownames(dge@meta.data),]$CellType2)
table(dge@meta.data$res.1,dge5@meta.data[rownames(dge@meta.data),]$CellType2)


### decide to use 13 clusters, double-check the number of clusters
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
pdf(file=paste0("Merged5subjects-NoC1-NoC7-No10pctMT_Centroid_norm_Seriation_",res,".pdf"))
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
dgeall=dge
save(dgeall,file="MonkeyMerged5-NoC1-NoC7-No10pctMT.Robj")


## double-check if I need to reverse the cluster ID orders
### use PRM2 as marker for the last cluster (elongating)
dge=SetAllIdent(dge,id="res.0.5order")
pdf("clusters_ordered1_PRM2.pdf",height=10,width=10)
VlnPlot(dge,c("ID4","PIWIL4","TCF3","FGFR3","MORC1","ZCWPW1","ACRV1","PRM2"),cols.use=myBrewerPalette,nCol=2,point.size.use=-1)
dev.off()
pdf("clusters_ordered1_BEND2.pdf",height=10,width=5)
VlnPlot(dge,c("NR2F2","BEND2","SSX1","PRDM9"),cols.use=myBrewerPalette,nCol=1,point.size.use=-1)
dev.off()


## after reverse cluster order, repeat the above plots
### compare with previous 13 clusters of merged 4 subjects after removing cluster 5, before removing cluster 1
table(dge@ident,dge5@meta.data[rownames(dge@meta.data),]$CellType2)
table(dge@ident,dge5@meta.data[rownames(dge@meta.data),]$res.0.5order)
Somatic Spermatogonia Spermatocyte RoundSpermatid Elongating
  1     1667             0            0              0          0
  2      431             1            0              0          0
  3        0           463            1              0          0
  4        0           687            4              0          0
  5        0             0         2003              0          0
  6        0             0         1057             10          0
  7        0             1          648              6          6
  8        0             0            6            567          0
  9        0             0          217            998          0
  10       0             0           23           1097         37
  11       0             0            7             73       3115
  12       0             0            5              0       3753
  13       0             0            0              0       4691
        1    2    3    4    5    6    7    8    9   10   11   12
  1  1662    5    0    0    0    0    0    0    0    0    0    0
  2     8  423    1    0    0    0    0    0    0    0    0    0
  3     0    0  461    2    0    0    1    0    0    0    0    0
  4     0    0    4  683    2    0    2    0    0    0    0    0
  5     0    0    0    0 1519  419   65    0    0    0    0    0
  6     0    0    0    0    0 1047   10   10    0    0    0    0
  7     0    0    0    1    0    1  647    1    5    1    3    2
  8     0    0    0    0    0    2    4  567    0    0    0    0
  9     0    0    0    0    0    0  217  679  319    0    0    0
  10    0    0    0    0    0    0   23    1 1096   37    0    0
  11    0    0    0    0    0    0    7    0   73 2375  740    0
  12    0    0    0    0    0    0    5    0    0   18 3538  197
  13    0    0    0    0    0    0    0    0    0    0  312 4379

# change ordered0 in the file name to order1
knownmarkers=c("VIM","CD163","S100A4","TYROBP","LYZ","RGS1","VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN","MYH11","ACTA2","DCN","DLK1","NR2F2","CYP17A1","CLU","ID4","FGFR3", "MORC1", "PIWIL4", "TCF3","ZCWPW1","PRM2","ACRV1")
length(knownmarkers) # 26
pdf("knownmarkers.pdf",height=8,width=21)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
knownmarkers=c("ID4","PIWIL4","TCF3","FGFR3","MORC1","BEND2","ZCWPW1","SSX1","PRDM9","ACRV1","PRM2")
length(knownmarkers) # 11
pdf("knownmarkers2.pdf",height=8,width=12)
FeaturePlot(object = dge,reduction.use="umap", features.plot = knownmarkers,col=c("lightblue","red"),nCol=4)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("lightblue","red"),nCol=4)
dev.off()


###### assign cell types based on markers for each ordered clusters in each dataset
dge=SetAllIdent(dge,id="res.0.5order")
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Merged5subjects-NoC1-NoC7-No10pctMT_res.0.5order_markersall_mindiff0.2_logfc2fold_3.28.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)

   1    2    3    4    5    6    7    8    9   10   11   12   13 
1667  432  464  691 2003 1067  661  573 1215 1157 3195 3758 4691 
323 241 859 679 848 594  26 355 214  78  17   2  77

markers67=FindMarkers(dge,6,7,only.pos=FALSE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers67,"Merged5subjects-NoC1-NoC7-No10pctMT_res.0.5order_markers6Vs7_mindiff0.2_logfc2fold_3.28.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
markers87=FindMarkers(dge,8,7,only.pos=FALSE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers87,"Merged5subjects-NoC1-NoC7-No10pctMT_res.0.5order_markers8Vs7_mindiff0.2_logfc2fold_3.28.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
markers68=FindMarkers(dge,6,8,only.pos=FALSE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers68,"Merged5subjects-NoC1-NoC7-No10pctMT_res.0.5order_markers6Vs8_mindiff0.2_logfc2fold_3.28.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
length(markers67[which(markers67$avg_logFC>0),1]) # 4 
length(markers67[which(markers67$avg_logFC<0),1]) # 12
length(markers87[which(markers87$avg_logFC>0),1]) # 96
length(markers87[which(markers87$avg_logFC<0),1]) # 19
length(markers68[which(markers68$avg_logFC>0),1]) # 182
length(markers68[which(markers68$avg_logFC<0),1]) # 284


markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf("markerstop.pdf",height=12,width=21)
FeaturePlot(object = dge, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=7)
FeaturePlot(object = dge, reduction.use="umap", features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=7)
dev.off()
markers %>% group_by(cluster) %>% top_n(1, avg_logFC)  -> top1markers
pdf("markerstop.pdf",height=9,width=15)
FeaturePlot(object = dge,reduction.use="umap", features.plot = top1markers$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=5)
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
genes=top1markers$gene

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
# flip PC2, UMAP to be consistent with orientation of human testis plots
dge@dr$pca@cell.embeddings[,2]=-dge@dr$pca@cell.embeddings[,2]
dge@dr$umap@cell.embeddings[,1]=-dge@dr$umap@cell.embeddings[,1]

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
pos=c("topleft","bottomleft","bottomleft","bottomleft")
# dge5 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-16,70),c(-16,70),c(-43,49),c(-13,8))
ylims[[1]]=list(c(-28,48),c(-46,24),c(-45,42),c(-15,9))

dim=list(dge@dr$pca@cell.embeddings[,1:2],dge@dr$pca@cell.embeddings[,c(1,3)],dge@dr$tsne@cell.embeddings[,1:2],dge@dr$umap@cell.embeddings[,1:2])

### plot PCs and tSNE for each batch using the other batches as background
dge=SetAllIdent(dge,id="orig.ident")
sets=levels(dge@meta.data$orig.ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette[c(8:7,10,6:5,4:1)]
xlim=xlims[[1]]
ylim=ylims[[1]]

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


###### Rank correlation and Dissimilarity matrix for each normalized centroid using HVG
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

write.table(genecountsall,"Merged5Monkeysubjects-NoC1-NoC7-No10pctMT_13clusters_Centroid.txt",quote=F,sep="\t",row.names=T,col.names=T)
write.table(genecountsall,"Merged5Monkeysubjects-NoC1-NoC7-No10pctMT_12celltypes_Centroid.txt",quote=F,sep="\t",row.names=T,col.names=T)

cc=cor(as.matrix(genecountsall)[dge@var.genes,],method="spearman")
dim(cc)
min(cc) # 0.036

### labeling
data.use=cc[levels,levels]
colsep.use=cumsum(table(gsub("_.*","",levels)))
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels)))-table(gsub("_.*","",levels))/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)

### visualize for cell types
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep("white",length(sidecol[1,]))
sidecol[2,]=myBrewerPalette
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","CellType")
colnames(clab)=c("CellType","")

col.use=redblue100
pdf(file=paste0(dgefile,"5subjectsNoC1NoC5No10pctMT_12celltypes_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,7))
dev.off()
# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)-1):101]
pdf(file=paste0(dgefile,"5subjectsNoC1NoC5No10pctMT_12celltypes_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,7))
dev.off()


### Per-cell attributes
######### Analyze per-cell attributes 
###### calculate Gini index
library(ineq)
GiniAll=apply( dge@raw.data,2,function(x) ineq(x,type=c("Gini")) )
GiniNon0=apply( dge@raw.data,2,function(x) ineq(x[which(x!=0)],type=c("Gini")) )

dge=AddMetaData(dge,GiniAll,"GiniAll")
dge=AddMetaData(dge,GiniNon0,"GiniNon0")


pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.x", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.y", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.autosome", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="GiniAll", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="GiniNon0", nCol = 1,point.size.use =-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()


# 3.27.2019 Qianyi
# label cell types
### label major cell types
#1-2: Somatic - need to use somatic cells we used for somatic subclustering
#3-4: SPG 
#5-7: spermatocytes
#8-10: round spermatids
#11-13: elongating spermatids


somatic=1:2
spg=3:4
scyte=5:7
rs=8:10
es=11:13

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
2099           1155           3731           2945          11644 

table(dge@ident)
   1    2    3    4    5    6    7    8    9   10   11   12   13 
1667  432  464  691 2003 1067  661  573 1215 1157 3195 3758 4691

## label somatic cell types
somatic1=read.table("Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA_8celltypes_label.txt",row.names=1)
somatic=as.character(somatic1$V2)
names(somatic)=rownames(somatic1)

## There is 1 SPG cells labeled as somatic based on new global clustering
### I need to re-label it as SPG in order to keep consistency in cells used for somatic and SPG subclustering
tmp=names(id)[which(id=="Somatic")]
SPG1cell=tmp[which(!(tmp %in% names(somatic)))]
length(SPG1cell) # 1
SPG1cell # [1] "Monkey2_ATTGATCCTAGA"
id[SPG1cell] <- "Spermatogonia"
table(id)
Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
2098           1156           3731           2945          11644 


id=factor(id,levels=celltype2,order=T)
dge=AddMetaData(dge,id,"CellType2")
dge=SetAllIdent(dge,id="CellType2")
dge@ident=id
table(dge@ident)
Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
2098           1156           3731           2945          11644 

dge=SetAllIdent(dge,id="CellType2")
dge@ident=dge@meta.data$CellType2
names(dge@ident)=rownames(dge@meta.data)
germ=as.character(dge@ident)[which(dge@ident!="Somatic")]
names(germ)=names(dge@ident)[which(dge@ident!="Somatic")]
celltype=c(somatic,germ)
levels=c("Tcell","Macrophage","Endothelial","m-Pericyte","f-Pericyte","Myoid","ImmLeydig","DiffLeydig","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
celltype=factor(celltype,levels=levels,order=T)
table(celltype)
        Tcell     Macrophage    Endothelial     m-Pericyte     f-Pericyte 
            20             17            178            116             39 
         Myoid      ImmLeydig     DiffLeydig  Spermatogonia   Spermatocyte 
           362           1324             42           1156           3731 
RoundSpermatid     Elongating 
          2945          11644

dge=AddMetaData(dge,celltype,"CellType")
dge=SetAllIdent(dge,id="CellType")
dge@ident=celltype
table(dge@ident)

save(dge,file="MonkeyMerged5-NoC1-NoC7-No10pctMT.Robj")


### plot major cell types 
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


dge=SetAllIdent(dge,id="CellType2")
dge@ident=factor(dge@ident,levels=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"))

pdf("PCA_tSNE_UMAP_Merged5_5CellTypeswithMergedSomaticArm.pdf",width=11.3,height=8)
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

pdf("UMAP_col_Merged5_5CellTypeswithMergedSomaticArm.pdf",width=11.3,height=8)
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


### plot 12 major cell types
library(RColorBrewer)
col11 <- c(brewer.pal(12,"Paired")[c(12,11,10,6,8,7,9)],brewer.pal(8,"Set1")[8],brewer.pal(12,"Paired")[1:4])  #used this 1 
col22 <- c(brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3,8)],gg_color_hue(4))
col33 <- c(brewer.pal(8,"Set1")[c(7,6,4)],gg_color_hue(5),brewer.pal(12,"Paired")[1:4])
col44 <- c(brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3,8)],rainbow_hcl(4))
myBrewerPalette=col11

dge=SetAllIdent(dge,id="CellType")
dge@ident=dge@meta.data$CellType
names(dge@ident)=rownames(dge@meta.data)

pdf("PCA_tSNE_UMAP_Merged5_12CellTypes.pdf",width=11.5,height=8)
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

pdf("UMAP_col_Merged5_12CellTypes.pdf",width=11.3,height=8)
plot1=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col11)
plot2=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col22)
plot3=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col33)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=col44)
plot_grid(plot1,plot2,plot3,plot4,ncol = 2)
dev.off()


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



# 3.27.2019 Qianyi
### given there are many Ensembl IDs for monkey without annotated gene symbols
### convert Monkey Ensembl IDs without gene symbols to non-duplicated Human gene symbols with 1-1 orthologs
macGene=read.table("Monkey_EnsemblIDconvertedtoHumanGenes1to1Ortholog_Ensembl93.txt",stringsAsFactors=F)
macGene=unique(macGene)
dim(macGene) # [1] 19810     2

### check how many monkey Ensembl IDs do not have gene symbols
all=rownames(dge@data)
length(all) # 23029
length(grep("^ENSMMUG",all)) # 7333
miss=grep("^ENSMMUG",all,value=T)
keep=all[which(!grepl("^ENSMMUG",all))]
length(miss) # 7333
length(keep) # 15696

### Convert monkey Ensembl IDs without gene symbols to Human Gene symbols with 1-to-1 ortholog
length(which(miss %in% macGene[,1])) # 511
length(which(macGene[,1] %in% miss)) # 511
convertgene=macGene[which(macGene[,1] %in% miss),]
dim(convertgene) # [1] 511    2
anyDuplicated(convertgene[,1]) # 0

### remove the human gene symbol ortholog that are duplicated for already existing monkey Gene symbols
length(which(convertgene[,2] %in% keep)) # 43
convertgene=convertgene[which(!(convertgene[,2] %in% keep)),]
dim(convertgene) # [1] 468    2
anyDuplicated(convertgene[,1]) # 0

### Replace monkey Ensembl IDs without gene symbols by non-duplicated Human Gene symbols with 1-to-1 ortholog 
allc=all
names(allc)=all
for(i in 1:nrow(convertgene)){
allc[which(allc==convertgene[i,1])]<-convertgene[i,2]
}
length(allc) # [1] 23029
length(unique(allc)) # [1] 23029
anyDuplicated(allc) # 0
names(allc)=NULL
length(grep("^ENSMMUG",allc)) # 6865

hvg=dge@var.genes
length(grep("^ENSMMUG",hvg)) # 342
for(i in 1:nrow(convertgene)){
hvg[which(hvg==convertgene[i,1])]<-convertgene[i,2]
}
anyDuplicated(hvg) # 0
length(hvg) # 3284
length(grep("^ENSMMUG",hvg)) # 301

### save in Seurat object
rownames(dge@raw.data)=allc
rownames(dge@data)=allc
rownames(dge@scale.data)=allc
rownames(dge@dr$pca@gene.loadings.full)=allc
rownames(dge@hvg.info)=allc

dge@var.genes=hvg
rownames(dge@dr$pca@gene.loadings)=hvg

### save object for merged 5 Monkey subjects with Human gene symbol 1-to-1 ortholog
save(dge,file="MonkeyMerged5-NoC1-NoC7-No10pctMT_GeneNames.Robj")
dge # 23029 genes across 21574 samples



###### markers for each cell type against all others
dge=SetAllIdent(dge,id="CellType")
dge@ident=dge@meta.data$CellType
names(dge@ident)=rownames(dge@meta.data)
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"MonkeyMerged5subjects-NoC1-NoC7-No10pctMT_12celltypes_markersall_mindiff0.2_logfc2fold_3.27.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
         Tcell     Macrophage    Endothelial     m-Pericyte     f-Pericyte 
            20             17            178            116             39 
           260            380            297            344            304 
         Myoid      ImmLeydig     DiffLeydig  Spermatogonia   Spermatocyte 
           362           1324             42           1156           3731 
           347            299            404            846            822 
RoundSpermatid     Elongating 
          2945          11644 
           206            273 

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
write.table(centroid,"MonkeyMerged5-NoC1-NoC7-No10pctMT_12celltypes_centroids_allgenes.txt",row.names=T,col.names=T,quote=F,sep="\t")

### Genes Standardized Across Cell Types
# note: used this
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
write.table(centroid.std,"MonkeyMerged5-NoC1-NoC7-No10pctMT_12celltypes_centroid_allgenes_std.txt",row.names=T,col.names=T,quote=F,sep="\t")

centroid=read.table(paste0("MonkeyMerged5-NoC1-NoC7-No10pctMT_12celltypes_centroids_allgenes.txt"),header=T,row.names=1)
centroid.std=read.table(paste0("MonkeyMerged5-NoC1-NoC7-No10pctMT_12celltypes_centroid_allgenes_std.txt"),header=T,row.names=1)
markers=read.table("MonkeyMerged5subjects-NoC1-NoC7-No10pctMT_12celltypes_markersall_mindiff0.2_logfc2fold_3.27.2019.txt",header=T,row.names=1,stringsAsFactors=F)
table(dge@ident)
table(markers$cluster)

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
