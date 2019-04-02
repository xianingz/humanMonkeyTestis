### Global Clustering for Directly-merged all 5 monkey subjects after removing two doublets clusters
# by Qianyi on 3.19.2019
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
### directly-merge 5 monkey subjects after removing two doublets clusters (1/8 and 7/7)
#-> 13 clusters


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



# 3.18.2019 by Qianyi 

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

### load object for merged 5 Monkey subjects before removing doublets clusters
load(file="MonkeyMerged5.Robj")
dge5=dge

### remove two doublet clusters
table(dge@ident)
#          C1/8           C7/7        Somatic  Spermatogonia   Spermatocyte 
#           207             18           2337           1309           4015 
#RoundSpermatid     Elongating 
#          2751          11603 
cells.use=names(dge@ident)[which(!(dge@ident %in% c("C1/8","C7/7")))]
table(dge@ident[cells.use])
#          C1/8           C7/7        Somatic  Spermatogonia   Spermatocyte 
#             0              0           2337           1309           4015 
#RoundSpermatid     Elongating 
#          2751          11603 

dgedata2=dgedata[,cells.use]
dim(dgedata2)                   #[1] 23101 22015
nCellperGene <- rowSums(dgedata2>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use)              #[1] 23089
dgedata2=dgedata[genes.use,cells.use]
dim(dgedata2)                   #[1] 23089 22015

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project="MonkeyDirectlyMerged5-NoC1NoC7", min.cells=1, min.genes=1)
dge # 23089 genes across 22015 samples.

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
length(y) # 603
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
length(chrx.genes)  # 962
length(chry.genes)  # 48
length(mito.genes)  # 35
length(mito.genes2) # 28
dim(autosome.genes) # 22044  2
table(autosome.genes[,2])
1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
2265 1341 1361 1269 1003 1112 1521  859  903 1115 1263  763  814 1297  880 1201 
  17   18   19   20 
 465  357 1368  887

# %x expression
percent.x <- colSums(expm1(dge@data[chrx.genes, ]))/colSums(expm1(dge@data))
percent.y <- colSums(expm1(dge@data[chry.genes, ]))/colSums(expm1(dge@data))
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
percent.mito2 <- Matrix::colSums(dge@raw.data[mito.genes2, ])/Matrix::colSums(dge@raw.data)
percent.autosome <- colSums(expm1(dge@data[autosome.genes[,1], ]))/colSums(expm1(dge@data))
# named integer(0)

cells10pctMT=names(percent.mito)[which(percent.mito>=0.1)]
length(cells10pctMT)  # 441
write.table(cells10pctMT,"MonkeyCellsOver10pctMT.txt",col.names=F,row.names=F,quote=F,sep="\t")
# Identified cells with >=10% MT in the filtered cells; need to correct for the 37 MT gene list

tmp=dge5@meta.data[cells10pctMT,c("nGene","CellType2")]
germcells10pctMT=rownames(tmp)[which(tmp$CellType2 != "Somatic")]
length(germcells10pctMT) # 202
write.table(germcells10pctMT,"MonkeyGermCellsOver10pctMT.txt",col.names=F,row.names=F,quote=F,sep="\t")

pdf("percent.mito.pdf",height=4,width=4)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(percent.mito,percent.mito2,pch=16,cex=.5,col=rgb(0,0,0,0.2),xlab="%MT calculated using 35 genes",ylab="%MT calculated using 28 genes")
#points(percent.mito[which(percent.mito>=0.1)],percent.mito2[which(percent.mito>=0.1)],pch=16,cex=.5,col=rgb(1,0,0,0.3))
abline(v=0.1,col="red")
text(0.13,0.06,paste(length(which(percent.mito>=0.1)),"cells"),col="red")
text(0.13,0.055,"w/ >=10% MT",col="red")
text(0.06,0.01,paste(length(which(percent.mito<0.1)),"cells"))
text(0.06,0.005,"w/ <10% MT")
dev.off()
length(which(percent.mito>=0.1)) # 441
summary(dge@meta.data[cells10pctMT,]$percent.mito)
table(dge@meta.data$indiv)
table(dge@meta.data[cells10pctMT,]$indiv)
#Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
#   2085    2782    4966    6560    5622 
#     72      40     241      27      61 
table(dge5@ident)
table(dge5@ident[cells10pctMT])
#Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating
#2337           1309           4015          2751          11603 
#239            157             44            0              1
table(dge@ident)
table(dge@ident[cells10pctMT])
#  1   2   3   4   5   6   7   8   9  10  11  12  13 
#1878  459  550  764 1967 1064  745  567 1222 1166 2871 4411 4351
#209  30  85  72  11   0  33   0   0   0   0   0   1 
load(file="Monkey5Somatic-NoC1-NoC7_CCA.Robj")
table(testis@ident)
table(testis@ident[cells10pctMT])
#      Tcell  Macrophage Endothelial    Pericyte Subcluster2       Myoid 
#         20          17         170         131          43         516 
#          0           0           6          14           1          45 
#  ImmLeydig  DiffLeydig 
#       1396          44 
#        168           5
cells10pctMT=names(percent.mito)[which(percent.mito>=0.1)]
pdf("percent.mito_highlight.pdf",height=4,width=4)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
DimPlot(dge,reduction.use="umap",cex=.8,cols.use="grey",cells.highlight=cells10pctMT,cols.highlight = rgb(1,0,0,0.5))
DimPlot(testis,reduction.use="umap",cex=.8,cols.use="grey",cells.highlight=cells10pctMT,cols.highlight = rgb(1,0,0,0.5))
tmp=dge@dr$umap@cell.embeddings
plot(tmp,pch=16,cex=.5,col="grey")
points(tmp[cells10pctMT,],pch=16,cex=.5,col=rgb(1,0,0,0.5))
tmp=testis@dr$umap@cell.embeddings
plot(tmp,pch=16,cex=.5,col="grey")
points(tmp[cells10pctMT[cells10pctMT %in% rownames(tmp)],],pch=16,cex=.5,col=rgb(1,0,0,0.5))
dev.off()

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
print(length(dge@var.genes)) # 3412
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 23089 genes across 22015 samples.          
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
   2085    2782    4966    6560    5622
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
dgeall=dge
save(dgeall,file="MonkeyMerged5-NoC1-NoC7.Robj")
table(gsub("_.*","",names(dge@ident)))
Monkey1.1 Monkey1.2   Monkey2 Monkey3.1 Monkey3.2 Monkey4.1 Monkey4.2 Monkey5.1 
     1274       811      2782      2430      2536      2933      3627      2950 
Monkey5.2 
     2672

### PCA
print(Sys.time())   # [1] "2019-03-19 11:17:33 EDT" 
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
Sys.time() #[1] "2019-03-19 12:03:50 EDT"
dgeall=dge
save(dgeall,file="MonkeyMerged5-NoC1-NoC7.Robj")


### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.6)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
print(c( length(unique(dge@meta.data$res.1.3)),length(unique(dge@meta.data$res.1.4)),length(unique(dge@meta.data$res.1.5)),length(unique(dge@meta.data$res.1.6)),length(unique(dge@meta.data$res.1.7)),length(unique(dge@meta.data$res.1.8)),length(unique(dge@meta.data$res.1.9)),length(unique(dge@meta.data$res.2)),length(unique(dge@meta.data$res.2.1)),length(unique(dge@meta.data$res.2.2)),length(unique(dge@meta.data$res.2.3)),length(unique(dge@meta.data$res.2.4)),length(unique(dge@meta.data$res.2.5)),length(unique(dge@meta.data$res.2.6)),length(unique(dge@meta.data$res.2.7)),length(unique(dge@meta.data$res.2.8)),length(unique(dge@meta.data$res.2.9)),length(unique(dge@meta.data$res.3)) ))

table(dge@meta.data$res.0.5,dge@meta.data$res.0.6)
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
pdf(file=paste0("Merged5subjects-NoC1-NoC7_Centroid_norm_Seriation_",res,".pdf"))
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
save(dgeall,file="MonkeyMerged5-NoC1-NoC7.Robj")


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
     Somatic Spermatogonia Spermatocyte RoundSpermatid Elongating
  1     1878             0            0              0          0
  2      459             0            0              0          0
  3        0           549            1              0          0
  4        0           759            5              0          0
  5        0             0         1967              0          0
  6        0             0         1051             13          0
  7        0             1          729              9          6
  8        0             0            6            561          0
  9        0             0          219           1003          0
  10       0             0           23           1099         44
  11       0             0            7             66       2798
  12       0             0            7              0       4404
  13       0             0            0              0       4351
table(dge@ident,dge5@meta.data[rownames(dge@meta.data),]$res.0.5order)


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
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Merged5subjects-NoC1-NoC7_res.0.5order_markersall_mindiff0.2_logfc2fold_3.19.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)

   1    2    3    4    5    6    7    8    9   10   11   12   13 
1878  459  550  764 1967 1064  745  567 1222 1166 2871 4411 4351 
321 242 841 634 833 593  33 357 218  80  26   2  82 

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
pos=c("topright","bottomright","bottomright","topright")
# dge5 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-70,17),c(-70,17),c(-45,45),c(-8,14))
ylims[[1]]=list(c(-29,45),c(-46,24),c(-45,45),c(-15,10))

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

write.table(genecountsall,"Merged5subjects-NoC1-NoC7_13clusters_Centroid.txt",quote=F,sep="\t",row.names=T,col.names=T)

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
pdf(file=paste(dgename,"5subjects-NoC1-NoC7_13clusters_Centroid_RankedCorrelation_HVG.pdf",sep=""),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()


### Per-cell attributes
######### Analyze per-cell attributes 
###### calculate Gini index
library(ineq)
#for(i in 1:(length(dataset)+1)){
#  dge=dgelist[[i]]
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

