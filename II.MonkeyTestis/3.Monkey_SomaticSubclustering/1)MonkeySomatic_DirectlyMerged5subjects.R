### Somatic Subclustering for directly-merged 5 monkey subjects
# 1.4.2019 by Qianyi
### directly-merged 5 monkey subjects on 12/18/2018
#-> 12 clusters
#Cluster 1-2: somatic cells
#Cluster 3-4: SPG
#Cluster 5-8: Spermatocytes
#Cluster 9-10: Round spermatids
#Cluster 11-12: Elongating spermatids
### Extracted somatic cells and directly merged, did somatic subclustering

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


# 1.4.2019 by Qianyi 

# Subclustering for Somatic

### load object for merged 5 Monkey subjects
load(file="MonkeyMerged5.Robj")
dge5=dgeall
dge5 # 48655 genes across 22663 samples

### Extract Somatic cells (cluster 1-2 of res.0.5order)
somatic=1:2
cells.use=names(dge@ident)[which(dge@ident %in% somatic)]
table(dge@ident[cells.use])
table(dge@ident,dge@meta.data$indiv)
   1    2    3    4    5    6    7    8    9   10   11   12 
1901  661    0    0    0    0    0    0    0    0    0    0 
     Monkey1 Monkey2 Monkey3 Monkey4 Monkey5
  1      164     146    1350     119     122
  2      245      36     178     154      48
  3       73      79     193     161      46
  4      578      55      46      30      48
  5      232     267     365      47     620
  6       35     148     298      30     958
  7       81     127     174      98     535
  8       71     155     175     316     541
  9       52     294     128     541     478
  10     121     395     323     857     735
  11     274     633     835    1850    1001
  12     182     469     993    2422     513

dgedata=dge@raw.data[,cells.use]
dim(dgedata) # [1] 23101  2562
nCellperGene <- rowSums(dgedata>0)
length(which(nCellperGene==0))     #4731
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgedata2=dgedata[genes.use,]
dim(dgedata2) # [1] 18370  2562

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project="MonkeyMerged5Somatic", min.cells=1, min.genes=1)
MT <- c("ENSMMUG00000028704", "ENSMMUG00000028703", "ENSMMUG00000028702", "ENSMMUG00000028701", "ENSMMUG00000028700", "ENSMMUG00000028699", "ENSMMUG00000028698", "ENSMMUG00000028697", "ENSMMUG00000028696", "ENSMMUG00000028694", "ENSMMUG00000028693", "ENSMMUG00000028692", "ENSMMUG00000028691", "ENSMMUG00000028690", "ENSMMUG00000028688", "ENSMMUG00000028687", "ENSMMUG00000028686", "ENSMMUG00000028685", "ENSMMUG00000028684", "ENSMMUG00000028683", "ENSMMUG00000028681", "ENSMMUG00000028679", "ENSMMUG00000028678", "ENSMMUG00000028676", "ENSMMUG00000028675", "ENSMMUG00000028674", "ENSMMUG00000028671", "ENSMMUG00000028670", "ENSMMUG00000028669", "ENSMMUG00000028668", "ENSMMUG00000028695", "MT-CO1", "MT-CO3", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6")
MT=unique(MT)
mito.genes=MT[which(MT %in% rownames(dge@data))]
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("Monkey_Merged5Somatic_HVG0.2.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # 2330
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 18370 genes across 2562 samples      
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
    409     182    1528     273     170
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
dge5=dge
dge5=dge
save(dge,file="MonkeyMerged5Somatic.Robj")
table(dge@meta.data$indiv)
Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
    409     182    1528     273     170
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

### PCA Plot
pdf("PCA_Merged5Somatic.pdf",width=12,height=10)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
### Scree Plot for PCA
numPCs=6;i=1 # HVG
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
Sys.time() 
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,3,by=0.1), save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])
Sys.time() 
dge5=dge
save(dge,file="MonkeyMerged5Somatic.Robj")


# 1.7.2019
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
pdf(file=paste0("Merged5SomaticSubjects_Centroid_norm_Seriation_",res,".pdf"))
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


## double-check if I need to reverse the cluster ID orders
## after reverse cluster order, repeat the above plots
### compare with previous somatic clusters of merged 5 subjects 
dge5=dge
load(file="MonkeyMerged5.Robj")
dgeall 
dge=dgeall
table(dge5@ident,dgeall@ident[names(dge5@ident)])

### visualize known markers
knownmarkers=c("VIM",
  "CD163","S100A4","TYROBP","LYZ","RGS1",
  "VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN",
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","PDGFRA","PDGFRB","PDGFB","ALDH1",
  "DLK1","HSD17B3","STAR","CYP17A1","NR5A1","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "CLU","SOX9","AMH","AMHR2"  )
length(knownmarkers) # 36
knownmarkers[which(!(knownmarkers %in% rownames(dge@data)))] #[1] "PDGFRA" "ALDH1"  "INHBA"  "AMH"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(dge@data))]
length(knownmarkers) # 32

pdf("knownmarkers_Feature.pdf",height=15,width=21)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf("knownmarkers_Violin.pdf",height=10,width=22)
VlnPlot(dge,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()


# Convert Monkey Ensembl IDs without annotated gene symbols to non-duplicated human gene symbols with 1-1 orthlogues
### load monkey Ensembl IDs converted to monkey/human gene symbols
macGene=read.table("Monkey_EnsemblIDconvertedtoHumanGenes1to1Ortholog_Ensembl93.txt",stringsAsFactors=F)
macGene=unique(macGene)
dim(macGene) # [1] 19810     2
### load object for Somatic cells of merged 5 Monkey subjects
load(file="MonkeyMerged5Somatic.Robj")
dge5=dge # 18370 genes across 2562 samples.

### check how many monkey Ensembl IDs do not have gene symbols
all=rownames(dge@data)
length(all) # 18370
length(grep("^ENSMMUG",all)) # 4523
miss=grep("^ENSMMUG",all,value=T)
keep=all[which(!grepl("^ENSMMUG",all))]
length(miss) # 4523
length(keep) # 13847

### Convert monkey Ensembl IDs without gene symbols to Human Gene symbols with 1-to-1 ortholog
length(which(miss %in% macGene[,1])) # 403
length(which(macGene[,1] %in% miss)) # 403
convertgene=macGene[which(macGene[,1] %in% miss),]
dim(convertgene) # [1] 403    2
anyDuplicated(convertgene[,1]) # 0

### remove the human gene symbol ortholog that are duplicated for already existing monkey Gene symbols
length(which(convertgene[,2] %in% keep)) # 35
convertgene=convertgene[which(!(convertgene[,2] %in% keep)),]
dim(convertgene) # [1] 368    2
anyDuplicated(convertgene[,1]) # 0

### Replace monkey Ensembl IDs without gene symbols by non-duplicated Human Gene symbols with 1-to-1 ortholog 
allc=all
names(allc)=all
for(i in 1:nrow(convertgene)){
allc[which(allc==convertgene[i,1])]<-convertgene[i,2]
}
length(allc)         # [1] 18370
length(unique(allc)) # [1] 18370
anyDuplicated(allc)  # 0
names(allc)=NULL
length(grep("^ENSMMUG",allc)) # 4155

hvg=dge@var.genes
length(grep("^ENSMMUG",hvg)) # 185
for(i in 1:nrow(convertgene)){
hvg[which(hvg==convertgene[i,1])]<-convertgene[i,2]
}
anyDuplicated(hvg) # 0
length(hvg)
length(grep("^ENSMMUG",hvg)) # 140


### save in Seurat object
rownames(dge@raw.data)=allc
rownames(dge@data)=allc
rownames(dge@scale.data)=allc
rownames(dge@dr$pca@gene.loadings.full)=allc
rownames(dge@hvg.info)=allc

dge@var.genes=hvg
rownames(dge@dr$pca@gene.loadings)=hvg


### save object for merged 5 Monkey subjects with Human gene symbol 1-to-1 ortholog
save(dge,file="MonkeyMerged5Somatic_GeneNames.Robj")
dge # 48655 genes across 22663 samples

### save monkey Ensembl IDs converted to Human gene symbol 1-to-1 ortholog in my data
write.table(convertgene,"Monkey_EnsemblIDconvertedtoHumanGenes1to1Ortholog_Ensembl93_MonkeySomatic.txt",col.names=F,row.names=F,quote=F,sep="\t")


# 1.26.2019 Visualize IntProg/Myoid/Leydig markers in Object with Ensembl IDs converted to known human gene symbols with one-to-one orthologue
load(file="MonkeyMerged5Somatic_GeneNames.Robj")
dge

knownmarkers=c(
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","TCF21","PDGFRA","PDGFRB","PDGFB",
  "DLK1","HSD17B3","STAR","CYP17A1","GATA4","NR5A1","EGR2","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "THY1","ENG","VCAM1","ALCAM","CD44","NT5E","PDGFA","CASP3","CLU" )
length(knownmarkers) # 31
knownmarkers[which(!(knownmarkers %in% rownames(dge@data)))] 
#[1] "PDGFRA" "CASP3"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(dge@data))]
length(knownmarkers) # 29

pdf("knownmarkers_IntProg_Feature.pdf",height=15,width=21)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf("knownmarkers_IntProg_Violin.pdf",height=10,width=22)
VlnPlot(dge,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
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
knownmarkers[which(!(knownmarkers %in% rownames(dge@data)))] 
#[1] "HAND2"      "BHLHA9"     "PTF1A"      "HAND1"      "TCF7L1-IT1"
#[6] "FIGLA"      "FERD3L" 
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(dge@data))]
length(knownmarkers) # 17

pdf("knownmarkers_TCF21paralogs_Feature.pdf",height=9,width=21)
FeaturePlot(object = dge, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
dev.off()
pdf("knownmarkers_TCF21paralogs_Violin.pdf",height=6,width=22)
VlnPlot(dge,knownmarkers,cols.use=myBrewerPalette,nCol=7,point.size.use=-1)
dev.off()



###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Merged5Somaticsubjects_res.0.5order_markersall_mindiff0.2_logfc2fold_1.7.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)

  1   2   3   4   5   6   7   8 
207 534 755 604 186  48 191  37 
585   1  11  26 113  92 125 110

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
pdf("PCA_tSNE_UMAP_Merged5SomaticMonkey_Rep.pdf",width=10.5,height=8)
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
pdf("PCA_tSNE_UMAP_Merged5SomaticMonkey_Subject.pdf",width=10,height=8)
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
pdf("PCA_tSNE_UMAP_Merged5SomaticMonkey_res.0.5order.pdf",width=9,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
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
pos=c("bottomleft","bottomleft","topleft","bottomleft")
# dge5 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-39,8),c(-39,8),c(-9,10),c(-45,32))
ylims[[1]]=list(c(-35,8),c(-26,14),c(-8,17),c(-39,30))

dim=list(dge@dr$pca@cell.embeddings[,1:2],dge@dr$pca@cell.embeddings[,c(1,3)],dge@dr$umap@cell.embeddings[,1:2],dge@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[1]]
ylim=ylims[[1]]

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
  write.table(ncellscluster,"ncellspercluster_batch_top6PCs.txt",quote=F,row.names=T,col.names=T,sep="\t")


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
jpeg(file=paste0("OrigSet_8clusters_",label,"_bg.jpeg"),res=300,height=1200,width=3500)
multiplot(plotset,cols = 5)
dev.off()
pdf(file=paste0("OrigSet_8clusters_",label,"_bg.pdf"),height=6,width=15)
multiplot(plotset,cols = 5)
dev.off()




# 1.22.2019 Qianyi
###### label cell types and remove myoid-spermatid doublet cluster
dge=SetAllIdent(dge,id="res.0.5order")
table(dge@ident)
  1   2   3   4   5   6   7   8 
207 534 755 604 186  48 191  37

### remove somatic-Scyte/Spermatid doublet Cluster 1
dge=SubsetData(dge,ident.remove=1)
table(dge@ident)
  2   3   4   5   6   7   8 
534 755 604 186  48 191  37 

### label cell types
id=as.numeric(dge@ident)
table(id)
names(id)=names(dge@ident)
celltype=c("IntProgenitor","IntProgenitor","Myoid","Pericyte","Leydig","Endothelial","Macrophage")
for(i in 1:7){
id[which(id==i)]<-celltype[i]
}
celltype=rev(unique(celltype))
id=factor(id,levels=celltype,order=T)
table(id)
   Macrophage   Endothelial        Leydig      Pericyte 
           37           191            48           186 
        Myoid IntProgenitor 
          604          1289 

table(dge@ident)
#celltype=rev(unique(celltype))
#id=factor(id,levels=celltype,order=T)
dge=AddMetaData(dge,id,"CellType")
dge=SetAllIdent(dge,id="CellType")
dge@ident=id
table(dge@ident)
   Macrophage   Endothelial        Leydig      Pericyte         Myoid 
           37           191            48           186           604 
IntProgenitor 
         1289 

save(dge,file="MonkeyMerged5Somatic-NoC1.Robj")
write.table(dge@ident,"MonkeyMerged5Somatic-NoC1_6celltypes_label.txt",row.names=T,col.names=F,quote=F,sep="\t")

# try to match color scheme for corresponding mouse cell types
#levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown")
#myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7)] # 1

celltype=c("Macrophage","Endothelial","Leydig","Pericyte","Myoid","IntProgenitor")
library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(11,10,9,8,6,7)] 
# used this color scheme

pdf("PCA_tSNE_UMAP_Merged5SomaticMonkey_NoC1_6celltypes.pdf",width=10.5,height=8)
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


###### Rank correlation for all cells
nrho=cor(as.matrix(dge@data)[dge@var.genes,],method="spearman")
length(dge@var.genes) # 3518
testcor=nrho

###### Jaccard distance
testcor=as.matrix(dge@snn)


### order cells by clusters 
res="res.0.5order";j=1;resi=1;
   dge=SetAllIdent(dge,id=res)
   print(length(unique(dge@ident)))
   TSNEPlot(dge)

## order cell by cluster ID and randomly shuffle cells within each batch
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

jpeg(file=paste(dgename,"5subjectsSomatic_8clusters_RankCor_HVG_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgename,"5subjectsSomatic_8clusters_RankCor_HVG.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()

col.use2=redblue100[c(rep(c(41:100),each=10),rep(100,100000))]
length(col.use2)
jpeg(file=paste(dgename,"5subjectsSomatic_8clusters_snn_6k.jpeg",sep=""),height=6000,width=6000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=3,cexRow=3,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
jpeg(file=paste(dgename,"5subjectsSomatic_8clusters_snn.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=1.5,cexRow=1.5,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()


###### Rank correlation and Dissimilarity matrix for each normalized centroid using HVG
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

write.table(genecountsall,"Merged5Somaticsubjects_8clusters_Centroid.txt",quote=F,sep="\t",row.names=T,col.names=T)

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

pdf(file=paste(dgename,"5subjectsSomatic_8clusters_Centroid_RankedCorrelation_HVG.pdf",sep=""),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(5,3))
dev.off()

# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)):101]
pdf(file=paste(dgename,"5subjectsSomatic_8clusters_Centroid_RankedCorrelation_HVG_2.pdf",sep=""),height=6.2,width=6)
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
pdf(file=paste0(dgefile,"5subjectsSomaticNoC1_6celltypes_Centroid_RankedCorrelation_HVG.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()
# modify color scheme so that it is roughly blue to red from 0 to 1
col.use=redblue100[max(1,round(min(c(cc))*100)-1):101]
pdf(file=paste0(dgefile,"5subjectsSomaticNoC1_6celltypes_Centroid_RankedCorrelation_HVG_2.pdf"),height=6.2,width=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,5.5))
dev.off()


### Per-cell attributes
pdf(file=paste0(dgefile,"PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="percent.mito", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()

