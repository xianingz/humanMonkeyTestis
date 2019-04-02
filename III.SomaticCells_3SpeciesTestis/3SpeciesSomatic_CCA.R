# Subclustering for Somatic cells from all 3 species merged by CCA
# 2.7.2019 by Qianyi
### Compare somatic of 3 species
### treated mouse as 2 subjects (old 24batches vs INT6 NewSca1 performed by two different hands)

home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/"
setwd(home)
dgefile=dgename="3SpeciesComparison_plot/CCA_"
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(cowplot)
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2")[c(4,8,1)],brewer.pal(8,"Set2")[c(4,8,1)])
redblue100<-rgb(read.table('redblue100.txt',sep='\t',row.names=1,header=T))
gg_color_hue <- function(n) {
hues = seq(15, 375, length = n + 1)
hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette1=gg_color_hue(3)
source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/Rcode_multiplot.R")


species=c("Human","Monkey","Mouse")
speciestags=c("H","M","m")

## load Seurat object for 3 species
load(file=paste0("human/Human4-1235-NoC5SomaticNoC1_CCA.Robj"))
dgeH=pbmc
load(file=paste0("monkey/Monkey5Somatic-NoC1_CCA.Robj"))
dgeM=pbmc
load(file=paste0("MouseAdultST25Somatic.Robj"))
dgem=dge

dgelist=list(dgeH,dgeM,dgem)

## Each batch for each species
batches=list()
for(i in 1:2){
  batches[[i]]=unique(dgelist[[i]]@meta.data$orig.ident)
}
i=3;
  batches[[i]]=levels(dgelist[[i]]@data.info$orig.ident)
for(i in 1:3){
  print(length(batches[[i]]))
}
#[1] 10
#[1] 9
#[1] 25
batches
#[[1]]
# [1] "Human1.1" "Human1.2" "Human1.3" "Human1.4" "Human1.5" "Human2.1"
# [7] "Human2.2" "Human3.1" "Human3.2" "Human5"  
#[[2]]
#[1] "Monkey1.1" "Monkey1.2" "Monkey2"   "Monkey3.1" "Monkey3.2" "Monkey4.1"
#[7] "Monkey4.2" "Monkey5.2" "Monkey5.1"
#[[3]]
# [1] "ST1"  "ST2"  "ST3"  "ST4"  "ST5"  "ST6"  "ST7"  "ST8"  "SPG1" "SPG2"
#[11] "SPG3" "INT1" "INT2" "INT3" "INT4" "INT5" "INT6" "SER1" "SER2" "SER3"
#[21] "SER4" "SER5" "SER6" "SER7" "SER8"

## Each subject for each species
subjects=list()
for(i in 1:2){
  subjects[[i]]=unique(dgelist[[i]]@meta.data$protocol)
}
i=3;
  subjects[[i]]=c("Mouse24","MouseINT6")
for(i in 1:3){
  print(length(subjects[[i]]))
}
#[1] 4
#[1] 5
#[1] 2
subjects
#[[1]]
#[1] "Human1" "Human2" "Human3" "Human5"
#[[2]]
#[1] "Monkey1" "Monkey2" "Monkey3" "Monkey4" "Monkey5"
#[[3]]
# [1] "Mouse24" "MouseINT6"


## check 1-1-1 gene orthologues present in 3-species data
Orth1to1to1=read.table("humMacMouOrth.1to1to1.txt",stringsAsFactors=F)

allgenes=Orth1to1to1
for(i in 1:length(dgelist)){
  genes=rownames(dgelist[[i]]@data)
  allgenes=allgenes[which(allgenes[,i] %in% genes),]
}
dim(allgenes) # [1] 10765     3
for(i in 1:length(dgelist)){
  print(length(unique(allgenes[,i])))
}
#[1] 10765
#[1] 10765
#[1] 10765

### keep 1-1-1 ortholgue genes for 3 species data
geneslist=list()
for(i in 1:length(dgelist)){
  genes=rownames(dgelist[[i]]@data)
  genes=genes[which(genes %in% allgenes[,i])]
  geneslist[[i]]=genes
  print(length(genes))
}
#[1] 10765
#[1] 10765
#[1] 10765



## Extract Somatic cells from each subject of each species
datalist=list()
for(i in 1:length(dgelist)){
  datalist[[i]]=list()
  subject=subjects[[i]]
  dgeall=dgelist[[i]]
  genes.use=geneslist[[i]]
  for(j in 1:length(subject)){
### Extract somatic cells from each subject of each species
    if(i!=3){
      cells.use=names(dgeall@ident)[which(dgeall@meta.data$protocol == subject[j])]
    } else if(j==1) {
      cells.use=names(dgeall@ident)[which(dgeall@data.info$orig.ident != "INT6")]
    } else {
      cells.use=names(dgeall@ident)[which(dgeall@data.info$orig.ident == "INT6")]
    }
    dgedata=dgeall@raw.data[genes.use,cells.use]
    print(dim(dgedata))
### replace gene names in monkey and mouse by human gene names for data
    for(g in 1:nrow(dgedata)){
      rownames(dgedata)[g]=allgenes[which(allgenes[,i]==rownames(dgedata)[g]),1]
    }
    dge <- CreateSeuratObject(raw.data = dgedata)
    dge <- NormalizeData(object = dge)
    dge <- ScaleData(object = dge)
    dge <- FindVariableGenes(object = dge, do.plot = FALSE)
    dge@meta.data[,"protocol"] <- subject[j]
    dge@meta.data[,"species"] <- species[i]
    datalist[[i]][[j]]=dge
  }
}

datalist=unlist(datalist,recursive=FALSE)
subject=unlist(subjects)

# we will take the union of the top 2k variable genes in each dataset for
# alignment, note that we use 2k genes in the manuscript examples, you can
# try different number of top HVG here with negligible changes to the overall results
hvglist=list()
for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(dge)
  hvglist[[i]] <- rownames(x = head(x = dge@hvg.info, n = 2000))
  print(length(hvglist[[i]])) # 2000
}
hvg.union=unique(unlist(hvglist))
length(hvg.union) # 8577

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
length(hvg.union) 

for(i in 1:length(subject)){
  dge=datalist[[i]]
  print(length(which(hvg.union %in% rownames(x = head(x = dge@hvg.info, n = 2000)))))
}


# check the number of hvg overlapped between datasets
cc=matrix(0,length(subject),length(subject))
for(i in 1:length(subject)){
    for(j in 1:length(subject)){
        cc[i,j]=length(intersect(rownames(x = head(x = datalist[[i]]@hvg.info, n = 2000)),rownames(x = head(x = datalist[[j]]@hvg.info, n = 2000))))
    }
}
cc
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
 [1,] 2000  637  658  594  643  548  584  585  583   571   569
 [2,]  637 2000  712  705  635  515  651  580  573   624   705
 [3,]  658  712 2000  698  597  506  645  572  572   625   681
 [4,]  594  705  698 2000  617  481  674  579  568   603   701
 [5,]  643  635  597  617 2000  628  741  672  720   605   586
 [6,]  548  515  506  481  628 2000  595  626  581   502   494
 [7,]  584  651  645  674  741  595 2000  664  673   686   717
 [8,]  585  580  572  579  672  626  664 2000  659   597   601
 [9,]  583  573  572  568  720  581  673  659 2000   600   582
[10,]  571  624  625  603  605  502  686  597  600  2000   864
[11,]  569  705  681  701  586  494  717  601  582   864  2000


### run multi CCA
testis <- RunMultiCCA(object.list = datalist, genes.use = hvg.union, num.ccs=30)
save(testis,file=paste0("Somatic_3Species_11Subjects_CCA.Robj"))
print(length(testis@var.genes)) # [1] 5186

dgefile=dgename="3SpeciesComparison_plot/CCA_"

p1 <- DimPlot(object = testis, reduction.use = "cca",cols.use=myBrewerPalette, group.by = "protocol", pt.size = 1.5, 
    do.return = TRUE)
p2 <- VlnPlot(object = testis, features.plot = "CC1",cols.use=myBrewerPalette, group.by = "protocol", do.return = TRUE)
pdf(paste(dgefile,"dge_CCA_orig.pdf",sep=""),height=4.5,width=11)
plot_grid(p1, p2)
dev.off()

pdf(paste(dgefile,"dge_CCs_top20CCs.pdf",sep=""),height=5)
MetageneBicorPlot(testis, grouping.var = "protocol", dims.eval = 1:20, 
    display.progress = FALSE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 1:6, 
    do.balanced = TRUE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 7:12, 
    do.balanced = TRUE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 13:18, 
    do.balanced = TRUE)
dev.off()

pdf(paste(dgefile,"dge_CCs_top30CCs.pdf",sep=""),height=5)
MetageneBicorPlot(testis, grouping.var = "protocol", dims.eval = 1:30, 
    display.progress = FALSE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 1:6, 
    do.balanced = TRUE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 7:12, 
    do.balanced = TRUE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 13:18, 
    do.balanced = TRUE)
DimHeatmap(object = testis, reduction.type = "cca", cells.use = 500, dim.use = 19:24, 
    do.balanced = TRUE)
dev.off()

#Rescaling group 1
#Rescaling group 2
#Rescaling group 3
#Rescaling group 4
#Rescaling group 5
#Rescaling group 6
#Rescaling group 7
#Rescaling group 8
#Rescaling group 9
#Rescaling group 10
#Rescaling group 11
#`geom_smooth()` using method = 'loess' and formula 'y ~ x'
#`geom_smooth()` using method = 'loess' and formula 'y ~ x'


## Add cell types from individual species analysis
ident=levels=NULL
for(i in 1:length(speciestags)){
  tmp=paste(speciestags[i],as.character(dgelist[[i]]@ident),sep="_")
  names(tmp)=names(dgelist[[i]]@ident)
  ident=c(ident,tmp)
  tmp=paste(speciestags[i],levels(dgelist[[i]]@ident),sep="_")
  if(i==1){tmp=rev(tmp)}
  if(i==3){
    tmp=paste(speciestags[i],levels(dgelist[[i]]@data.info$CellType),sep="_")
  }
  levels=c(levels,tmp)
}  
ident=factor(ident,levels=levels)
table(ident)
             H_7              H_6              H_5              H_4 
              65              231              297              519 
             H_3              H_2              H_1              M_1 
            2121              410               79               37 
             M_2              M_3              M_4              M_5 
             175               42              180              521 
             M_6              M_7 m_InnateLymphoid     m_Macrophage 
            1382               18               64              139 
   m_Endothelial          m_Myoid        m_Unknown         m_Leydig 
             179               49             2205              314 
       m_Sertoli 
            2131 

testis=AddMetaData(testis,ident,"SingleSpeciesClusters")
testis=SetAllIdent(testis,"SingleSpeciesClusters")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$SingleSpeciesClusters))
table(testis@ident)
save(testis,file=paste0("Somatic_3Species_11Subjects_CCA.Robj"))


### Before we align the subspaces, we first search for cells whose expression profile cannot be well-explained by low-dimensional CCA, compared to low-dimensional PCA.
### because we would like to know if all cells match between two batches
numCCs=17

testis <- CalcVarExpRatio(object = testis, reduction.type = "pca", grouping.var = "protocol", 
    dims.use = 1:numCCs)
summary(testis@meta.data$var.ratio.pca)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03826 0.55232 0.73917 0.71282 0.87395 2.87517 

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
 
testis.all.save <- testis
testis <- SubsetData(object = testis, subset.name = "var.ratio.pca", accept.low = 0.5)
# For illustrative purposes, we can also look at the discarded cells. Note that discarded cells tend to have lower gene counts, but a subset also express high levels of PF4. This is because megakaryocytes are only present in the 10X dataset, and thus are correctly discarded as ‘dataset-specific’]
testis.discard <- SubsetData(object = testis.all.save, subset.name = "var.ratio.pca", 
    accept.high = 0.5)

table(testis.all.save@meta.data$protocol)
table(testis.discard@meta.data$protocol)
table(testis@meta.data$protocol)

table(testis.all.save@meta.data$SingleSpeciesClusters)
table(testis.discard@meta.data$SingleSpeciesClusters)
table(testis@meta.data$SingleSpeciesClusters)

> table(testis.all.save@meta.data$SingleSpeciesClusters)
            H_7              H_6              H_5              H_4 
              65              231              297              519 
             H_3              H_2              H_1              M_1 
            2121              410               79               37 
             M_2              M_3              M_4              M_5 
             175               42              180              521 
             M_6              M_7    m_Endothelial m_InnateLymphoid 
            1382               18              179               64 
        m_Leydig     m_Macrophage          m_Myoid        m_Unknown 
             314              139               49             2205 
       m_Sertoli 
            2131 
> table(testis.discard@meta.data$SingleSpeciesClusters)
             H_7              H_6              H_5              H_4 
               3               18               72              221 
             H_3              H_2              H_1              M_1 
             530               83               28               14 
             M_2              M_3              M_4              M_5 
              98               34               68              227 
             M_6              M_7    m_Endothelial m_InnateLymphoid 
             648               15                4                1 
        m_Leydig     m_Macrophage          m_Myoid        m_Unknown 
               1                0                9               99 
       m_Sertoli 
               0 
> table(testis@meta.data$SingleSpeciesClusters)
             H_7              H_6              H_5              H_4 
              62              213              225              298 
             H_3              H_2              H_1              M_1 
            1591              327               51               23 
             M_2              M_3              M_4              M_5 
              77                8              112              294 
             M_6              M_7    m_Endothelial m_InnateLymphoid 
             734                3              175               63 
        m_Leydig     m_Macrophage          m_Myoid        m_Unknown 
             313              139               40             2106 
       m_Sertoli 
            2131 



testis <- testis.all.save
# do not discard cells



###### Now we align the CCA subspaces, which returns a new dimensional reduction called cca.aligned
numCCs=17 # used this after rescaling group

testis <- AlignSubspace(object = testis, reduction.type = "cca", grouping.var = "protocol", 
    dims.align = 1:numCCs)

### Visualize the aligned CCA and perform integrated analysis
p1 <- VlnPlot(object = testis, features.plot = "CC1", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p2 <- VlnPlot(object = testis, features.plot = "CC2", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p3 <- VlnPlot(object = testis, features.plot = "CC3", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p4 <- VlnPlot(object = testis, features.plot = "CC4", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p11 <- VlnPlot(object = testis, features.plot = "ACC1", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p22 <- VlnPlot(object = testis, features.plot = "ACC2", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p33 <- VlnPlot(object = testis, features.plot = "ACC3", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
p44 <- VlnPlot(object = testis, features.plot = "ACC4", group.by = "species", 
    do.return = TRUE,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)
pdf(paste(dgefile,numCCs,"CCs_dge_ACCs_top.pdf",sep=""),height=6,width=15)
plot_grid(p1, p2,p3,p4,p11,p22,p33,p44,ncol=4)
dev.off()

###### Now we can run a single integrated analysis on all cells!

testis <- RunTSNE(object = testis, reduction.use = "cca.aligned", dims.use = 1:numCCs, 
    do.fast = TRUE)
testis <- RunUMAP(testis, reduction.use = "cca.aligned", dims.use = 1:numCCs)

###### Visualize the original batch identity in ACC and tSNE spaces
plotlist=list()
plotlist[[1]]=DimPlot(testis,group.by = "protocol",reduction.use = "cca.aligned",1,2,do.return=TRUE)
plotlist[[2]]=DimPlot(testis,group.by = "protocol",reduction.use = "cca.aligned",1,3,do.return=TRUE)
plotlist[[3]]=DimPlot(testis,group.by = "protocol",reduction.use = "umap",do.return=TRUE)
plotlist[[4]]=TSNEPlot(testis,group.by = "protocol",do.return=TRUE)
pdf(paste(dgefile,"ACCs_tSNE_orig.pdf",sep=""),height=8,width=10.5)
multiplot(plotlist,cols = 2)
dev.off()

###### Louvain-jaccard clustering
testis <- FindClusters(object = testis, reduction.type = "cca.aligned", dims.use = 1:numCCs, 
    resolution=seq(0.1,3,by=0.1),save.SNN = TRUE)
save(testis,file=paste0("Somatic_3Species_11Subjects_CCA_",numCCs,"CCs.Robj"))
testis3=testis

print(c( length(unique(testis@meta.data$res.0.1)),length(unique(testis@meta.data$res.0.2)),length(unique(testis@meta.data$res.0.3)),length(unique(testis@meta.data$res.0.4)),length(unique(testis@meta.data$res.0.5)),length(unique(testis@meta.data$res.0.6)),length(unique(testis@meta.data$res.0.7)),length(unique(testis@meta.data$res.0.8)),length(unique(testis@meta.data$res.0.9)),length(unique(testis@meta.data$res.1)) ))
print(c( length(unique(testis@meta.data$res.1.1)),length(unique(testis@meta.data$res.1.2)),length(unique(testis@meta.data$res.1.3)),length(unique(testis@meta.data$res.1.4)),length(unique(testis@meta.data$res.1.5)),length(unique(testis@meta.data$res.1.6)),length(unique(testis@meta.data$res.1.7)),length(unique(testis@meta.data$res.1.8)),length(unique(testis@meta.data$res.1.9)),length(unique(testis@meta.data$res.2)) ))


ACCPlot <- function(object, ...) {
  return(DimPlot(object = object, reduction.use = "cca.aligned", label.size = 4, ...))
}

res=paste0("res.",c(paste0("0.",1:9),1,paste0("1.",1:9))) # res.0.8 for 10 clusters
plotlist=list()
pdf(paste(dgefile,numCCs,"CCs_clusters.pdf",sep=""),height=7.5,width=9)
for(resi in 1:length(res)){
testis=SetAllIdent(testis,id=res[resi])
plotlist[[1]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(testis,reduction.use = "umap",cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[4]]=TSNEPlot(testis,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)
}
dev.off()



table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.3)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.4)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.5)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.6)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.7)
table(testis@meta.data$res.0.7,testis@meta.data$res.0.6)
table(testis@meta.data$res.0.7,testis@meta.data$res.0.8)
table(testis@meta.data$res.0.7,testis@meta.data$res.0.9)
table(testis@meta.data$res.0.7,testis@meta.data$res.1)
table(testis@meta.data$res.0.7,testis@meta.data$res.1.1)


  load(file=paste0("Somatic_3Species_11Subjects_CCA.Robj"))
  testis3=testis

###### order clusters for each dataset
testis=testis3;
numCCs=17;res=paste0("res.0.",3:7);
for(resi in 1:length(res)){


testis=SetAllIdent(testis,id=res[resi])
#testis <- BuildClusterTree(testis, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs)
#testis@meta.data[,res[resi]]=testis@ident
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

library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

pdf(file=paste(dgename,"Centroid_RankedCorrelation_",numCCs,"CCs_",res[resi],".pdf",sep=""),height=5.5,width=5)
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

ordered=paste0(res[resi],"order")

testis=AddMetaData(testis,cells.ident.ordered,ordered)
testis@meta.data[,ordered]=factor(testis@meta.data[,ordered])
testis=SetAllIdent(testis,ordered)


}
# save the testis file
save(testis,file=paste0("Somatic_3Species_11Subjects_CCA_",numCCs,"CCs.Robj"))
save(testis,file=paste0("Somatic_3Species_11Subjects_CCA.Robj"))
testis3=testis

# go above to re-plot PCA and tSNE

res=paste0("res.0.",3:7,"order");
plotlist=list()
pdf(paste(dgefile,numCCs,"CCs_clusters_order.pdf",sep=""),height=7.5,width=9)
for(resi in 1:length(res)){
testis=SetAllIdent(testis,id=res[resi])
plotlist[[1]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,2,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[3]]=DimPlot(testis,reduction.use = "cca.aligned",cols.use=myBrewerPalette,1,3,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[2]]=DimPlot(testis,reduction.use = "umap",cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
plotlist[[4]]=TSNEPlot(testis,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
multiplot(plotlist,cols = 2)
}
dev.off()


table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.3order)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.4order)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.5order)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.6order)
table(testis@meta.data$SingleSpeciesClusters,testis@meta.data$res.0.7order)
table(testis@meta.data$res.0.7order,testis@meta.data$res.0.6order)



### visualize known markers
knownmarkers=c("VIM",
  "CD163","S100A4","TYROBP","LYZ","RGS1",
  "VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN",
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","PDGFRA","PDGFRB","PDGFB",
  "DLK1","HSD17B3","STAR","CYP17A1","NR5A1","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "CLU","SOX9","AMH"  )
length(knownmarkers) # 34
knownmarkers[which(!(knownmarkers %in% rownames(testis@data)))] 
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(testis@data))]
length(knownmarkers) # 34
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
length(knownmarkers) # 16

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
knownmarkers[which(!(knownmarkers %in% rownames(testis@data)))] #[1] "PDGFRA" "STAR"   "CASP3"
knownmarkers=knownmarkers[which(knownmarkers %in% rownames(testis@data))]
length(knownmarkers) # 28

testis=SetAllIdent(testis,id="species")
pdf(paste0(dgefile,"knownmarkers_IntProg_Feature.pdf"),height=12,width=21)
for(i in species){
  testis1=SubsetData(testis,ident.use=i)
  FeaturePlot(object = testis1, features.plot = knownmarkers,col=c("grey80","red"),nCol=7)
}
dev.off()


pdf(paste0(dgefile"knownmarkers_IntProg_Violin.pdf"),height=10,width=22)
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


### PCA and tSNE plot
testis=SetAllIdent(testis,id="orig.ident")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-1235-NoC5Somatic_Rep.pdf"),width=10.5,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,colors.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
testis=SetAllIdent(testis,id="species")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-1235-NoC5Somatic_species.pdf"),width=10,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
testis=SetAllIdent(testis,id="protocol")
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-1235-NoC5Somatic_Subject.pdf"),width=10,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
testis=SetAllIdent(testis,id=res[resi])
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Merged4-1235-NoC5Somatic_orderedclusters.pdf"),width=9.5,height=8)
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
pos=c("bottomleft","bottomleft","bottomright","bottomright")
# testis4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-2.7,3.2),c(-2.8,3.2),c(-7,6),c(-37,42))
ylims[[1]]=list(c(-4.1,2.9),c(-4.9,4.3),c(-7,6.5),c(-40,39))

dim=list(testis@dr$cca.aligned@cell.embeddings[,1:2],testis@dr$cca.aligned@cell.embeddings[,c(1,3)],testis@dr$umap@cell.embeddings[,c(1,2)],testis@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

### plot PCs and tSNE for each batch using the other batches as background
testis=SetAllIdent(testis,id="protocol")
sets=levels(testis@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette

pdf(paste0(dgefile,"OrigSubject_ACC_UMAP_tSNE.pdf"),height=6.9,width=9.2)
par(mfrow=c(3,4),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(testis@ident)
names(ident)=names(testis@ident)
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
testis=SetAllIdent(testis,id="species")
sets=levels(testis@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

pdf(paste0(dgefile,"OrigSpecies_ACC_UMAP_tSNE.pdf"),height=2.3,width=6.9)
par(mfrow=c(1,3),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
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
}
dev.off()

### plot PCs and tSNE for each batch using the other batches as background
### visualize 11 clusters
res=paste0("res.0.",3:7,"order");
resi=3
testis=SetAllIdent(testis,id=res[resi])
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


source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_testis/Rcode_multiplot.R")
source("C:/Users/qzm/Desktop/malek/Scripts_Clustering/Rcode_multiplot.R")
cols=myBrewerPalette1
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())



res=paste0("res.0.",3:7,"order");
resi=3
testis=SetAllIdent(testis,id=res[resi])
table(testis@ident)
sets=species

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
cells.use=names(testis@ident)[which(testis@meta.data$species==set)]
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

pdf(file=paste0(dgefile,"OrigSpecies_11clusters_",label,"_bg.pdf"),height=3,width=9)
multiplot(plotset,cols = 3)
dev.off()

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
cells.use=names(testis@ident)[which(testis@meta.data$species==set)]
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

pdf(file=paste0(dgefile,"OrigSpecies_11clusters_",label,"_bg.pdf"),height=3,width=9)
multiplot(plotset,cols = 3)
dev.off()



### Visualize somatic cell types from single-species somatic subclustering
testis=SetAllIdent(testis,id="SingleSpeciesClusters")
testis@ident=factor(testis@ident,levels(testis@meta.data$SingleSpeciesClusters))
table(testis@ident)
sets=species
which(rownames(data)!=names(testis@ident))


cols=myBrewerPalette1
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())


source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_testis/Rcode_multiplot.R")
source("C:/Users/qzm/Desktop/malek/Scripts_Clustering/Rcode_multiplot.R")

j=3
label="UMAP"
data=data_bg=as.data.frame(testis@dr$umap@cell.embeddings)
data$ident=testis@ident
data %>% dplyr::group_by(ident) %>% summarize(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)) -> centers3
label.size=4

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
centers=centers3[((i-1)*7+1):(i*7),]
# plot each cell type separately, adding all others as background
cells.use=names(testis@ident)[which(testis@meta.data$species==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=UMAP1,y=UMAP2,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=myBrewerPalette,limits=levels(data$ident)[((i-1)*7+1):(i*7)],drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}

pdf(file=paste0(dgefile,"OrigSpecies_SingleSpecies7clusters_",label,"_bg.pdf"),height=3,width=9)
multiplot(plotset,cols = 3)
dev.off()

j=4
label="tSNE"
data=data_bg=as.data.frame(testis@dr$tsne@cell.embeddings)
data$ident=testis@ident
data %>% dplyr::group_by(ident) %>% summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) -> centers3
label.size=4

plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
centers=centers3[((i-1)*7+1):(i*7),]
# plot each cell type separately, adding all others as background
cells.use=names(testis@ident)[which(testis@meta.data$species==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=myBrewerPalette,limits=levels(data$ident)[((i-1)*7+1):(i*7)],drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}

pdf(file=paste0(dgefile,"OrigSpecies_SingleSpecies7clusters_",label,"_bg.pdf"),height=3,width=9)
multiplot(plotset,cols = 3)
dev.off()


