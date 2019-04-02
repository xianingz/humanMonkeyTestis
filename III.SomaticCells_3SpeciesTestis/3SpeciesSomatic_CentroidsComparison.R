# Compare somatic cell type centroids from adult testis of 3 species - human, monkey and mouse
# 3.27.2019 by Qianyi
### 4 types of somatic cell type centroids:
#1) somatic cell type centroids from each of 3 species
#2) somatic cell type centroids centered by mean for each of 3 species
#3) somatic cell type centroids centered by median for each of 3 species
#4) somatic cell type centroids centered by mean and standardized for each of 3 species
### 4 sets of genes for somatic cells:
#1) union of top 6k highly-variable genes for somatic subset from each of 3 species
#2) intercept of top 6k highly-variable genes for somatic subset from each of 3 species
#3) highly-variable genes used in CCA for somatic cells of 3 species
#4) union of top 50 markers for each somatic cell type of 3 species


home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/"
setwd(home)
dgefile="3SpeciesComparison_plot/Centroids_"
redblue100<-rgb(read.table('redblue100.txt',sep='\t',row.names=1,header=T))
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(cowplot)

## load Seurat object for CCA-merged somatic cells of 3 species
### corrected for batch effect for each individual subject
load(file=paste0("Somatic_3Species_11Subjects_CCA.Robj"))
dgeall=pbmc

## load Seurat object for 3 species
load(file=paste0("human/Human4-1235-NoC5SomaticNoC1_CCA.Robj"))
dgeH=pbmc
load(file=paste0("monkey/Monkey5Somatic-NoC1-NoC7-No10pctMT_CCA.Robj"))
dgeM=testis
load(file=paste0("MouseAdultST25Somatic.Robj"))
dgem=dge

dgelist=list(dgeH,dgeM,dgem)

## load cluster centroids for 3 species
centroidH=read.table("human/CCA4-1235-NoC5SomaticNoC1_7celltypes_Centroid.txt",header=T,row.names=1)
centroidM=read.table("monkey/Merged5Somatic-NoC1-NoC7-No10pctMT_CCA_8celltypes_Centroid.txt",header=T,row.names=1)
centroidm=read.table("mouse/MouseAdultST25_11celltypes_centroid_UMI20cell15genes.txt",header=T,row.names=1)
centroidH=centroidH[,c(1:3,6,7,5,4)]
centroidm=centroidm[,c(1:4,7,5)] # removed sertoli cell type as it is not present in human or monkey
colnames(centroidH)=paste0("H_",colnames(centroidH))
colnames(centroidM)=paste0("M_",colnames(centroidM))
colnames(centroidm)=paste0("m_",colnames(centroidm))
colnames(centroidm)[5]="m_IntProg"
centroids=list(centroidH,centroidM,centroidm)


## check 1-1-1 gene orthologues present in 3-species data
Orth1to1to1=read.table("humMacMouOrth.1to1to1.txt",stringsAsFactors=F)

allgenes=Orth1to1to1
for(i in 1:length(centroids)){
  centroid=centroids[[i]]
  allgenes=allgenes[which(allgenes[,i] %in% rownames(centroid)),]
}
dim(allgenes) # [1] 10768      3
for(i in 1:length(centroids)){
  print(length(unique(allgenes[,i])))
}


### keep 1-1-1 ortholgue genes for 3 species cluster centroids data
for(i in 1:length(centroids)){
  centroid=centroids[[i]]
  centroid=centroid[which(rownames(centroid) %in% allgenes[,i]),]
  centroids[[i]]=centroid
  print(nrow(centroid))
}

### replace gene names in monkey and mouse by human gene names for cluster centroids data
for(i in 2:3){
	centroid=centroids[[i]]
  	for(g in 1:nrow(centroid)){
  		rownames(centroid)[g]=allgenes[which(allgenes[,i]==rownames(centroid)[g]),1]
  	}
  	centroids[[i]]=centroid
}	


### merge cluster centroids data for 3 species together
set=list()
for(i in 1:length(centroids)){
    set[[i]]=data.frame(GENE=rownames(centroids[[i]]),centroids[[i]])
}
allcentroids=Reduce(function(x,y) merge(x,y,all=TRUE), set)
rownames(allcentroids)=allcentroids[,1]
allcentroids=allcentroids[,-1]
dim(allcentroids)   # [1] 10768    21


### center cluster centroids for each species by mean separately
centroids.cen=list()
for(i in 1:3){
    centroid=centroids[[i]]
    centroid.cen=centroid-apply(centroid,1,mean)
    centroids.cen[[i]]=centroid.cen
} 

set=list()
for(i in 1:length(centroids.cen)){
    set[[i]]=data.frame(GENE=rownames(centroids.cen[[i]]),centroids.cen[[i]])
}
allcentroids.cen=Reduce(function(x,y) merge(x,y,all=TRUE), set)
rownames(allcentroids.cen)=allcentroids.cen[,1]
allcentroids.cen=allcentroids.cen[,-1]
dim(allcentroids.cen)   


### center cluster centroids for each species by median separately
centroids.med=list()
for(i in 1:3){
    centroid=centroids[[i]]
    centroid.med=centroid-apply(centroid,1,median)
    centroids.med[[i]]=centroid.med
} 

set=list()
for(i in 1:length(centroids.med)){
    set[[i]]=data.frame(GENE=rownames(centroids.med[[i]]),centroids.med[[i]])
}
allcentroids.med=Reduce(function(x,y) merge(x,y,all=TRUE), set)
rownames(allcentroids.med)=allcentroids.med[,1]
allcentroids.med=allcentroids.med[,-1]
dim(allcentroids.med)   

### center and standardize cluster centroids for each species 
centroids.std=list()
for(i in 1:3){
    centroid=centroids[[i]]
    centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
    centroids.std[[i]]=centroid.std
} 

set=list()
for(i in 1:length(centroids.std)){
    set[[i]]=data.frame(GENE=rownames(centroids.std[[i]]),centroids.std[[i]])
}
allcentroids.std=Reduce(function(x,y) merge(x,y,all=TRUE), set)
rownames(allcentroids.std)=allcentroids.std[,1]
allcentroids.std=allcentroids.std[,-1]
dim(allcentroids.std)   # [1] 10768    21
allcentroids.std[is.na(allcentroids.std)]<-0



## Top 6k Highly-variable genes selected from each of 3 species somatic
hvglist=list()
    dgelist[[3]]@calc.params=dgelist[[1]]@calc.params
for(i in 1:length(dgelist)){
  dge=dgelist[[i]]
  dge <- FindVariableGenes(object = dge, do.plot = FALSE)
  hvglist[[i]] <- rownames(x = head(x = dge@hvg.info, n = 6000))
  print(length(hvglist[[i]])) 
}
#[1] 6000
#[1] 6000
#[1] 6000

### keep 1-1-1 ortholgue genes for 3 species hvg
for(i in 1:length(hvglist)){
  hvg=hvglist[[i]]
  hvg=hvg[which(hvg %in% allgenes[,i])]
  hvglist[[i]]=hvg
  print(length(hvg))
}
#[1] 1992
#[1] 4122
#[1] 2996

### replace gene names in monkey and mouse by human gene names for hvg
for(i in 2:3){
  hvg=hvglist[[i]]
  for(g in 1:length(hvg)){
    hvg[g]=allgenes[which(allgenes[,i]==hvg[g]),1]
  }
  hvglist[[i]]=hvg
} 

### Using union of Highly-variable genes from each of 3 species somatic
hvg.union=unique(unlist(hvglist))
hvg.union=intersect(hvg.union,rownames(allcentroids))
print(length(hvg.union)) 
#[1] 6124


### Using intersect of Highly-variable genes from each of 3 species somatic
hvg.intersect=intersect(hvglist[[1]],hvglist[[2]])
hvg.intersect=intersect(hvg.intersect,hvglist[[3]])
hvg.intersect=unique(hvg.intersect)
print(length(hvg.intersect)) 
# 673


## Highly-variable genes selected for CCA-combined all 3 species somatic 
hvg=dgeall@var.genes[which(dgeall@var.genes %in% rownames(allcentroids))]


## Top Markers from Each of 3 species somatic
markerlist=list()
markersH=read.table("human/NoC5NoC1CCA4Somatic_7celltypes_markersall_mindiff0.2_logfc2fold_3.13.2019.txt",header=T,row.names=1)
markersM=read.table("monkey/NoC1NoC7No10pctMTCCA5Somatic_8celltypes_markersall_mindiff0.2_logfc2fold_3.27.2019.txt",header=T,row.names=1)
markersm=read.table("mouse/MouseAdultST25Somatic_MarkersAll_7celltypes_pct0.2_diffpct0.2_thresh2fold_1.8.2018.txt",header=T,row.names=1)
markersm=markersm[,c(1:4,1,5:6)]
colnames(markersm)=colnames(markersH)
markersall=list(markersH,markersM,markersm)

top50=top50markers=list()
for(i in 1:length(markersall)){
  markersall[[i]] %>% group_by(cluster) %>% top_n(50, p_val)  -> top50[[i]]
  top50markers[[i]]=as.character(top50[[i]]$gene)
  print(table(top50[[i]]$cluster))
  print(length(unique(top50[[i]]$gene)))
  print(length(unique(top50markers[[i]])))
}
#Endothelial  f-Pericyte   ImmLeydig  m-Pericyte  Macrophage       Myoid 
#         50          50          42          50          50          32 
#      Tcell 
#         50 
#[1] 303

# DiffLeydig Endothelial  f-Pericyte   ImmLeydig  m-Pericyte  Macrophage 
#         50          50          50          40          50          50 
#      Myoid       Tcell 
#         31          50 
#[1] 333
#[1] 333

#   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
#            50             50             50             50             50 
#       Sertoli        Unknown 
#            50             50 
#[1] 337


### keep 1-1-1 ortholgue genes for 3 species somatic marker
for(i in 1:length(top50markers)){
  marker=top50markers[[i]]
  marker=marker[which(marker %in% allgenes[,i])]
  top50markers[[i]]=marker
  print(length(marker))
}
#[1] 215
#[1] 268
#[1] 245

### replace gene names in monkey and mouse by human gene names for marker
for(i in 2:3){
  marker=top50markers[[i]]
  for(g in 1:length(marker)){
    marker[g]=allgenes[which(allgenes[,i]==marker[g]),1]
  }
  top50markers[[i]]=marker
} 

### Using union of markers from each of 3 species somatic
markers=unique(unlist(top50markers))
markers=intersect(markers,rownames(allcentroids))
length(markers) 
#[1] 575

length(which(markers %in% hvg.intersect)) # 188



## match colors of similar cell types
colnames(centroidH)
colnames(centroidM)
colnames(centroidm)
#[1] "H_Tcell"       "H_Macrophage"  "H_Endothelial" "H_m.Pericyte" 
#[5] "H_f.Pericyte"  "H_Myoid"       "H_ImmLeydig"  
#[1] "M_Tcell"       "M_Macrophage"  "M_Endothelial" "M_m.Pericyte"   
#[5] "M_f.Pericyte" "M_Myoid"       "M_ImmLeydig"   "M_DiffLeydig" 
#[1] "m_InnateLymphoid" "m_Macrophage"     "m_Endothelial"    "m_Myoid"         
#[5] "m_IntProg"        "m_Leydig" 

library(RColorBrewer)
colH=brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3)]
colM=brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3,8)]
colm=brewer.pal(8,"Set1")[c(7,6,4,5,3,8)]
#colH=myBrewerPalette[c(2:1,4,10,6,5,9)]
#colM=myBrewerPalette[c(1,4,10,5,6,8,7)]
#colm=myBrewerPalette[c(2:1,4,6,8,10,12)]
cols=list(colH,colM,colm)
ncluster=c(ncol(centroidH),ncol(centroidM),ncol(centroidm))

col.use=redblue100

## Sets of centroids data and genes
cens=list(allcentroids,allcentroids.cen,allcentroids.med,allcentroids.std)
tag1s=c("","_SpeciesCenteredbyMean","_SpeciesCenteredbyMed","_SpeciesScaled")
genes=list(hvg.union,hvg.intersect,hvg,markers)
tag2s=c("HVGunion","HVGintersect","HVG","Top50markers")

library(rgl)
library(car)

for(i in 1:length(cens)){
  cen=cens[[i]]
  tag1=tag1s[i]
  for(j in 1:length(genes)){
    genes.use=genes[[j]]
    tag2=tag2s[j]

## PCA
pp=prcomp(t(cen[genes.use,]))


## Rank correlation for cluster centroids using HVG for 3 species
cc=cor(as.matrix(cen)[genes.use,],method="spearman")
dim(cc) # [1] 21 21
min(cc) # [1] 0.1859652
levels=colnames(cen)

### labeling
data.use=cc[levels,levels]
colsep.use=cumsum(table(gsub("_.*","",levels))[unique(gsub("_.*","",levels))])
row.lab=gsub(".*_","",levels)
col.lab=gsub("_.*","",levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=c(rep("white",ncluster[1]),rep("grey50",ncluster[2]),rep("black",ncluster[3]))
sidecol[2,]=unlist(cols)
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("Species","Cell Type")
colnames(clab)=c("Cell Type","Species")


## plot rank correlation
pdf(file=paste0(dgefile,"nrho_3SpeciesSomatic",tag1,"_",tag2,".pdf"),height=6,width=6)

### set diagonal values to NA
data.use0=data.use
diag(data.use0) <- NA 
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0)) # ,bg="grey"
heatmap.3(data.use0,na.col="grey",dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,6))

par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0),bg="white")
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=.9,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,6))

dev.off()


## visualize PCA

pdf(file=paste0(dgefile,"PCA0_3SpeciesSomatic",tag1,"_",tag2,".pdf"),height=5,width=5)
### color each cell type by different colors; use different symbol for different species
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],2],xlim=c(min(pp$x[,1])-2.5,max(pp$x[,1])+2),ylim=c(min(pp$x[,2])+.2,max(pp$x[,2])+.2),xlab="PC1",ylab="PC2",pch=15,col=colH)
points(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),2],pch=16,col=colM)
points(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),2],pch=17,col=colm)
text(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],2]+.2,labels=colnames(allcentroids)[1:ncluster[1]],col=colH)
text(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),2]+.2,labels=colnames(allcentroids)[(ncluster[1]+1):sum(ncluster[1:2])],col=colM)
text(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),2]+.2,labels=colnames(allcentroids)[(sum(ncluster[1:2])+1):sum(ncluster)],col=colm)

### use different colors for different species
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],2],xlim=c(min(pp$x[,1])-2.5,max(pp$x[,1])+2),ylim=range(pp$x[,2]),xlab="PC1",ylab="PC2",pch=16,col=rgb(0,0,0,0.5))
points(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),2],pch=16,col=rgb(0,0,1,0.5))
points(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),2],pch=16,col=rgb(1,0,0,0.5))
text(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],2],labels=colnames(allcentroids)[1:ncluster[1]])
text(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),2],labels=colnames(allcentroids)[(ncluster[1]+1):sum(ncluster[1:2])],col="blue")
text(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),2],labels=colnames(allcentroids)[(sum(ncluster[1:2])+1):sum(ncluster)],col="red")

dev.off()


### modify PCA plots to aid better visualization
pdf(file=paste0(dgefile,"PCA1_3SpeciesSomatic",tag1,"_",tag2,".pdf"),height=3.5,width=3.5)
par(mar=c(4,4,1,1),mgp=c(2, .5, 0))
### color each cell type by different colors; 
### use different symbol for different species
plot(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],2],xlim=c(min(pp$x[,1])-2.5,max(pp$x[,1])+2),ylim=c(min(pp$x[,2])+.2,max(pp$x[,2])+.2),xlab="PC1",ylab="PC2",pch=15,col=colH)
points(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),2],pch=16,col=colM)
points(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),2],pch=17,col=colm)

plot(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],3],xlim=c(min(pp$x[,1])-2.5,max(pp$x[,1])+2),ylim=c(min(pp$x[,3])+.2,max(pp$x[,3])+.2),xlab="PC1",ylab="PC3",pch=15,col=colH)
points(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),3],pch=16,col=colM)
points(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),3],pch=17,col=colm)

plot(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],2],xlim=c(min(pp$x[,1])-2.5,max(pp$x[,1])+2),ylim=c(min(pp$x[,2])+.2,max(pp$x[,2])+.2),xlab="PC1",ylab="PC2",pch=0,col=colH)
points(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),2],pch=1,col=colM)
points(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),2],pch=2,col=colm)

### minimize text: label each cell type only once for all 3 species
data=data.frame(pp$x[,1:2])
data$species=c(rep(1,7),rep(2,8),rep(3,6))
data$ident=rownames(data)
data$celltype=gsub(".*_","",rownames(data))
data$celltype[which(data$celltype=="InnateLymphoid")]<-"Tcell"
data$celltype[which(data$celltype=="Myoid" & data$species==3)]<-"m.Pericyte"
data$celltype[which(data$celltype=="ImmLeydig")]<-"IntProg/ImmLeydig"
data$celltype[which(data$celltype=="IntProg")]<-"IntProg/ImmLeydig"
data$celltype[which(data$celltype=="DiffLeydig")]<-"Leydig"
data$celltype=factor(data$celltype,levels=c("Tcell","Macrophage","Endothelial","m.Pericyte","f.Pericyte","Myoid","IntProg/ImmLeydig","Leydig"))
table(data$celltype)
data %>% dplyr::group_by(celltype) %>% summarize(PC1 = median(PC1), PC2 = median(PC2)) -> centers

plot(pp$x[1:ncluster[1],1],pp$x[1:ncluster[1],2],xlim=c(min(pp$x[,1])-2.5,max(pp$x[,1])+2),ylim=c(min(pp$x[,2])+.2,max(pp$x[,2])+.2),xlab="PC1",ylab="PC2",pch=15,col=alpha(colH,0.5))
points(pp$x[(ncluster[1]+1):sum(ncluster[1:2]),1],pp$x[(ncluster[1]+1):sum(ncluster[1:2]),2],pch=16,col=alpha(colM,0.5))
points(pp$x[(sum(ncluster[1:2])+1):sum(ncluster),1],pp$x[(sum(ncluster[1:2])+1):sum(ncluster),2],pch=17,col=alpha(colm,0.5))
for(i in 1:nrow(centers)){
text(as.numeric(centers[i,2]),as.numeric(centers[i,3]),levels(data$celltype)[i],col=colM[i])
}

### add ellipse
dataEllipse(data[,1],data[,2],xlab="PC1",ylab="PC2",col=colM,
  groups=data$celltype,levels=0.5, fill=TRUE, fill.alpha=0.1)
#scatter3d(data[,1],data[,2],data[,3], groups=data$celltype, surface=FALSE,surface.col=colM, ellipsoid = TRUE)

dataEllipse(data[,1],data[,2],xlab="PC1",ylab="PC2",col=colM,
  groups=data$celltype,group.labels=NULL,
  center.cex=0,
  levels=0.5, fill=TRUE, fill.alpha=0.1)
#pch=c(15:17)[data$species],

dataEllipse(data$PC1,data$PC2,xlab="PC1",ylab="PC2",col=colM,
  groups=data$celltype,group.labels=NULL,
  center.cex=0,cex=0,
  levels=0.5, fill=TRUE, fill.alpha=0.1)
points(data[,1],data[,2],pch=c(2:4)[data$species],col=colM[as.numeric(data$celltype)])

dataEllipse(data$PC1,data$PC2,xlab="PC1",ylab="PC2",col=colM,
  groups=data$celltype,group.labels=NULL,
  center.cex=0,cex=0,
  levels=0.5, fill=TRUE, fill.alpha=0.1)
points(data[,1],data[,2],pch=c(15:17)[data$species],col=alpha(colM[as.numeric(data$celltype)],0.8))


dev.off()
}
}



i=2;j=2;
  cen=cens[[i]]
  tag1=tag1s[i]
    genes.use=genes[[j]]
    tag2=tag2s[j]
pp=prcomp(t(cen[genes.use,]))
eigs <- pp$sdev^2
v1=eigs[1] / sum(eigs)
v2=eigs[2] / sum(eigs)
xlab=paste0("PC1 (",sprintf("%1.1f%%", 100*v1),")")
ylab=paste0("PC2 (",sprintf("%1.1f%%", 100*v2),")")

data=data.frame(pp$x[,1:2])
data$species=c(rep(1,7),rep(2,8),rep(3,6))
data$ident=rownames(data)
data$celltype=gsub(".*_","",rownames(data))
data$celltype[which(data$celltype=="InnateLymphoid")]<-"Tcell"
data$celltype[which(data$celltype=="Myoid" & data$species==3)]<-"m.Pericyte"
data$celltype[which(data$celltype=="ImmLeydig")]<-"IntProg/ImmLeydig"
data$celltype[which(data$celltype=="IntProg")]<-"IntProg/ImmLeydig"
data$celltype[which(data$celltype=="DiffLeydig")]<-"Leydig"
data$celltype=factor(data$celltype,levels=c("Tcell","Macrophage","Endothelial","m.Pericyte","f.Pericyte","Myoid","IntProg/ImmLeydig","Leydig"))
table(data$celltype)
data %>% dplyr::group_by(celltype) %>% summarize(PC1 = median(PC1), PC2 = median(PC2)) -> centers

write.table(data,paste0("3SpciesSomaticCentroids",tag1,"_",tag2,".txt"),row.names=T,col.names=T,sep="\t",quote=F)

data=read.table(paste0("3SpciesSomaticCentroids",tag1,"_",tag2,".txt"),row.names=1,header=T)
data$celltype=factor(data$celltype,levels=c("Tcell","Macrophage","Endothelial","m.Pericyte","f.Pericyte","Myoid","IntProg/ImmLeydig","Leydig"))
### modify xlim and ylim

### modify xlim and ylim
pdf(file=paste0(dgefile,"PCA2_3SpeciesSomatic",tag1,"_",tag2,".pdf"),height=3.5,width=3.5)
par(mar=c(4,4,1,1),mgp=c(2, .5, 0))
dataEllipse(data$PC1,data$PC2,xlab=xlab,ylab=ylab,col=colM,
  xlim=c(-15,11),ylim=c(-9,8.5),
  groups=data$celltype,group.labels=NULL,
  center.cex=0,cex=0,
  levels=0.5, fill=TRUE, fill.alpha=0.5)
points(data[,1],data[,2],pch=c(15:17)[data$species],col=colM[as.numeric(data$celltype)])
dataEllipse(data$PC1,data$PC2,xlab=xlab,ylab=ylab,col=colM,
  xlim=c(-15,11),ylim=c(-9,8.5),
  groups=data$celltype,group.labels=NULL,
  center.cex=0,cex=0,lwd=0.1,
  levels=0.5, fill=TRUE, fill.alpha=0.5)
points(data[,1],data[,2],pch=c(15:17)[data$species],col=colM[as.numeric(data$celltype)])
dataEllipse(data$PC1,data$PC2,xlab=xlab,ylab=ylab,col=colM,
  xlim=c(-15,11),ylim=c(-9,8.5),
  groups=data$celltype,group.labels=NULL,
  center.cex=0,cex=0,lwd=0,
  levels=0.5, fill=TRUE, fill.alpha=0.5)
points(data[,1],data[,2],pch=c(15:17)[data$species],col=colM[as.numeric(data$celltype)])
dataEllipse(data$PC1,data$PC2,xlab=xlab,ylab=ylab,col=colM,
  xlim=c(-15,11),ylim=c(-9,8.5),
  groups=data$celltype,group.labels=NULL,
  center.cex=0,cex=0,
  levels=0.5, fill=TRUE, fill.alpha=0.3)
points(data[,1],data[,2],pch=c(15:17)[data$species],col=alpha(colM[as.numeric(data$celltype)],0.8))
dataEllipse(data$PC1,data$PC2,xlab=xlab,ylab=ylab,col=colM,
  xlim=c(-15,11),ylim=c(-9,8.5),
  groups=data$celltype,group.labels=NULL,
  center.cex=0,cex=0,
  levels=0.5, fill=TRUE, fill.alpha=0.3)
points(data[,1],data[,2],pch=c(2,3,4)[data$species],col=colM[as.numeric(data$celltype)])
dev.off()



plot.new()
legend("topleft",cex=2,pch=15:17,c("Human","Monkey","Mouse"))

plot.new()
library(RColorBrewer)
display.brewer.pal(name="Paired",12)

plot.new()
legend("topleft",cex=1,pch=16,col=colM,legend=levels(data$celltype))
