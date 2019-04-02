### Monkey Pericyte subclustering
# 3.6.2019 by Qianyi
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering 
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1 cells after removing C1 doublets
### Removed Cluster 7/7 of SCyte-Somatic doublets from somatic subclustering by CCA
### Re-did CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing C1 and C7 doublets
### Extract Pericyte, directly-merged and do subclustering
#-> 2 subclusters (Pericytes + Subcluster2)
# note: here I do not use CCA for Pericyte subclustering as there are too few cells for each subject

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


# setup Seurat objects 
### remove 7/7 from previous CCA and somatic subclustering
table(testis1@ident)
#Macrophage/Tcell      Endothelial         Pericyte            Myoid 
#              37              170              174              516 
#       ImmLeydig       DiffLeydig 
#            1396               44 
testis=SubsetData(testis1,ident.use="Pericyte")
table(testis@ident)
#Pericyte 
#     174 
table(testis@ident,testis@meta.data$protocol)
#           Monkey1 Monkey2 Monkey3 Monkey4 Monkey5
#  Pericyte     140       6      15      10       3
nCellperGene <- rowSums(testis@data>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # [1] 10761

datalist=list()
for(i in 1:length(subject)){
### Extract Pericyte for each monkey subject
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
 10761 genes across 140 samples.

[[2]]
An object of class seurat in project SeuratProject 
 10761 genes across 6 samples.

[[3]]
An object of class seurat in project SeuratProject 
 10761 genes across 15 samples.

[[4]]
An object of class seurat in project SeuratProject 
 10761 genes across 10 samples.

[[5]]
An object of class seurat in project SeuratProject 
 10761 genes across 3 samples.
# note: too few cells for each subject, so cannot use CCA


dgedata2=testis@raw.data[genes.use,names(testis@ident)]
dim(dgedata2) # [1] 10761   174

### create new Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project="MonkeyPericyte_DirectMerged5Subjects", min.cells=1, min.genes=1)
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf("Monkey_Pericyte_DirectMerged5Subjects_HVG0.2.pdf",height=7.5,width=11)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) # [1] 2187
### Scale data 
dge <- ScaleData(object = dge)
dge                    # 10761 genes across 174 samples     
### Add a column of subject ID in meta data sheet
ident<-gsub("\\..*","",as.character(dge@ident))
names(ident)=names(dge@ident)
ident=as.factor(ident)
table(ident)
#Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
#    140       6      15      10       3 
dge=AddMetaData(dge,ident,"indiv")
dge=SetAllIdent(dge,"indiv")
save(dge,file="MonkeyPericyte_DirectMerged5Subjects.Robj")
table(dge@meta.data$indiv)
#Monkey1 Monkey2 Monkey3 Monkey4 Monkey5 
#    140       6      15      10       3
table(gsub("_.*","",names(dge@ident)))
#Monkey1.1 Monkey1.2   Monkey2 Monkey3.1 Monkey3.2 Monkey4.1 Monkey4.2 Monkey5.2 
#       81        59         6        10         5         7         3         3 

### PCA
print(Sys.time())  # "2019-03-06 16:27:17 EST" 
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)

### PCA Plot
pdf("PCA_Pericyte_DirectMerged5Subjects.pdf",width=12,height=10)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()
### Scree Plot for PCA
numPCs=3;i=1 # HVG
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
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,1,by=0.1), save.SNN = T)
dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
dge <- RunUMAP(dge, reduction.use = "pca", dims.use = 1:numPCs[i])
Sys.time()  # [1] "2019-03-06 16:30:30 EST"

save(dge,file="MonkeyPericyte_DirectMerged5Subjects.Robj")


### check the number of clusters
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
# 2 2 2 2 2 2 3 3 4 4

table(dge@meta.data$res.0.1,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.2,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.3,dge@meta.data$res.0.5)
table(dge@meta.data$res.0.4,dge@meta.data$res.0.5)

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
pdf(file=paste0("Pericyte_DirectMerged5SubjectsSubjects_Centroid_norm_Seriation_",res,".pdf"))
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

### label cell types for Pericyte subclustering
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident2[which(ident2==1)]<-"Pericyte"
ident2[which(ident2==2)]<-"Subcluster2"
dge=AddMetaData(dge,ident2,"CellType2")
dge=SetAllIdent(dge,id="CellType2")
save(dge,file="MonkeyPericyte_DirectMerged5Subjects.Robj")

pdf(paste0(dgefile,"PCA_tSNE_UMAP_Pericyte2.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

### Add these two clusters back to All Somatic cells of Monkey
testis=testis1
ident=as.character(testis1@ident)
names(ident)=names(testis1@ident)
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident[names(ident2)]<-ident2
levels=c(levels(testis1@ident)[1:3],"Subcluster2",levels(testis1@ident)[4:6])
table(ident)[levels]
#Macrophage/Tcell      Endothelial         Pericyte      Subcluster2 
#              37              170              131               43 
#           Myoid        ImmLeydig       DiffLeydig 
#             516             1396               44 
ident=factor(ident,levels=levels)
testis=AddMetaData(testis,ident,"CellType2")

testis=SetAllIdent(testis,id="CellType2")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType2))
testis1=testis
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA.Robj"))

### Visualize the additional Pericyte subcluster in Somatic subset
library(RColorBrewer)
myBrewerPalette1 <- c(brewer.pal(7,"Set1"))[c(1:3,7,4:6)]
myBrewerPalette2 <- c(brewer.pal(12,"Paired"))[c(11,10,6,7,9,5)] 
# used this color scheme
pdf(paste0(dgefile,"PCA_tSNE_UMAP_Somatic6celltypes1_Pericyte2.pdf"),width=11,height=8)
plot2=DimPlot(testis,reduction.use = "cca.aligned",1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot3=DimPlot(testis,reduction.use = "cca.aligned",1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette1)
plot4=DimPlot(testis,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette1)
plot5=TSNEPlot(testis,do.return=TRUE,pt.size = .8,do.label=T,colors.use=myBrewerPalette1)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

###### Visualize individual batch and subject
table(dge@meta.data$orig.ident,dge@ident)
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlims=ylims=list()
pos=c("bottomleft","topleft","topleft","bottomleft")
# dge4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-1.4,3.8),c(-1.4,3.8),c(-8,11),c(-9,12))
ylims[[1]]=list(c(-5.5,1.9),c(-3.5,4.5),c(-4.9,9.5),c(-9,7))
dim=list(dge@dr$cca.aligned@cell.embeddings[,1:2],dge@dr$cca.aligned@cell.embeddings[,c(1,3)],dge@dr$umap@cell.embeddings[,c(1,2)],dge@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

sets=subject
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())

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
cells.use=names(dge@ident)[which(gsub("\\..*","",gsub("_.*","",names(dge@ident)))==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
#scale_colour_manual(values=myBrewerPalette,limits=levels(data$ident),drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}

### plot
source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_dge/Rcode_multiplot.R")
pdf(file=paste0(dgefile,"Pericyte_OrigSubject_2clusters_",label,"_bg.pdf"),height=6,width=9)
multiplot(plotset,cols = 3)
dev.off()



### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlims=ylims=list()
pos=c("bottomleft","topleft","topleft","bottomleft")
# testis4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-1.4,3.8),c(-1.4,3.8),c(-8,11),c(-31,28))
ylims[[1]]=list(c(-5.5,1.9),c(-3.5,4.5),c(-4.9,9.5),c(-41,35))

dim=list(testis@dr$cca.aligned@cell.embeddings[,1:2],testis@dr$cca.aligned@cell.embeddings[,c(1,3)],testis@dr$umap@cell.embeddings[,c(1,2)],testis@dr$tsne@cell.embeddings[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

### plot PCs and tSNE for each cluster using the other clusters as background
sets=levels(testis@ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

pdf(paste0(dgefile,"Somatic_WithPericyte2_IndivCluster_PCtSNE.pdf"),height=4.6,width=9.2)
par(mfrow=c(2,4),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
j=4 #for(j in 1:length(dim)){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(testis@ident)
names(ident)=names(testis@ident)
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.4,col=c("grey70",alpha(cols[seti],0.99))[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
#plot.new()
#}
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


###### assign cell types based on markers for each ordered clusters in each dataset
markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Pericyte2_res.0.5order_markersall_mindiff0.2_logfc2fold_3.6.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
   Pericyte Subcluster2 
        131          43
         75          53


pdf(file=paste0("Pericyte2_PerCellAttributes_ViolinPlot.pdf"),height=4,width=4)
VlnPlot(dge, features.plot="nGene", nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(dge, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()


markers=FindAllMarkers(testis,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,"Somatic_WithPericyte2_markersall_mindiff0.2_logfc2fold_3.6.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
table(testis@ident)
table(markers$cluster)
Macrophage/Tcell      Endothelial         Pericyte      Subcluster2 
              37              170              131               43 
           Myoid        ImmLeydig       DiffLeydig 
             516             1396               44 

Macrophage/Tcell      Endothelial         Pericyte      Subcluster2 
             110              128              150               57 
           Myoid        ImmLeydig       DiffLeydig 
              25               29               90 

pdf(file=paste0("Somatic_Pericyte2_PerCellAttributes_ViolinPlot.pdf"),height=4,width=8)
VlnPlot(testis, features.plot="nGene", nCol = 1,point.size.use=-1,cols.use=myBrewerPalette1)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
VlnPlot(testis, features.plot="nUMI", y.log=T, nCol = 1,point.size.use=-1,cols.use=myBrewerPalette1)+geom_boxplot(width=0.1,outlier.size = -1)#,cols.use=myBrewerPalette)
dev.off()

markers %>% group_by(cluster) %>% top_n(5, avg_logFC)  -> top5
data.frame(top5)
markers %>% group_by(cluster) %>% top_n(2, avg_logFC)  -> top2
pdf("markerstop.pdf",height=12,width=12)
FeaturePlot(object = testis, features.plot = top2$gene, min.cutoff = "q9", cols.use = c("lightblue", 
    "red"), pt.size = .8,nCol=4)
dev.off()


markers=FindMarkers(testis,"Subcluster2","ImmLeydig",only.pos=FALSE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,paste0(dgefile,"Pericyte2VsImmLeydig_markers_mindiff0.2_logfc2fold_3.6.2019.txt"),col.names=T,row.names=T,quote=F,sep="\t")
length(which(markers$avg_logFC>0)) # 77
length(which(markers$avg_logFC<0)) # 86


