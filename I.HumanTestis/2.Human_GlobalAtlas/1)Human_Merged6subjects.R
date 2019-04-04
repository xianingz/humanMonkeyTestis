### Global Clustering for Directly-merged all 6 human subjects
# by Qianyi on 11.30.2018
### visualized each individual batch in global PCs and t-SNE, Examined batch effect and markers



# 11.30.2018 by Qianyi
# clustering for merged 6 human subject, order clusters
### set paths and load library
home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/human"
setwd(home)
library(Seurat)
library(RColorBrewer)
myBrewerPalette=brewer.pal(12,"Paired")
redblue100<-rgb(read.table('../redblue100.txt',sep='\t',row.names=1,header=T))
datainfo=read.table("datainfo_human")
datainfo
dataset=as.character(datainfo[,2])
subject=unique(gsub("-.*","",dataset))
### load object for merged 6 subjects
load(file="HumanMerged6.Robj")
dge6=dgeall
dge=dgeall
dge6 # 48655 genes across 22663 samples
### Number of highly-variable genes
print(length(dge@var.genes)) # 3771
### Add per-cell attributes
##### percent.mito 
mito.genes <- grep("^MT-", rownames(dge@data), value = T)
length(mito.genes) # 32 dgeall, 27 dge25
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

### PCA
print(Sys.time())   # [1] "2018-11-29 17:30:31 EST"
dge <- PCA(dge, do.print = TRUE, pcs.print = 1, genes.print = 5)
save(dge,file="HumanMerged6.Robj")
print(Sys.time())   # [1] "2018-11-29 17:30:31 EST"
### PCA Plot
pdf("PCA_Merged6Human.pdf",width=15,height=12)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
dev.off()
### Scree Plot for PCA
numPCs=23;i=1 # HVG
pdf("dge_PCA_Variablel_variation.pdf",height=5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:200],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
plot(dge@pca.obj[[1]]$sdev[1:60],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
eigenvalue # 1.909231
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),xlim=c(0,6),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
plot(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()
### Louvain-Jaccard clustering and tSNE
Sys.time() # [1] "2018-11-30 09:18:27 EST"
dge=RunTSNE(dge,dims.use = 1:numPCs[i],do.fast=T)    # max_iter=2000
Sys.time() # [1] "2018-11-30 10:13:47 EST"
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(0.1,3,0.1), print.output = 0, save.SNN = T)
Sys.time() 
dgeall=dge
save(dgeall,file="HumanMerged6.Robj")


# 12.5.2018 Qianyi
# clustering for merged 6 human subject, order clusters
### load object for merged 6 subjects
###### top 23 PCs
load(file="HumanMerged6.Robj")
dge6=dgeall
dge=dgeall
dge6 # 48655 genes across 22663 samples
### Number of highly-variable genes
print(length(dge@var.genes)) # 3771
### check the number of hvg overlapped among subjects
    for(j in 1:6){
        print(length(intersect(dge6@var.genes,dgealllist[[j]]@var.genes)))
    }
[1] 1538
[1] 1397
[1] 1865
[1] 1396
[1] 1592
[1] 1722


### Scree Plot for PCA
numPCs=14;i=1 # HVG
pdf(paste0("dge_PCA_Variablel_variation_",numPCs,"PCs.pdf"),height=5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:200],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
plot(dge@pca.obj[[1]]$sdev[1:60],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
eigenvalue # 1.909231
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),xlim=c(0,6),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
plot(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()
### tSNE
Sys.time() # [1] "2018-12-05 10:22:22 EST"
dge=RunTSNE(dge,dims.use = 1:numPCs[i],do.fast=T)    # max_iter=2000
Sys.time() # [1] "2018-12-05 10:26:12 EST"
dgeall=dge
save(dgeall,file="HumanMerged6_14PCs.Robj")


library(dplyr)
### PCA and tSNE plot
dge=SetAllIdent(dge,id="orig.ident")
pdf("PCA_tSNE_Merged6Human_Rep.pdf",width=11,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
plot2=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,6,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,7,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,8,do.return = TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
dev.off()
dge=SetAllIdent(dge,id="indiv")
pdf("PCA_tSNE_Merged6Human_Subject.pdf",width=11,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
plot2=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,6,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,7,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,8,do.return = TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
dev.off()
dge=SetAllIdent(dge,id="res.0.2")
pdf("PCA_tSNE_Merged6Human_12clusters.pdf",width=10,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
plot2=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,6,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,7,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,8,do.return = TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
dev.off()

### plot PCs and tSNE for each batch using the other batches as background
xlims=ylims=list()
pos=c("topright","bottomright","topleft","topleft")
# dge6 used top 23 PCs for tSNE
i=1
xlims[[1]]=list(c(-38,22),c(-38,22),c(-38,22),c(-45,50))
ylims[[1]]=list(c(-32,16),c(-23,22),c(-19,14),c(-48,50))
# dge6 used top 14 PCs for tSNE
i=2
xlims[[2]]=list(c(-38,22),c(-38,22),c(-38,22),c(-48,48))
ylims[[2]]=list(c(-32,16),c(-23,22),c(-19,14),c(-49,50))


### plot PCs and tSNE for each batch using the other batches as background

dge=SetAllIdent(dge,id="orig.ident")
sets=levels(dge@data.info$orig.ident)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette

dim=list(dge@pca.rot[,1:2],dge@pca.rot[,c(1,3)],dge@pca.rot[,c(1,4)],dge@tsne.rot[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

pdf("OrigSet_PCtSNE.pdf",height=7,width=11.5)
par(mfrow=c(3,5),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(dge@ident)
names(ident)=names(dge@ident)
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.8,col=c("grey70",cols[seti])[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
  plot.new()
  plot.new()
}
dev.off()

### plot PCs and tSNE for each subject using the other subjects as background
dge=SetAllIdent(dge,id="indiv")
sets=levels(dge@data.info$indiv)
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette

dim=list(dge@pca.rot[,1:2],dge@pca.rot[,c(1,3)],dge@pca.rot[,c(1,4)],dge@tsne.rot[,1:2])
xlim=xlims[[i]]
ylim=ylims[[i]]

pdf("OrigSubject_PCtSNE.pdf",height=4.6,width=7)
par(mfrow=c(2,3),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:4){
for(seti in 1:length(sets)){
set=sets[seti]
ident=as.character(dge@ident)
names(ident)=names(dge@ident)
ident[which(!(ident %in% set))] <- "Others"
ident=factor(ident,levels=c("Others",set))
ident=sort(ident)
tmp=dim[[j]][names(ident),]
plot(tmp,pch=16,cex=0.8,col=c("grey70",cols[seti])[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}
}
dev.off()

### plot PCs and tSNE for each batch using the other batches as background
### visualize 12 clusters
dge=SetAllIdent(dge,id="res.0.2")
table(dge@data.info$indiv,dge@data.info$res.0.2)

  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@data.info$orig.ident,dge@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,"ncellspercluster_batch_top23PCs.txt",quote=F,row.names=T,col.names=T,sep="\t")


sets=levels(dge@data.info$orig.ident)
which(rownames(data)!=names(dge@ident))

j=1
label="PC1-2"
data=data_bg=dge@pca.rot[,1:2]
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
scale_colour_manual(values=myBrewerPalette,limits=levels(data$ident),drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}


j=4
label="tSNE"
data=data_bg=dge@tsne.rot
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
jpeg(file=paste0("OrigSet_12clusters_",label,"_bg.jpeg"),res=300,height=2500,width=4300)
MultiPlotList(plotset,cols = 5)
dev.off()
pdf(file=paste0("OrigSet_12clusters_",label,"_bg.pdf"),height=12,width=20)
MultiPlotList(plotset,cols = 5)
dev.off()

### end
### Determined Human1-3 are highly consistent; Human6 had steroid drug usage; only Human5 had Leydig cells; Human4 had batch effect from Human5; decided to remove Human4 and Human6 (renamed Human5 as Human4 in the manuscript)
