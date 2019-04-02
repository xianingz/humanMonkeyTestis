### Project somatic cells of the other 3 subjects to PCA of Human5
# 1.9.2019 by Qianyi
### Human5 had the largest amount of somatic cells and most comprehensive representation of somatic cell types

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


# 1.9.2019 by Qianyi 

# project the other 3 Human subjects to PCA of Human5 NoC5 Somatic Subset

### load object for directly merged somatic cells of 4 human subjects after removing Cluster 5
load(file="HumanMerged4-1235-NoC5Somatic.Robj")
dge4=dge
dge4 # 35063 genes across 3765 samples.

### load Existing PCA for Human5 (1 subject) NoCluster5 Somatic subset
load(file="Human5-NoC5Somatic.Robj")
dge5=dge
dge5 # 31253 genes across 2278 samples

### Extract gene expression data for the somatic cells of 3 human subjects (removed C5 doublets)
cells.use=names(dge4@ident)[which(dge4@meta.data$indiv !="Human5")]
table(dge4@ident,dge4@meta.data$indiv)
table(dge4@ident[cells.use])
     Human1 Human2 Human3 Human5
  1       2     10     11     45
  2      26     84     45     86
  3     116     44    145      6
  4      16     38     14    140
  5      13     36     24    231
  6       1      1      3    394
  7       1      0     10    643
  8       0      0      5    724
  9      29    299    346      5
  10      8     86     36      0
  11      5     24      9      4
  1   2   3   4   5   6   7   8   9  10  11 
 23 155 305  68  73   5  11   5 674 130  38 

# extract genes used in existing PCA for somatic arm of Human subject 5
# extract cells in somatic subset of Human subjects 1-2-3
dgedata=dge4@raw.data[rownames(dge5@data),cells.use]
dim(dgedata) # [1] 31253  1487
dge <- CreateSeuratObject(raw.data = dgedata, project="Human123-NoC5Somatic_Human5Genes", min.cells=0, min.genes=1)
dge <- NormalizeData(object = dge)
dge <- ScaleData(object = dge)
dge3=dge 
dge3 # 31253 genes across 1487 samples.                

###### project to somatic cells of Human5 (1 subject)
### use the same set of genes as used in old PCA
dgeold=dge5
dgenew=dge3
dgefile="ProjectHuman123toPCAofHuman5/predict_"


norm=dgenew@data
scale=dgenew@scale.data
data.use=norm;genes.use=rownames(dgeold@data);do.center=T; do.scale=T;scale.max=10
# use the same center and scale as the old PCA for 6 ST
            scale.data <- matrix(NA, nrow = length(genes.use), ncol = ncol(data.use))
            dimnames(scale.data) <- dimnames(data.use)
              bin.size <- 1000
              max.bin <- floor(length(genes.use)/bin.size) + 1
              cat("Scaling data matrix", file = stderr())
              pb <- txtProgressBar(min = 0, max = max.bin, style = 3, file = stderr())
              for(i in 1:max.bin) {
                my.inds <- ((bin.size * (i - 1)):(bin.size * i - 1))+1
                my.inds <- my.inds[my.inds <= length(genes.use)]
                #print(my.inds)
                new=as.matrix(data.use[genes.use[my.inds], ])
                old=as.matrix(dgeold@data[genes.use[my.inds], ])
                new.data <- (new-apply(old,1,mean))/apply(old,1,sd)
                new.data[new.data>scale.max] <- scale.max
                scale.data[genes.use[my.inds], ] <- new.data
                setTxtProgressBar(pb, i)
              }
              close(pb)
            
scalebyold=scale.data

# Using predict() function from stats package
oldpca=prcomp(t(dge5@scale.data[dge5@var.genes,]),center=F,scale=F)
oldpca1=oldpca
oldpca1$x[1:5,1:3]
                           PC1      PC2       PC3
Human5_TCTCAATTTCGT  -4.908917 20.83404 -2.017548
Human5_CGGAGGGTCAAC  -5.089969 17.18946 -2.028244
Human5_TGTTTGAAGATC  -4.632483 16.01297 -1.887121
Human5_AATCACTTATTC   2.850502 17.61249  2.082322
Human5_CTAATTCCTACA -18.103534 12.54379  7.701427
dgeold@dr$pca@cell.embeddings[1:5,1:3]
                          PC1       PC2       PC3
Human5_TCTCAATTTCGT  4.908917 -20.83404 -2.017548
Human5_CGGAGGGTCAAC  5.089969 -17.18946 -2.028244
Human5_TGTTTGAAGATC  4.632483 -16.01297 -1.887121
Human5_AATCACTTATTC -2.850502 -17.61249  2.082322
Human5_CTAATTCCTACA 18.103534 -12.54379  7.701427
# note: these are the same, except some are flipped 


# scaled and centered by old only
oldpca=oldpca1
tmp2=predict(oldpca,newdata=t(scalebyold[dgeold@var.genes,]))

oldpca$rotation=dgeold@dr$pca@gene.loadings
oldpca$x[,1:40]=dgeold@dr$pca@cell.embeddings
tmp1=predict(oldpca,newdata=t(scalebyold[dgeold@var.genes,]))

tmp1[1:5,1:3]
                             PC1       PC2      PC3
Human1.1_TGAATTTGCCTT 23.0545848 -7.764598 8.661360
Human1.1_AGCCATGCTTCC -0.1710823 -4.555276 8.858004
Human1.1_TTCCTCAGCCTC  0.7313332 -5.006140 7.960026
Human1.1_TGACCTTCCGCT -1.3719078 -2.168412 9.433789
Human1.1_TCATCTTCTGCT  0.6821735 -3.707660 8.767539
tmp2[1:5,1:3]
                              PC1      PC2      PC3
Human1.1_TGAATTTGCCTT -23.0545848 7.764598 8.661360
Human1.1_AGCCATGCTTCC   0.1710823 4.555276 8.858004
Human1.1_TTCCTCAGCCTC  -0.7313332 5.006140 7.960026
Human1.1_TGACCTTCCGCT   1.3719078 2.168412 9.433789
Human1.1_TCATCTTCTGCT  -0.6821735 3.707660 8.767539
# again, just flipped for some PCs

### save the projection result
  pca.obj <- new(
    Class = "dim.reduction",
    gene.loadings = matrix(),
    cell.embeddings = tmp1,
    sdev = numeric(),
    key = "PC"
  )
dgenew@dr$pca<- pca.obj
dge=dgenew
dge3=dge
save(dge,file="Human123ProjectedtoHuman5Somatic.Robj")

### add ident
dge=dgenew

dge@meta.data$indiv=gsub("\\..*","",dge@meta.data$orig.ident)
dge=SetAllIdent(dge,id="indiv")
table(dge@ident) 

ident=dge4@ident
dge=AddMetaData(dge,ident,"Allres.0.6order")
dge=SetAllIdent(dge,id="Allres.0.6order")
table(dge@ident) 

ident=pbmc@ident
dge=AddMetaData(dge,ident,"AllCCA")
dge=SetAllIdent(dge,id="AllCCA")
table(dge@ident) 

dgenew=dge
dge3=dge
save(dge,file="Human123ProjectedtoHuman5Somatic.Robj")


### 2D density plot
library(gplots)


setsall="Somatic"
### somatic cells
j=1
set=setsall[j]
dgename=dgefile="ProjectHuman123toPCAofHuman5/project_Somatic_"

###### plot PCA
### plot for new dataset
dge=dgenew

dge=SetAllIdent(dge,id="orig.ident")
pdf("predict_PCA_123-NoC5Somatic_Rep.pdf",width=10.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:9,2:1,4:3,10)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

dge=SetAllIdent(dge,"indiv")
pdf("predict_PCA_123-NoC5Somatic_Subject.pdf",width=10,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

dge=SetAllIdent(dge,id="Allres.0.6order")
pdf("predict_PCA_123-NoC5Somatic_AllDirectlyMergedres.0.6order.pdf",width=9.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

dge=SetAllIdent(dge,id="AllCCA")
pdf("predict_PCA_123-NoC5Somatic_AllCCA.pdf",width=9.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

###### Visualize individual batch and subject
### plot PCs and tSNE for each batch using the other batches as background
library(scales)
xlims=ylims=list()
pos=c("bottomright","bottomright","bottomright","bottomleft")
# dge4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-8,25),c(-8,25),c(-42,18),c(-38,46))
ylims[[1]]=list(c(-22,10),c(-12.5,12.5),c(-28,9),c(-48,46))

dge=dgenew
dim=list(dge@dr$pca@cell.embeddings[,1:2],dge@dr$pca@cell.embeddings[,c(1,3)])
dimold=list(dgeold@dr$pca@cell.embeddings[,1:2],dgeold@dr$pca@cell.embeddings[,c(1,3)])
xlim=xlims[[i]]
ylim=ylims[[i]]

### plot PCs for each subject using the other subjects as background
dge=SetAllIdent(dge,id="indiv")
sets=subject
plot2set=plot3set=plot4set=plottset=NULL
cols=myBrewerPalette1

pdf("OrigSubject_PCtSNE.pdf",height=4.6,width=4.6)
par(mfrow=c(2,2),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:2){
for(seti in 1:3){
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

seti=4
set="Human5"
  ident=as.factor(rep(set,length(dgeold@ident)))
  names(ident)=names(dgeold@ident)
  tmp=dimold[[j]][names(ident),]
  plot(tmp,pch=16,cex=0.4,col=alpha(cols[seti],0.8)[ident],xlim=xlim[[j]],ylim=ylim[[j]])
legend(pos[j],pch=16,set,col=cols[seti])
}

dev.off()


### plot PCs for each subject using the other subjects as background
### Visualize the clusters instead of subjects
dgenew=SetAllIdent(dgenew,id="AllCCA")
dgeold=SetAllIdent(dgeold,id="AllCCA")
sets=subject
cols=myBrewerPalette
source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/Rcode_multiplot.R")
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())


pdf(file=paste0("OrigSet_AllCCA7clusters_bg.pdf"),height=6,width=12)

j=1
plotset=NULL
dge=dgenew
data=data_bg=as.data.frame(dim[[j]])
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(PC1 = median(PC1), PC2 = median(PC2)) -> centers
label.size=6
for(i in 1:3){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(dge@ident)[which(dge@meta.data$indiv==set)]
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
i=4
set=sets[i]
dge=dgeold
data=as.data.frame(dimold[[j]])
data_bg=as.data.frame(dim[[j]])
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(PC1 = median(PC1), PC2 = median(PC2)) -> centers
label.size=6
cells.use=names(dge@ident)[which(dge@meta.data$indiv==set)]
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


j=2
plotset2=NULL
dge=dgenew
data=data_bg=as.data.frame(dim[[j]])
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(PC1 = median(PC1), PC3 = median(PC3)) -> centers
label.size=6
for(i in 1:3){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(dge@ident)[which(dge@meta.data$indiv==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=PC1,y=PC3,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=cols,limits=levels(data$ident),drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset2[[i]]=p3
}
i=4
set=sets[i]
dge=dgeold
data=as.data.frame(dimold[[j]])
data_bg=as.data.frame(dim[[j]])
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(PC1 = median(PC1), PC3 = median(PC3)) -> centers
label.size=6
cells.use=names(dge@ident)[which(dge@meta.data$indiv==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=PC1,y=PC3,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=cols,limits=levels(data$ident),drop=FALSE) +
coord_cartesian(ylim=ylims[[1]][[j]],xlim=xlims[[1]][[j]]) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset2[[i]]=p3


multiplot(c(plotset,plotset2),cols = 4)
dev.off()

### end




###### Heatmap ranges for PCA
numPCs=c(1,2)
range(dgenew@dr$pca@cell.embeddings[,numPCs[1]])
range(dgeold@dr$pca@cell.embeddings[,numPCs[1]])
range(dgenew@dr$pca@cell.embeddings[,numPCs[2]])
range(dgeold@dr$pca@cell.embeddings[,numPCs[2]])
xyrange=c(-8,25,-22,10)

numPCs=c(1,3)
range(dgenew@dr$pca@cell.embeddings[,numPCs[2]])
range(dgeold@dr$pca@cell.embeddings[,numPCs[2]])
xyrange=c(-8,25,-12.5,12.5)


numPCs=c(1,5)
range(dgenew@dr$pca@cell.embeddings[,numPCs[2]])
range(dgeold@dr$pca@cell.embeddings[,numPCs[2]])
xyrange=c(-60,45,-45,50)
xyrange=c(-72,72,-120,18)

numPCs=c(1,6)
range(dgenew@dr$pca@cell.embeddings[,numPCs[2]])
range(dgeold@dr$pca@cell.embeddings[,numPCs[2]])
xyrange=c(-60,45,-32,90)
xyrange=c(-70,72,-120,16)

#setwd("C:/Users/qzm/Desktop/DropSeq/plot/")
#dgePC16s=read.table("figJun2017_MouseAdultoldSca1/oldsertoli10_PC1-6.txt")
#dgePC12s=read.table("figJun2017_MouseAdultoldSca1/oldSca1_PC1-6.txt")
#dgePC11i=read.table("figJun2017_MouseAdultoldSca1/oldinterstitial5_PC1-6.txt")

##### old
h6PCs=dgeold@dr$pca@cell.embeddings[,numPCs]

jpeg(file=paste0(dgefile,"_hist2d_24ST_samerange_PC",paste(numPCs,collapse="PC"),"raw.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
h6=hist2d_range(h6PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=set,cex.main=2)
dev.off()
# log
h6log=log(h6$counts+1)
jpeg(file=paste0(dgefile,"_hist2d_samerange_PC",paste(numPCs,collapse="PC"),"log.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h6log,col=redblue100,main=set,axes=F,an=F,cex.main=2)
dev.off()

write.table(h6$counts,paste0(dgefile,"hist2d_50by50_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(h6log,paste0(dgefile,"hist2d_50by50_24ST",set,"_samerange_PC",paste(numPCs,collapse="PC"),"log.txt"),row.names=T,col.names=T,quote=F,sep="\t")


##### for all batches in each experiment
### new
h2PCs=dgenew@dr$pca@cell.embeddings[,numPCs]
set=setsall[j]

### plot together
jpeg(file=paste0(dgefile,"hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"all.jpeg"),res=300,height=1900,width=650)
par(mfrow=c(3,1),mar=c(0,1,1,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=set,cex.main=2)
# log
h2log=log(h2$counts+1)
image(h2log,col=redblue100,main=set,axes=F,an=F,cex.main=2)
# ratio with symmetric color scale
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=set,axes=F,an=F,cex.main=2)
dev.off()

### plot separately
jpeg(file=paste0(dgefile,"_hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=set,cex.main=2)
dev.off()
# log
h2log=log(h2$counts+1)
jpeg(file=paste0(dgefile,"_hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"log.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2log,col=redblue100,main=set,axes=F,an=F,cex.main=2)
dev.off()
# ratio
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
jpeg(file=paste0(dgefile,"_hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"by6STratio.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=set,axes=F,an=F,cex.main=2)
dev.off()

### write table
write.table(h2$counts,paste0(dgefile,"hist2d_50by50_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(h2by6STratio,paste0(dgefile,"hist2d_50by50_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"logratio.txt"),row.names=T,col.names=T,quote=F,sep="\t")


### plot heatmap for each individual subject: 
for(set in subject[1:3]){
h2PCs=dgenew@dr$pca@cell.embeddings[which(grepl(set,dgenew@meta.data$indiv)),numPCs]
h6PCs=dgeold@dr$pca@cell.embeddings[,numPCs] #[which(!grepl(set, rownames(dgeold@dr$pca@cell.embeddings))),numPCs]
h6=hist2d_range(h6PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=set,cex.main=2)

### plot together
jpeg(file=paste0(dgefile,"hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"all.jpeg"),res=300,height=1900,width=650)
par(mfrow=c(3,1),mar=c(0,1,1,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=set,cex.main=2)
# log
h2log=log(h2$counts+1)
image(h2log,col=redblue100,main=set,axes=F,an=F,cex.main=2)
# ratio with symmetric color scale
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=set,axes=F,an=F,cex.main=2)
dev.off()
}

### plot separately
jpeg(file=paste0(dgefile,"_hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=set,cex.main=2)
dev.off()
# log
h2log=log(h2$counts+1)
jpeg(file=paste0(dgefile,"_hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"log.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2log,col=redblue100,main=set,axes=F,an=F,cex.main=2)
dev.off()
# ratio
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
jpeg(file=paste0(dgefile,"_hist2d_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"by6STratio.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=set,axes=F,an=F,cex.main=2)
dev.off()

### write table
write.table(h2$counts,paste0(dgefile,"hist2d_50by50_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(h2by6STratio,paste0(dgefile,"hist2d_50by50_",set,"_samerange_PC",paste(numPCs,collapse="PC"),"logratio.txt"),row.names=T,col.names=T,quote=F,sep="\t")


##### for individual batch
for(i in 1:length(sets)){
seti=sets[i]
h2PCs=dge@dr$pca@cell.embeddings[which(gsub("_.*","",rownames(dge@dr$pca@cell.embeddings))==seti),numPCs]
### plot together
jpeg(file=paste0(dgefile,"hist2d_",set,"_",seti,"_samerange_PC",paste(numPCs,collapse="PC"),"all.jpeg"),res=300,height=3600,width=1100)
par(mfrow=c(3,1),mar=c(0,1,2,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(h2PCs,nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=seti,cex.main=2)
# log
h2log=log(h2$counts+1)
image(h2log,col=redblue100,main=seti,axes=F,an=F,cex.main=2)
# ratio
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=seti,axes=F,an=F,cex.main=2)
dev.off()

### plot separately
jpeg(file=paste0(dgefile,"_hist2d_",set,"_",seti,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
h2=hist2d_range(dge@dr$pca@cell.embeddings[which(gsub("_.*","",rownames(dge@dr$pca@cell.embeddings))==seti),1:2],nbin=50,range=xyrange,same.scale=F,col=redblue100,axes=F,main=seti,cex.main=2)
dev.off()
# log
h2log=log(h2$counts+1)
jpeg(file=paste0(dgefile,"_hist2d_",set,"_",seti,"_samerange_PC",paste(numPCs,collapse="PC"),"log.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2log,col=redblue100,main=seti,axes=F,an=F,cex.main=2)
dev.off()
# ratio
h2by6STratio=(log(h2$counts+1)/sum(log(h2$counts+1)))-(log(h6$counts+1)/sum(log(h6$counts+1)))
z=max(abs(h2by6STratio))
jpeg(file=paste0(dgefile,"_hist2d_",set,"_",seti,"_samerange_PC",paste(numPCs,collapse="PC"),"by6STratio.jpeg"),res=300,height=1200,width=1100)
par(mar=c(1,1,2,1),mgp=c(1.5,0.5,0))
image(h2by6STratio,zlim=c(-z,z),col=redblue100,main=seti,axes=F,an=F,cex.main=2)
dev.off()

### write table
write.table(h2$counts,paste0(dgefile,"hist2d_50by50_",set,"_",seti,"_samerange_PC",paste(numPCs,collapse="PC"),"raw.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(h2by6STratio,paste0(dgefile,"hist2d_50by50_",set,"_",seti,"_samerange_PC",paste(numPCs,collapse="PC"),"logratio.txt"),row.names=T,col.names=T,quote=F,sep="\t")
}

