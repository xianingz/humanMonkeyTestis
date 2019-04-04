### project the other 4 subjects to PCA of Monkey3
# 1.30.2019 by Qianyi


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
datainfo=read.table("datainfo_monkey")
datainfo
dataset=as.character(datainfo[,2])
subject=unique(gsub("-.*","",dataset))
length(dataset) # 9
length(subject) # 5


### load object for Monkey3 somatic 
load(file="Monkey3Somatic.Robj")
dge3=dge
dgeold=dge
dge # 35063 genes across 3765 samples.


### load object for CCA of 4 subjects after removing Cluster 5
load(file="Monkey4Somatic_CCA.Robj")
pbmc


# 1.30.2019 by Qianyi 

# project the other 3 Monkey subjects to PCA of Monkey3 Somatic Subset

### load object for Directly-merged 5 Monkey subjects
##### Ensembl ID converted to human gene symbols with 1-1 orthologues 
load(file="MonkeyMerged5Somatic_GeneNames.Robj")
dge5=dge # 18370 genes across 2562 samples.
dge 
##### labeld with cell types
load(file="MonkeyMerged5Somatic-NoC1.Robj")
dge5celltype=dge # 18370 genes across 2562 samples.
dge 

### load Existing PCA for Monkey3 (1 subject) NoCluster5 Somatic subset
load(file="Monkey3Somatic.Robj")
dge3=dge
dge3 # 17013 genes across 1528 samples.

### Extract gene expression data for the somatic cells of 3 monkey subjects (removed C5 doublets)
cells.use=names(dge5@ident)[which(dge5@meta.data$indiv !="Monkey3")]
table(dge5@ident,dge5@meta.data$indiv)
table(dge5@ident[cells.use])
    Monkey1 Monkey2 Monkey3 Monkey4 Monkey5
  1      19      21      83      62      22
  2      39      29     311      63      92
  3      39      62     625      10      19
  4      80      52     409      51      12
  5     154       6      12      11       3
  6      15       4       7      20       2
  7      59       6      62      50      14
  8       4       2      19       6       6
  1   2   3   4   5   6   7   8 
124 223 130 195 174  41 129  18  

# extract genes used in existing PCA for somatic arm of Monkey subject 3
# extract cells in somatic subset of Monkey subjects 1-2-4-5
dgedata=dge5@raw.data[rownames(dge3@data),cells.use]
dim(dgedata) # [1] 17013  1034
dge <- CreateSeuratObject(raw.data = dgedata, project="Monkey1245Somatic_Monkey3Genes", min.cells=0, min.genes=1)
dge <- NormalizeData(object = dge)
dge <- ScaleData(object = dge)
dge4=dge 
dge4 # 17013 genes across 1034 samples.               

###### project to somatic cells of Monkey3 (1 subject)
### use the same set of genes as used in old PCA
dgeold=dge3
dgenew=dge4
dgefile="ProjectMonkey1245toPCAofMonkey3Somatic/predict_"


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

#(1) Using predict() function from stats package
oldpca=prcomp(t(dge3@scale.data[dge3@var.genes,]),center=F,scale=F)
oldpca1=oldpca
oldpca1$x[1:2,1:3]
dgeold@dr$pca@cell.embeddings[1:2,1:3]
                              PC1        PC2        PC3
Monkey3.1_CCTGATGAATGA  -1.470474 -4.1439950 -1.7580303
Monkey3.1_ATAGCCTCACTA -47.662684  0.9271215 -0.2786888
                              PC1        PC2       PC3
Monkey3.1_CCTGATGAATGA  -1.470474 -4.1439950 1.7580303
Monkey3.1_ATAGCCTCACTA -47.662684  0.9271215 0.2786888
# note: these PCs are the same, except some are flipped 


# scaled and centered by old only
oldpca=oldpca1
tmp2=predict(oldpca,newdata=t(scalebyold[dgeold@var.genes,]))

oldpca$rotation=dgeold@dr$pca@gene.loadings
oldpca$x[,1:40]=dgeold@dr$pca@cell.embeddings
tmp1=predict(oldpca,newdata=t(scalebyold[dgeold@var.genes,]))

tmp1[1:5,1:3]
tmp2[1:5,1:3]
                              PC1         PC2       PC3
Monkey1.1_ATTGAACATCTG -42.516001 -0.86679489 4.5592978
Monkey1.1_AGCCAGCCTTCA   2.111140 -4.45914887 0.6539861
Monkey1.1_ACCTAATGGTAC -41.949557 -0.68053403 3.8857985
Monkey1.1_CATACTTGGCTT   2.919509  0.05221273 1.9322264
Monkey1.1_TTCCCCAAACTT   3.684806  2.02610537 0.2331531
                              PC1         PC2        PC3
Monkey1.1_ATTGAACATCTG -42.516001 -0.86679489 -4.5592978
Monkey1.1_AGCCAGCCTTCA   2.111140 -4.45914887 -0.6539861
Monkey1.1_ACCTAATGGTAC -41.949557 -0.68053403 -3.8857985
Monkey1.1_CATACTTGGCTT   2.919509  0.05221273 -1.9322264
Monkey1.1_TTCCCCAAACTT   3.684806  2.02610537 -0.2331531
# again, same values, just flipped for some PCs

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
dge4=dge
save(dge,file="Monkey1245ProjectedtoMonkey3Somatic.Robj")

### add ident
dge=dgenew

dge@meta.data$indiv=gsub("\\..*","",dge@meta.data$orig.ident)
dge=SetAllIdent(dge,id="indiv")
table(dge@ident) 

ident=dge5@ident
dge=AddMetaData(dge,ident,"Allres.0.5order")
dge=SetAllIdent(dge,id="Allres.0.5order")
table(dge@ident) 

ident=dge5celltype@ident
dge=AddMetaData(dge,ident,"Allcelltype")
dge@meta.data$Allcelltype=factor(dge@meta.data$Allcelltype,levels=levels(dge5celltype@ident))
dge=SetAllIdent(dge,id="Allcelltype")
dge@ident=factor(dge@ident,levels=levels(dge5celltype@ident))
table(dge@ident) 

ident=pbmc@ident
dge=AddMetaData(dge,ident,"AllCCA")
dge=SetAllIdent(dge,id="AllCCA")
table(dge@ident) 

dgenew=dge
dge4=dge
save(dge,file="Monkey1245ProjectedtoMonkey3Somatic.Robj")


setsall="Somatic"
### somatic cells
j=1
set=setsall[j]
dgename=dgefile="ProjectMonkey1245toPCAofMonkey3Somatic/project_Somatic_"

###### plot PCA
### plot for new dataset
dge=dgenew

dge=SetAllIdent(dge,id="orig.ident")
pdf("predict_PCA_1245Somatic_Rep.pdf",width=10.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:6,9,2:1,4:3)])
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:6,9,2:1,4:3)])
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:6,9,2:1,4:3)])
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette[c(5:6,9,2:1,4:3)])
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

dge=SetAllIdent(dge,"indiv")
pdf("predict_PCA_1245Somatic_Subject.pdf",width=10,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

dge=SetAllIdent(dge,id="Allres.0.5order")
pdf("predict_PCA_1245Somatic_AllDirectlyMergedres.0.5order.pdf",width=9.5,height=8)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=T,cols.use=myBrewerPalette)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = .8,do.label=F,cols.use=myBrewerPalette)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

dge=SetAllIdent(dge,id="AllCCA")
pdf("predict_PCA_1245Somatic_AllCCA.pdf",width=9.5,height=8)
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
pos=c("bottomleft","topleft","bottomright","bottomleft")
# dge4 used top 9 PCs for tSNE
i=1
xlims[[1]]=list(c(-49,5),c(-49,6),c(-48,6),c(-38,46))
ylims[[1]]=list(c(-33,5),c(-7,17),c(-28,9),c(-48,46))

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

pdf("OrigSubject_PCtSNE.pdf",height=4.6,width=6.9)
par(mfrow=c(2,3),mar=c(2.8,2.5,0.5,0.5),mgp=c(1.5, 0.5, 0))
for(j in 1:2){
for(seti in c(1,2,4,5)){
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

seti=3
set="Monkey3"
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

dgenew=SetAllIdent(dgenew,id="Allcelltype")
dgenew@ident=factor(dgenew@ident,levels=levels(dge5celltype@ident))
dgeold=SetAllIdent(dgeold,id="Allcelltype")
dgeold@ident=factor(dgeold@ident,levels=levels(dge5celltype@ident))

dgenew=SetAllIdent(dgenew,id="Allres.0.5order")
dgeold=SetAllIdent(dgeold,id="Allres.0.5order")

sets=subject
cols=myBrewerPalette
source("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/Rcode_multiplot.R")
gg.xax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.x = element_text(face=z, colour=y, size=x),
                                                              axis.text.x  = element_text(angle=90, vjust=0.5, size=x2)))
gg.yax=function(x=16,y="#990000",z="bold",x2=12) return(theme(axis.title.y = element_text(face=z, colour=y, size=x),
                                                              axis.text.y  = element_text(angle=90, vjust=0.5, size=x2)))
no.legend.title=theme(legend.title=element_blank())



j=1
plotset=NULL
dge=dgenew
data=data_bg=as.data.frame(dim[[j]])
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(PC1 = median(PC1), PC2 = median(PC2)) -> centers
label.size=6
for(i in c(1,2,4,5)){
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
i=3
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
for(i in c(1,2,4,5)){
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
i=3
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


pdf(file=paste0("OrigSet_AllCCA7clusters_bg.pdf"),height=6,width=12)

pdf(file=paste0("OrigSet_AllDirectlyMergedres.0.5order_bg.pdf"),height=6,width=9)

multiplot(plotset,cols = 3)
multiplot(plotset2,cols = 3)
dev.off()

### end


