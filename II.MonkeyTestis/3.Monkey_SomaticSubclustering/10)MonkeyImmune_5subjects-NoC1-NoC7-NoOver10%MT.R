### Monkey Immune cells subclustering
# 3.25.2019 by Qianyi
### Extracted somatic cells and directly merged, did somatic subclustering
### Removed Cluster 1/8 of SCyte/STids-Somatic doublets from somatic subclustering 
### Corrected for batch effect by CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing doublets
### Removed Cluster 7/7 of SCyte-Somatic doublets from somatic subclustering by CCA
### Re-do CCA and did Subclustering for Somatic-NoC1-NoC7 cells after removing C1 and C7 doublets
### Extract Macrophage/Tcell cluster, directly-merged and do subclustering
#-> Decided to use Macrophage/Tcell clusters from 3-species somatic clustering, as that is more accurate
### No cell with >=10% MT genes in the Macrophage/Tcell immune subset 

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

dgefile=dgename="CCA5Subjects-Somatic-NoC1-NoC7/"
### load object for Immune subclustering
load(file="MonkeyImmune_DirectMerged5Subjects.Robj")
dge # 6713 genes across 37 samples.
### load object for 3-species somatic clustering 
load(file=paste0("../Somatic_3Species_11Subjects_CCA.Robj"))
testisall=testis
testis # 10765 genes across 11158 samples.



### load object of CCA and somatic subclustering after removing cells with >=10% MT genes, Cluster 1/8 and 7/7 of spermatid-myoid doublets for somatic subclustering 
dgefile=dgename="CCA5Subjects-Somatic-NoC1-NoC7/"
load(file="Monkey5Somatic-NoC1-NoC7_CCA-No10pctMT.Robj")
testis # 17251 genes across 2098 samples.
testis14=testis

table(testis14@ident)
Macrophage/Tcell      Endothelial       m-Pericyte       f-Pericyte 
              37              178              117               38 
           Myoid        ImmLeydig       DiffLeydig 
             362             1324               42
testis=SubsetData(testis14,ident.use="Macrophage/Tcell")
table(testis@ident)
#Macrophage/Tcell 
#              37

### Comparison with separeted Macrophage and T cell for Monkey from CCA of 3-species Somatic
load(file=paste0("../Somatic_3Species_11Subjects_CCA.Robj"))
testisall=testis
testis # 10765 genes across 11158 samples.

table(testisall@ident[names(dge@ident)],dge@ident)
#      1  2
#  9  17  0
#  10  7 13


pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()


## label cell types for Immune subclustering
### 1. label cell types from clustering of 3-species somatic merged by CCA
# note: decided to keep this cluster solution as final separation for monkey Tcell/macrophag
table(testisall@ident[names(dge@ident)])
# 9 10  
#17 20  
ident3=as.character(testisall@ident[names(dge@ident)])
names(ident3)=names(dge@ident)
ident3[which(ident3==9)]<-"Macrophage"
ident3[which(ident3==10)]<-"Tcell"
dge=AddMetaData(dge,ident3,"CellType1")
dge=SetAllIdent(dge,id="CellType1") # used this
save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")


pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune1.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Immune1_3Species_markersall_mindiff0.2_logfc2fold_3.11.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
#Macrophage      Tcell 
#        17         20 
#        41         13 

### 2. label cell types from subclustering by Louvain-Jaccard clustering
table(dge@meta.data$res.1order)
# 1  2 
#24 13
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident2[which(ident2==1)]<-"Macrophage"
ident2[which(ident2==2)]<-"Tcell"
dge=AddMetaData(dge,ident2,"CellType2")
dge=SetAllIdent(dge,id="CellType2")
save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")

pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune2.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Immune2_res.1order_markersall_mindiff0.2_logfc2fold_3.11.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
# 1  2 
#24 13
#32 32 

### 3. label cell types from k-means of UMAP
kclust=kmeans(dge@dr$umap@cell.embeddings,centers=2)
table(kclust$cluster)
# 1  2 
#19 18 

ident3=kclust$cluster
ident3[which(ident3==2)]<-"Macrophage"
ident3[which(ident3==1)]<-"Tcell"
dge=AddMetaData(dge,ident3,"CellType3")
dge=SetAllIdent(dge,id="CellType3")
save(dge,file="MonkeyImmune_DirectMerged5Subjects.Robj")

pdf(paste0(dgefile,"PCA_tSNE_UMAP_Immune3.pdf"),width=9.5,height=7)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = .8,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = .8,do.label=F)
plot4=DimPlot(dge,reduction.use="umap",do.return = TRUE,pt.size = .8,do.label=T)
plot5=TSNEPlot(dge,do.return=TRUE,pt.size = .8,do.label=T)
plot_grid(plot2,plot3,plot4,plot5,ncol = 2)
dev.off()

markers=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
write.table(markers,"Immune3_kmeans2_markersall_mindiff0.2_logfc2fold_3.11.2019.txt",col.names=T,row.names=T,quote=F,sep="\t")
table(dge@ident)
table(markers$cluster)
#Macrophage      Tcell 
#        18         19 
#        42         11 



### Add these two clusters back to All Somatic cells of Monkey
dge=SetAllIdent(dge,id="CellType1") # used this
testis=testis14
testis=SetAllIdent(testis,id="CellType2")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType2))
ident=as.character(testis@ident)
names(ident)=names(testis@ident)
ident2=as.character(dge@ident)
names(ident2)=names(dge@ident)
ident[names(ident2)]<-ident2
levels=c(levels(dge@ident)[2:1],levels(testis@ident)[2:7])
table(ident)[levels]
#      Tcell  Macrophage Endothelial  m-Pericyte  f-Pericyte       Myoid 
#         20          17         178         117          38         362 
#  ImmLeydig  DiffLeydig 
#       1324          42 

ident=factor(ident,levels=levels)
testis=AddMetaData(testis,ident,"CellType3")

testis=SetAllIdent(testis,id="CellType3")
testis@ident=factor(testis@ident,levels=levels(testis@meta.data$CellType3))
testis14=testis
save(testis,file=paste0("Monkey5Somatic-NoC1-NoC7_CCA-No10pctMT.Robj"))

table(testis1@ident[names(testis@ident)],testis@ident)
              Tcell Macrophage Endothelial m-Pericyte f-Pericyte Myoid
  Tcell          20          0           0          0          0     0
  Macrophage      0         17           0          0          0     0
  Endothelial     0          0         162          0          0     0
  Pericyte        0          0           0        116          0     0
  Subcluster2     0          0           3          1         32     0
  Myoid           0          0           0          0          0   355
  ImmLeydig       0          0          13          0          6     7
  DiffLeydig      0          0           0          0          0     0
             
              ImmLeydig DiffLeydig
  Tcell               0          0
  Macrophage          0          0
  Endothelial         2          0
  Pericyte            1          0
  Subcluster2         6          0
  Myoid             113          3
  ImmLeydig        1201          1
  DiffLeydig          1         38
