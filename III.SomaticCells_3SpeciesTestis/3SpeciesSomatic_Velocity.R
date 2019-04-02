### Estimating RNA Velocity for Somatic Cells for Each of 3 species
# 4.1.2019 by Qianyi

path <- "/scratch/junzli_flux/qzm/Dropseq_analysis/HumanMonkey/Somatic/"
setwd(path)

###### Data loading
library(velocyto.R)
library(h5)

### load loom matrix
species=c("Human","Monkey","Mouse")

set="Somatic"
ss=species[1]
ss=species[2]
ss=species[3]

for(ss in species){

ldat <- read.loom.matrices(paste0("merged",ss,set,".loom"))
emat <- ldat$spliced; nmat <- ldat$unspliced
#colnames(emat)=gsub("-",".",colnames(emat))
colnames(emat)=gsub(":","_",colnames(emat))
colnames(emat)=gsub("x","",colnames(emat))
#colnames(nmat)=gsub("-",".",colnames(nmat))
colnames(nmat)=gsub(":","_",colnames(nmat))
colnames(nmat)=gsub("x","",colnames(nmat))
which(colnames(emat)=="SPG2_TAAGACGTGCCA-TAAGACGTGCCA")
colnames(emat)[which(colnames(emat)=="SPG2_TAAGACGTGCCA-TAAGACGTGCCA")]<-"SPG2_TAAGACGTGCCA"
colnames(nmat)[which(colnames(nmat)=="SPG2_TAAGACGTGCCA-TAAGACGTGCCA")]<-"SPG2_TAAGACGTGCCA"
emat=emat[!duplicated(rownames(emat)),]
nmat=nmat[!duplicated(rownames(nmat)),]

# emat <- emat[,rownames(r$counts)]; nmat <- nmat[,rownames(r$counts)]; # restrict to cells that passed p2 filter
ematall=emat;nmatall=nmat

### double-check the number of cells from each batch
print(ss)
print(dim(emat))
print(table(gsub("_.*","",colnames(emat))))

[1] "Human"
[1] 56638  3722
Human1.1 Human1.2 Human1.3 Human1.4 Human1.5 Human2.1 Human2.2 Human3.1
      38       62       42       38       31      327      273      216
Human3.2   Human5
     422     2273

[1] "Monkey"
[1] 27121  2337
Monkey1.1 Monkey1.2   Monkey2 Monkey3.1 Monkey3.2 Monkey4.1 Monkey4.2 Monkey5.1
      226       160       160       656       780       120        88        59
Monkey5.2
       88

[1] "Mouse"
[1] 45810  5081
INT1 INT2 INT3 INT4 INT5 INT6 SER1 SER2 SER3 SER4 SER5 SER6 SER7 SER8 SPG1 SPG2
 106  134  102  495  164 1453   24   33  118  218  518  531  198  411    3    1
SPG3  ST1  ST2  ST3  ST4  ST5  ST6  ST7  ST8
  91   48   69   71   91   73   55   32   42


###### Read in cell cluster assignment and tSNE embedding 
tsne=read.table(paste0("../../data_DGE/",ss,"_",set,"_tSNE.txt"),row.names=1,header=T)
ident=read.table(paste0("../../data_DGE/",ss,"_",set,"_celltypes.txt"),stringsAsFactors=F,row.names=1)
levels=read.table(paste0("../../data_DGE/",ss,"_",set,"_celltype_order.txt"),stringsAsFactors=F)[,1]
cluster=factor(ident[,1],levels=levels)
names(cluster)=rownames(ident)
ncluster=length(levels(cluster))

library(RColorBrewer)
colH=brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3)]
colM=brewer.pal(8,"Set1")[c(7,6,4,5,1,2,3,8)]
colm=brewer.pal(8,"Set1")[c(7,6,4,5,3,8)]
cols=list(colH,colM,colm)
names(cols)=species

###### Prepare matrices and clustering data:
### select cells
cells.use=names(cluster)
### take cluster labels
cluster.label <- as.numeric(cluster)
names(cluster.label) <- names(cluster)
cell.colors <- cols[[ss]][as.numeric(cluster)]
names(cell.colors) <- names(cluster)
### take embedding
emat=ematall;nmat=nmatall
tsne <- as.matrix(tsne)[cells.use,]
if(ss!="Human"){
  pca=read.table(paste0("../../data_DGE/",ss,"_",set,"_PCA.txt"),row.names=1,header=T)
  emb <- as.matrix(pca)[cells.use,1:2]
  emb3 <- as.matrix(pca)[cells.use,c(1,3)]
  emb4 <- as.matrix(pca)[cells.use,c(1,4)]
  emb5 <- as.matrix(pca)[cells.use,c(1,5)]
  emb6 <- as.matrix(pca)[cells.use,c(1,6)]
  emb7 <- as.matrix(pca)[cells.use,c(1,7)]
} 
if(ss!="Mouse"){
  acc=read.table(paste0("../../data_DGE/",ss,"_",set,"_ACC.txt"),row.names=1,header=T)
  acc2 <- as.matrix(acc)[cells.use,1:2]
  acc3 <- as.matrix(acc)[cells.use,c(1,3)]
  umap=read.table(paste0("../../data_DGE/",ss,"_",set,"_UMAP.txt"),row.names=1,header=T)
  umap <- as.matrix(umap)[cells.use,]
}
emat <- emat[,cells.use]
nmat <- nmat[,cells.use]
which(! (rownames(tsne)%in% names(cell.colors) )) #integer(0)
colnames(emat)[which(!(colnames(emat) %in% cells.use))] # character(0)
#[1] "SPG2_TAAGACGTGCCA-TAAGACGTGCCA"
ematall=emat;nmatall=nmat

###### Distribution of reads per cell and per gene
pdf(file=paste0(ss,"_",set,"_exp.pdf"),width=6,height=4)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
### Spliced expression magnitude distribution across genes:
hist(log10(rowSums(emat)+1),col='wheat',xlab='log10(#spliced UMIs/gene + 1)',main='#Spliced UMIs/gene')
### Spliced Cell size
hist(log10(colSums(emat)),col='wheat',xlab='log10(#spliced UMIs/cell)',main='Spliced Cell Size')
### Unspliced expression magnitude distribution across genes:
hist(log10(rowSums(nmat)+1),col='wheat',xlab='log10(#unspliced UMIs/gene + 1)',main='#Unspliced UMIs/gene')
### Unspliced Cell size
hist(log10(colSums(nmat)),col='wheat',xlab='log10(#unspliced UMIs/cell)',main='Unspliced Cell Size')
dev.off()


######### Velocity estimation
######### 1. Velocity estimation using Gene-Relative Model
###### Gene filtering
### filtering genes to leave those that exceed some pre-defined g to the average expression magnitude
### Filter genes based on the minimum average expresion magnitude (in at least one of the clusters), output total number of resulting valid genes:
enfilters=list(c(0.015,0.01),c(0.01,0.005),c(0.05,0.015))
names(enfilters)=species

enfilter=enfilters[[ss]]
pdf(file=paste0(ss,"_",set,"_genefilter.pdf"),width=8,height=3)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
emat <- filter.genes.by.cluster.expression(ematall,cell.colors,min.max.cluster.average = enfilter[1],do.preview=TRUE)
nmat <- filter.genes.by.cluster.expression(nmatall,cluster.label,min.max.cluster.average = enfilter[2],do.preview=TRUE)
dev.off()
print(dim(emat)[1]) # 11025  11606  10451 
print(dim(nmat)[1]) # 11878  11995   10495
print(length(intersect(rownames(emat),rownames(nmat))))
# [1] 7922   9141    7517
# look at the resulting gene set
str(intersect(rownames(emat),rownames(nmat)))
# chr [1:3449] "Ppp1r42" "Cops5" "Arfgef1" "Tram1" "4930444P10Rik" "Rpl7" ...
filteredgenes=intersect(rownames(emat),rownames(nmat))

###### Estimate RNA velocity (using gene-relative model with k=20 cell kNN pooling and using top/bottom 2% quantiles for gamma fit):
velname="rvel_"
kCells=10 # putty
fit.quantile <- 0.02 # gamma fit is based on the top/bottom 2% of cells by spliced expression magnitude
rvel <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells=kCells,fit.quantile=fit.quantile)
# cell.dist=cell.dist,
#calculating cell knn ... done
#calculating convolved matrices ... done
#fitting gamma coefficients ... done. 
print(dim(rvel$current)) 
#[1] 3517 3722 
#[1] 3190 2098
#[1] 5020 5081

# Human, k=10
succesfful fit for 7922 genes
filtered out 3688 out of 7922 genes due to low nmat-emat correlation
filtered out 717 out of 4234 genes due to low nmat-emat slope
calculating RNA velocity shift ... done
calculating extrapolated cell state ... done

# Monkey, k=10
succesfful fit for 9141 genes
filtered out 5265 out of 9141 genes due to low nmat-emat correlation
filtered out 686 out of 3876 genes due to low nmat-emat slope
calculating RNA velocity shift ... <sparse>[ <logic> ] : .M.sub.i.logical() maybe inefficient
done

# Mouse, k=10
succesfful fit for 7517 genes
filtered out 1727 out of 7517 genes due to low nmat-emat correlation
filtered out 770 out of 5790 genes due to low nmat-emat slope
calculating RNA velocity shift ... done
calculating extrapolated cell state ... done


class(rvel)
#[1] "list"
names(rvel)
# [1] "cellKNN"        "conv.nmat.norm" "conv.emat.norm" "gamma"
# [5] "projected"      "current"        "deltaE"         "deltaT"
# [9] "ko"             "mult"           "cellKNN"        "kCells"
# current: observed current normalized exp state
# projected: extrapolated future state
# deltaE: unscaled transcriptional change
# ko / sfit: fit results
# gamma <- ko$g; offset <- ko$o

save(rvel,file=paste0(ss,"_",set,"_rvel_kCells",kCells,".Robj"))


###### Visualize velocity on the t-SNE embedding, using velocity vector fields:
velname="rvel_";kCells=10;load(file=paste0(ss,"_",set,"_rvel_kCells",kCells,".Robj"))

vel=rvel; alpha=0.4; cex=1;n=300
# need to set n>=100; otherwise the arrows point at all directions

arrow.scale=16;
pdf(file=paste0(ss,"_",set,"_",velname,'tsne.shift.plots.pdf'),height=5,width=4.5)
# when cell numbers are large, show.grid.flow=TRUE
  show.velocity.on.embedding.cor(tsne,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()

arrow.scale=8;
pdf(file=paste0(ss,"_",set,"_",velname,'UMAP.shift.plots.pdf'),height=5,width=4.5)
  show.velocity.on.embedding.cor(umap,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()

arrow.scale=4;
pdf(file=paste0(ss,"_",set,"_",velname,'ACC.shift.plots.pdf'),height=5,width=4.5)
  show.velocity.on.embedding.cor(acc2,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  show.velocity.on.embedding.cor(acc3,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()

if(ss=="Mouse" | ss=="Monkey"){
arrow.scale=8;
pdf(file=paste0(ss,"_",set,"_",velname,'PCA.shift.plots.pdf'),height=5,width=4.5)
  show.velocity.on.embedding.cor(emb,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  show.velocity.on.embedding.cor(emb3,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  show.velocity.on.embedding.cor(emb4,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  show.velocity.on.embedding.cor(emb5,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  show.velocity.on.embedding.cor(emb6,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  show.velocity.on.embedding.cor(emb7,vel,n=n,scale='sqrt',cell.colors=ac(cell.colors,alpha=alpha),cex=cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
dev.off()
}

}
}

### record log below
delta projections ... sqrt knn ... transition probs ... donele,show.gcalculating arrows ... grid estimates ... 
# Human
grid.sd= 1.254924  min.arrow.size= 0.02509849  max.grid.arrow.length= 0.0283761  done
grid.sd= 0.2564367  min.arrow.size= 0.005128735  max.grid.arrow.length= 0.0283761  done
grid.sd= 0.1137608  min.arrow.size= 0.002275217  max.grid.arrow.length= 0.0283761  done
grid.sd= 0.1174762  min.arrow.size= 0.002349524  max.grid.arrow.length= 0.0283761  done
# Monkey
grid.sd= 1.113016  min.arrow.size= 0.02226033  max.grid.arrow.length= 0.0283761  done
grid.sd= 0.2811162  min.arrow.size= 0.005622324  max.grid.arrow.length= 0.0283761  done
grid.sd= 0.131616  min.arrow.size= 0.002632321  max.grid.arrow.length= 0.0283761  done
grid.sd= 0.1247098  min.arrow.size= 0.002494196  max.grid.arrow.length= 0.0283761  done
# Mouse
grid.sd= 1.517896  min.arrow.size= 0.03035791  max.grid.arrow.length= 0.0283761  done



###### Visualize a fit for a particular gene (we reuse rvel to save on calcualtions here):
vel=rvel
# known somatic markers 
markers=c("VIM",
  "CD163","S100A4","TYROBP","LYZ","RGS1",
  "VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN",
  "NR2F2","MYH11","ACTA2","PTCH1","PTCH2",
  "DCN","TCF21","PDGFRA","PDGFRB","PDGFB",
  "DLK1","HSD17B3","STAR","CYP17A1","GATA4","NR5A1","EGR2","IGF1","IGFBP5","IGFBP3","INHBA","VIT",
  "THY1","ENG","VCAM1","ALCAM","CD44","NT5E","PDGFA","CASP3",
  "CLU","SOX9","AMH",
  "Id2","Il7r","Rora","Thy1","Ccl5","Cd52","Cd2","Adgre1","Dab2","Apoe","Mrc1",
  "Tie1","Vwf","Tek",
  "Nr2f2","Myh11","Acta2","Acta1","Pcna","Ptch1","Ptch2",
  "Dcn","Tcf21","Pdgfra","Pdgfrb","Pdgfb","Ly6a","Thra","Thrb","Arx",
  "Dlk1","Hsd17b3","Star","Cyp17a1","Gata4","Nr5a1","Egr2","Igf1","Igfbp5","Igfbp3","Inhba","Vit",
  "Thy1","Eng","Vcam1","Alcam","Cd44","Nt5e","Pdgfa","Casp3",
  "Sox9","Clu","Vim" )

markers=c("NR2F2","MYH11","ACTA2","PTCH1","PTCH2","DCN","TCF21","PDGFRA","PDGFRB","PDGFB","DLK1","HSD17B3","STAR","CYP17A1","GATA4","SF1","NR5A1","EGR2","IGF1","IGFBP5","IGFBP3","INHBA","VIT","THY1","ENG","VCAM1","ALCAM","CD44","NT5E","PDGFA","CASP3","CLU","ALDH1A1","INSL3","PTGDS","CDC25C","ID4","SFRP1","KLF6","RSPO1","LCN2","ITLN1","RGS5","MCAM","CSPG4","MYL9","FRZB","CD36","NOTCH3","ADIRF","CRIP1","GATA4","TAGLN","VIM","CD163","S100A4","TYROBP","LYZ","RGS1","VWF","EPAS1","TGFBR2","NOSTRIN","PALMD","POSTN","SOX9","AMH")

markers=markers[which(markers %in% rownames(emat))]
markers=markers[which(markers %in% rownames(rvel$current))]
length(markers) # 43 30 40

emb=tsne
pdf(file=paste0(ss,"_",set,"_",velname,'fit.tsne.gene.kCells10.pdf'),height=3.2,width=12)
for(gene in markers){
gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 10,kGenes=1,fit.quantile=fit.quantile,cell.emb=emb,cell.colors=cell.colors,show.gene=gene,old.fit=vel,do.par=T)
# cell.dist=cell.dist,
}
dev.off()
# calculating convolved matrices ... done

emb=umap
pdf(file=paste0(ss,"_",set,"_",velname,'fit.umap.gene.kCells10.pdf'),height=3.2,width=12)
for(gene in markers){
gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 10,kGenes=1,fit.quantile=fit.quantile,cell.emb=emb,cell.colors=cell.colors,show.gene=gene,old.fit=vel,do.par=T)
# cell.dist=cell.dist,
}
dev.off()
