### Match 1-1-1 Orthologues for Human-Monkey-Mouse gene names
# 2.5.2019 by Qianyi
### note: I've already replaced Monkey Ensembl ID by non-duplicated Human Gene Symbol with one-to-one orthologue for single species analysis for Monkey

setwd("C:/Users/qzm/Desktop/HumanMonkeyST/3SpeciesComparison_plot")

# Orthologue Genes
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

### Human-Macaque
library(biomaRt) 
listEnsemblArchives()
listMarts() # # Ensembl Gene 95
# Use Ensembl Gene 94
ensembl <- useMart(host='http://oct2018.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL')
#ensembl = useMart("ENSEMBL_MART_ENSEMBL") # Ensembl Gene 95
#ensembl = useEnsembl(biomart="ensembl",version=94)
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
#listFilters(ensembl) #265 with_mmulatta_homolog
#searchAttributes(mart = ensembl, pattern = "mmulatta")

attributes = c("ensembl_gene_id","external_gene_name","mmulatta_homolog_ensembl_gene","mmulatta_homolog_associated_gene_name","mmulatta_homolog_perc_id","mmulatta_homolog_perc_id_r1","mmulatta_homolog_goc_score","mmulatta_homolog_wga_coverage","mmulatta_homolog_dn","mmulatta_homolog_ds","mmulatta_homolog_orthology_confidence","mmulatta_homolog_orthology_type")

#Get all Human-Macaque gene pairs
humMac <- getBM(attributes = attributes,mart = ensembl)
dim(humMac)      # [1] 66867    12

#Remove non-homologs
humMacOrth <- humMac[humMac$mmulatta_homolog_ensembl_gene!='',]
dim(humMacOrth)  # [1] 23871    12

#Confidence
table(humMacOrth[,11]) # "mmulatta_homolog_orthology_confidence"
#    0     1 
# 2976 20895 
sum(humMacOrth$mmulatta_homolog_orthology_confidence == 1)
#20895  High confidence
sum(humMacOrth$mmulatta_homolog_orthology_confidence == 0)
#2976  Low confidence


### Human-Mouse
#searchAttributes(mart = ensembl, pattern = "mmusculus")
attributes = c("ensembl_gene_id","external_gene_name","mmusculus_homolog_ensembl_gene","mmusculus_homolog_associated_gene_name", "mmusculus_homolog_perc_id","mmusculus_homolog_perc_id_r1","mmusculus_homolog_goc_score","mmusculus_homolog_wga_coverage","mmusculus_homolog_dn","mmusculus_homolog_ds","mmusculus_homolog_orthology_confidence","mmusculus_homolog_orthology_type")

humMou <- getBM(attributes = attributes,mart = ensembl)
dim(humMou)  # [1] 71550    12

#Remove non-homologs
humMouOrth <- humMou[humMou$mmusculus_homolog_ensembl_gene!='',]
dim(humMouOrth)  # [1] 26312    12

#Confidence
table(humMouOrth[,11]) # "mmulatta_homolog_orthology_confidence"
#    0     1 
# 8840 17472
sum(humMouOrth$mmusculus_homolog_orthology_confidence == 1)
#17472  High confidence
sum(humMouOrth$mmusculus_homolog_orthology_confidence == 0)
#8840  Low confidence


### Mouse-Macaque
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
#searchAttributes(mart = ensembl, pattern = "mmulatta")

attributes = c("ensembl_gene_id","external_gene_name","mmulatta_homolog_ensembl_gene","mmulatta_homolog_associated_gene_name", "mmulatta_homolog_perc_id","mmulatta_homolog_perc_id_r1","mmulatta_homolog_goc_score","mmulatta_homolog_wga_coverage","mmulatta_homolog_dn","mmulatta_homolog_ds","mmulatta_homolog_orthology_confidence","mmulatta_homolog_orthology_type")

#Get all Human-Macaque gene pairs
mouMac <- getBM(attributes = attributes,mart = ensembl)
dim(mouMac)      # [1] 61306    12

#Remove non-homologs
mouMacOrth <- mouMac[mouMac$mmulatta_homolog_ensembl_gene!='',]
dim(mouMacOrth)  # [1] 25112    12

#Confidence
table(mouMacOrth[,11]) # "mmulatta_homolog_orthology_confidence"
#    0     1 
# 11273 13839
sum(mouMacOrth$mmulatta_homolog_orthology_confidence == 1)
#13839  High confidence
sum(mouMacOrth$mmulatta_homolog_orthology_confidence == 0)
#11273  Low confidence


### Human-Macaque-Mouse
humMacMouOrth <- merge(humMacOrth,humMouOrth,by = "ensembl_gene_id")


### replace Ensembl Gene IDs with no Gene Symbol in macaque with non-duplicated human Gene Symbols
emprows = which(humMacMouOrth$mmulatta_homolog_associated_gene_name == '' & !(humMacMouOrth$external_gene_name.x %in% humMacMouOrth$mmulatta_homolog_associated_gene_name))
length(emprows) # [1] 3073
humMacMouOrth$mmulatta_homolog_associated_gene_name[emprows] = humMacMouOrth$external_gene_name.x[emprows]


#### One-to-one-to-one genes
humMacMouOrth.1to1 <- humMacMouOrth[humMacMouOrth$mmulatta_homolog_orthology_type=="ortholog_one2one",]
humMacMouOrth.1to1 <- humMacMouOrth.1to1[humMacMouOrth.1to1$mmusculus_homolog_orthology_type=="ortholog_one2one",]
humMacMouOrth.1to1 <- humMacMouOrth.1to1[humMacMouOrth.1to1$mmulatta_homolog_associated_gene_name!="",]

colnames(humMacMouOrth.1to1)
# [2] "external_gene_name.x"                  
# [4] "mmulatta_homolog_associated_gene_name" 
# [13] "external_gene_name.y" 
# [15] "mmusculus_homolog_associated_gene_name"
dim(humMacMouOrth.1to1) # [1] 14743    23
which(humMacMouOrth.1to1$external_gene_name.x != humMacMouOrth.1to1$external_gene_name.y)
# integer(0)

Orth1to1to1=humMacMouOrth.1to1[,c(2,4,15)]
Orth1to1to1=unique(Orth1to1to1)
dim(Orth1to1to1) # [1] 14748     
length(unique(Orth1to1to1[,1])) # [1] 14718
length(unique(Orth1to1to1[,2])) # [1] 14398
length(unique(Orth1to1to1[,3])) # [1] 14742

Orth1to1to1[anyDuplicated(Orth1to1to1[,3]),3]  # Gm26265
Orth1to1to1[which(Orth1to1to1[,3]=="Gm26265"),] 
humMac[which(humMac[,2]=="SNORA50A"),c(2,4)]

write.table(Orth1to1to1,"humMacMouOrth.1to1to1.txt",row.names=F,col.names=F,quote=F,sep="\t")
