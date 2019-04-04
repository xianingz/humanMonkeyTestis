### Convert Monkey Ensembl IDs without gene symbols to nonduplicated Human Gene Symbols with 1-1 orthologues
# 1.28.2019 by Qianyi
### note:
### Monkey Ensembl IDs are not well annotated with monkey gene symbols 
### A lot of Monkey Ensembl Gene IDs do not have annotated gene symbol
### Decided to Replace Monkey Ensembl ID without gene symbols by non-duplicated Human Gene Symbols with one-to-one orthologue to monkey


# Orthologue Genes
home="/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/monkey"
setwd(home)
library(biomaRt) 
listEnsemblArchives()
listMarts() # # Ensembl Gene 95

### Human-Macaque
# Use Ensembl Gene 93 for Monkey
### note: I used Ensembl gene 93 to map and annotate genes for Monkey
ensembl <- useMart(host='http://jul2018.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL')
#ensembl = useMart("ENSEMBL_MART_ENSEMBL") # Ensembl Gene 95
#ensembl = useEnsembl(biomart="ensembl",version=94)
listDatasets(ensembl)
ensembl = useDataset("mmulatta_gene_ensembl",mart=ensembl)
#listFilters(ensembl) #265 with_mmulatta_homolog
#searchAttributes(mart = ensembl, pattern = "mmulatta")

attributes = c("ensembl_gene_id","external_gene_name","hsapiens_homolog_ensembl_gene","hsapiens_homolog_associated_gene_name","hsapiens_homolog_perc_id","hsapiens_homolog_perc_id_r1","hsapiens_homolog_goc_score","hsapiens_homolog_wga_coverage","hsapiens_homolog_dn","hsapiens_homolog_ds","hsapiens_homolog_orthology_confidence","hsapiens_homolog_orthology_type")

#Get all Monkey-Human gene pairs
humMac <- getBM(attributes = attributes,mart = ensembl)
dim(humMac)      # [1] 34372    12

#Remove non-homologs
humMacOrth <- humMac[humMac$hsapiens_homolog_ensembl_gene!='',]
dim(humMacOrth)  # [1] 23907    12

#Confidence
table(humMacOrth[,11]) # "hsapiens_homolog_orthology_confidence"
#    0     1 
# 2935 20972 
sum(humMacOrth$hsapiens_homolog_orthology_confidence == 1)
#20972  High confidence
sum(humMacOrth$hsapiens_homolog_orthology_confidence == 0)
#2935  Low confidence


### replace Ensembl Gene IDs with no Gene Symbol in macaque by human Gene Symbols with one-to-one orthologue
humMacOrth1to1=humMacOrth[humMacOrth$hsapiens_homolog_orthology_type=="ortholog_one2one",]
id = which(humMacOrth1to1$external_gene_name == '')
macGene=humMacOrth1to1[,1:2]
macGene[id,2]<-humMacOrth1to1[id,4]
macGene[id,2]
dim(macGene) # [1] 19810     2
length(id)   # converted 560 monkey Ensembl IDs without gene symbols to Human gene symbols with 1-to-1 ortholog
write.table(macGene,"Monkey_EnsemblIDconvertedtoHumanGenes1to1Ortholog_Ensembl93.txt",col.names=F,row.names=F,quote=F,sep="\t")

### double-check the example below
### note: some human gene symbols with 1-1 orthologues to monkey Ensembl IDs are duplicated from existing Monkey gene symbols annotated for another Ensembl ID
### my default: use monkey gene symbols for monkey Enseml IDs if there is any duplicated gene symbols from human
humMac[which(humMac[,1]=="ENSMMUG00000004552"),c(1:4,12)]
humMac[which(humMac[,1]=="ENSMMUG00000004554"),c(1:4,12)]
         ensembl_gene_id external_gene_name
13969 ENSMMUG00000004552                DDT
      hsapiens_homolog_ensembl_gene
13969                              
      hsapiens_homolog_associated_gene_name
13969                                      
      hsapiens_homolog_orthology_type
13969                                

         ensembl_gene_id external_gene_name
14353 ENSMMUG00000004554                   
      hsapiens_homolog_ensembl_gene
14353               ENSG00000099977
      hsapiens_homolog_associated_gene_name
14353                                   DDT
      hsapiens_homolog_orthology_type
14353                ortholog_one2one



