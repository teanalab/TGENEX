rm(list = ls())
source('~/CLIGEN_tgit/R/1-TensorDefinition/createBoolMutation.R')
boolMutation <- createBoolMutation(fileName = "/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_mutations_extended.txt")


#write.csv2(data.frame(names(boolMutation)),
#           file = "temp/genes.csv",row.names = FALSE)
#get genes from https://www.genenames.org/cgi-bin/symbol_checker
gene4genenames <-read.table(file="temp/1-geneFromGeneNames.txt", sep = "\t",
                            quote = "",
                            header = TRUE, stringsAsFactors = FALSE,
                            skipNul = TRUE, #skip = 1,
                            blank.lines.skip = TRUE,
                            fill=TRUE)

genes2map <- names(boolMutation)
genes2del <- which(gene4genenames$Match.type == "Unmatched")
genes2map <- genes2map[-genes2del]



#from my gist
library("org.Hs.eg.db")
Sym2Entrez <- get("org.Hs.egSYMBOL2EG")
genes2map <- gene4genenames$Approved.symbol
genes2mapI <- gene4genenames$Input
genes2del <- which(gene4genenames$Match.type == "Unmatched")
genes2map <- genes2map[-genes2del]
genes2mapI <- genes2mapI[-genes2del]
gene4genenames <- gene4genenames[!duplicated(gene4genenames$Input),]

genesEntrez <- rep(0, length(genes2map))
for (Gi in seq_along(genes2map) )
{
  newInstruction <- try( genesEntrez[Gi] <- as.list(Sym2Entrez[genes2map[Gi]])[1],
                         silent=TRUE )
  if( inherits(newInstruction, "try-error") )
  {
    try( genesEntrez[Gi] <- as.list(Sym2Entrez[genes2mapI[Gi]])[1] ,
         silent=FALSE )
  }
}

genesEntrezL <- unlist(genesEntrez)
boolMutationRedu <- boolMutation[,genes2mapI]
#names(boolMutationRedu) <- genes2map
binaMutation <- as.data.frame(lapply(boolMutationRedu, as.numeric), stringsAsFactors = FALSE)
write.csv2(t(binaMutation), file = "temp/10-binaMredu_OVA.csv",row.names = FALSE)
write.csv2(genesEntrezL, file = "temp/10-genes_entrez_OVA.csv",row.names = FALSE)
write.csv2(genes2map, file = "temp/10-genes_redu_OVA.csv",row.names = FALSE)







############ 2018 -07-17 --------


######## trying here same ------
library("org.Hs.eg.db")
entrez2Sym <- get("org.Hs.eg.dbSYMBOL")
entrez2Sym <- get("org.Hs.egSYMBOL")
Sym2Entrez <- get("org.Hs.egSYMBOL2EG")
genes2map <- genes2map[-c("RP11-323I15.2")]
genes2map <- names(boolMutation)
genes2map <- genes2map[-c("RP11-323I15.2")]
genes2map <- genes2map[-"RP11-323I15.2"]
genes2del <- which(genes2map == "RP11-323I15.2")
genes2del
genes2map <- genes2map[-genes2del]
xx <- as.list(Sym2Entrez[genes2map])
genes2del <- which(genes2map %in% c("RP11-323I15.2", "KIRREL")
xx <- as.list(Sym2Entrez[genes2map])
genes2del <- which(genes2map %in% c("RP11-323I15.2", "KIRREL"))
genes2map <- genes2map[-genes2del]
xx <- as.list(Sym2Entrez[genes2map])
for (Gi in seq_along(genes2map) )
{
Sym2Entrez[Gi]
}
for (Gi in seq_along(genes2map) )
{
print ( Sym2Entrez[genes2map[Gi]] )
}
for (Gi in seq_along(genes2map) )
{
try ( Sym2Entrez[genes2map[Gi]] )
}
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "\t", quote = "",
header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "\t", quote = "",
header = TRUE, stringsAsFactors = FALSE, skipNul = TRUE, #skip = 1,
blank.lines.skip = TRUE)
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "|", quote = "",
header = TRUE, stringsAsFactors = FALSE, skipNul = TRUE, #skip = 1,
blank.lines.skip = TRUE)
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "|",
quote = "",
header = TRUE, stringsAsFactors = FALSE,
skipNul = TRUE, #skip = 1,
blank.lines.skip = TRUE,
fill=TRUE)
View(gene4genenames)
which(gene4genenames$Match.type == "Unmatched")
gene4genenames[which(gene4genenames$Match.type == "Unmatched"),1]




####### 2018-07-15 ---------

write.csv2(patients, file = "temp/10-patients.csv",row.names = FALSE)
row.names(boolMutation)
col.names(boolMutation)
names(boolMutation)
write.csv2(names(boolMutation), file = "temp/10-genes.csv",row.names = FALSE)
load(file= "data/boolClinical.RData")
load("data/tensorClinical.RData")
View(tensorClinical)
clinical <- tensorClinical
mutation <-read.table(file="data/brca4mCBioportal/data_mutations_extended.txt", sep = "\t", quote = "",
header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)
#uuid = NO A GOOD IDENTIFIER
#Tumor_Sample_Barcode  = GOOD IDENTIFIER, BUT WE NEED TO CUT
#TCGA-06-5411-01A-01D-1696-08
mutationBarcodes <- substring(tolower(mutation$Tumor_Sample_Barcode), 1, 12)
mutation$mutationBarcodes <- mutationBarcodes
patients <<- intersect(mutationBarcodes, clinical$patient.bcr_patient_barcode)
######### start here ----
rm(list = ls())
load("data/tensorClinical.RData")
clinical <- tensorClinical
mutation <-read.table(file="data/brca4mCBioportal/data_mutations_extended.txt", sep = "\t", quote = "",
header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)
mutation <-read.table(file="data/brca4mCBioportal/data_mutations_extended.txt", sep = "\t", quote = "",
header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)
#uuid = NO A GOOD IDENTIFIER
#Tumor_Sample_Barcode  = GOOD IDENTIFIER, BUT WE NEED TO CUT
#TCGA-06-5411-01A-01D-1696-08
mutationBarcodes <- substring(tolower(mutation$Tumor_Sample_Barcode), 1, 12)
mutation$mutationBarcodes <- mutationBarcodes
patients <<- intersect(mutationBarcodes, clinical$patient.bcr_patient_barcode)
View(clinical)
write.csv2(clinical, file = "temp/10-clinical.csv",row.names = FALSE)
names(boolMutation)
load(file="data/boolMutation.RData")

######## trying here NOT GOOD ------
names(boolMutation)

library("org.Hs.eg.db")
entrez2Sym <- get("org.Hs.egSYMBOL")
Sym2Entrez <- get("org.Hs.egSYMBOL2EG")
mapped_genes <- mappedkeys(Sym2Entrez)
xx <- as.list(Sym2Entrez[mapped_genes])

genes2map <- genes2map[-c("RP11-323I15.2")]

genes2map <- names(boolMutation)
genes2map <- genes2map[-c("RP11-323I15.2")]
genes2map <- genes2map[-"RP11-323I15.2"]
genes2del <- which(genes2map == "RP11-323I15.2")
genes2del
genes2map <- genes2map[-genes2del]
xx <- as.list(Sym2Entrez[genes2map])
genes2del <- which(genes2map %in% c("RP11-323I15.2", "KIRREL")
xx <- as.list(Sym2Entrez[genes2map])
genes2del <- which(genes2map %in% c("RP11-323I15.2", "KIRREL"))
genes2map <- genes2map[-genes2del]
xx <- as.list(Sym2Entrez[genes2map])
for (Gi in seq_along(genes2map) )
{
Sym2Entrez[Gi]
}
for (Gi in seq_along(genes2map) )
{
print ( Sym2Entrez[genes2map[Gi]] )
}
for (Gi in seq_along(genes2map) )
{
try ( Sym2Entrez[genes2map[Gi]] )
}
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "\t", quote = "",
header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "\t", quote = "",
header = TRUE, stringsAsFactors = FALSE, skipNul = TRUE, #skip = 1,
blank.lines.skip = TRUE)
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "|", quote = "",
header = TRUE, stringsAsFactors = FALSE, skipNul = TRUE, #skip = 1,
blank.lines.skip = TRUE)
gene4genenames <-read.table(file="data/1-geneFromGeneNames", sep = "|",
quote = "",
header = TRUE, stringsAsFactors = FALSE,
skipNul = TRUE, #skip = 1,
blank.lines.skip = TRUE,
fill=TRUE)







