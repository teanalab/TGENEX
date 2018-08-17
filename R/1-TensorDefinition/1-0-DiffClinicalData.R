# what happened with these:?
#   #cols3 <- cols3[-9] # icd_o_3_site
#   #cols3 <- cols3[-8] # repited with icd_o_3_site
#   #cols3 <- cols3[-7] # histological_type
#
# Those are GBM vaiables

  rm(list = ls())
tcga_pub_clinical<-read.table(file="../../12Datasets/BRCA/brca4mCBioportal/brca_tcga_pub_clinical_data.tsv", sep = "\t",
                     header = TRUE, stringsAsFactors = FALSE)

data_clinical<-read.table(file="../../12Datasets/BRCA/brca4mCBioportal/data_clinical.txt", sep = "\t",
                              header = TRUE, stringsAsFactors = FALSE)

clinical<-read.table(file="../../12Datasets/BRCA/clinical.txt", sep = "\t",
                     header = TRUE, stringsAsFactors = FALSE)


names(tcga_pub_clinical)

names(data_clinical)
