# version 3 (8/8/2018)
#1) patients with fewer than 10 mutations were discarded initially, but after trimming, remaining patients have three or more mutations
#2) only important variables
#3) using cBio data
#4) using firebrowse and cBiolite data

##### OVARIAN ##########
rm(list = ls())
#version 2 (Aug/2/2018)

######### mutation ----
source('~/CLIGEN_tgit/R/1-TensorDefinition/createBoolMutation.R')

boolMutation <- createBoolMutation(fileName = "/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_mutations_extended.txt",
                                   boolOrBin="bina")
#description boolean table
genes <- names(boolMutation)
patients <- row.names(boolMutation)
length(genes)
length(patients)
dim(boolMutation)
#save(boolMutation, patients, file="temp/boolMutation_OVA.RData")
#save(boolMutation, patients, file="temp/binaMutation_OVA.RData")
load("temp/binaMutation_OVA.RData")


##### compare different sources of OV Data -----

load("loadls.RData")

# firebrowse
load( file="temp/binaMutation_OVA.RData")
mutatFirebrowse <- boolMutation
patientFirebrowse <- patients
sum(mutatFirebrowse)

# data from Matlab Smoothed network----
NBSboolMutation <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_sampleXgenes.csv",
                       sep = ",", quote = "",
                       header = FALSE, stringsAsFactors = FALSE)
NBSpatients <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_output_sample_id4.csv",
                       sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)

NBSgenes <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_genes.csv",
                          sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)


sum(NBSboolMutation)
max(NBSboolMutation)
boolMutation <- NBSboolMutation
row.names(boolMutation)  <- as.character(unlist(NBSpatients))
names(boolMutation) <- as.character(unlist(NBSgenes))

source("R/1-TensorDefinition/fromNumeric2Logic.R")
boolMutation <- fromNumeric2Logic(boolMutation) #takes time



#Number of patients with mutation x
pat4genes <- apply(boolMutation,2,sum)
#Genes to remove
if ( length(which(pat4genes<=0)) > 0 ){
  boolMutation <- boolMutation[,-which(pat4genes<=0)]
}

x<-apply(boolMutation,2,sum)
min(x)
save(boolMutation, file="data/boolMutationSmooth.RData")


source("R/1-TensorDefinition/fromLogic2Numeric.R")
binaMutation <- fromLogic2Numeric(boolMutation)
save(binaMutation, file="data/binaMutationSmooth.Rd")

# data from Matlab raw ----
NBSboolMutation <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_sampleXgenes_raw.csv",
                              sep = ",", quote = "",
                              header = FALSE, stringsAsFactors = FALSE)
NBSpatients <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_sample_raw.csv",
                          sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)

NBSgenes <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_genes_raw.csv",
                       sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)


sum(NBSboolMutation)
max(NBSboolMutation)
boolMutation <- NBSboolMutation
row.names(boolMutation)  <- as.character(unlist(NBSpatients))
names(boolMutation) <- as.character(unlist(NBSgenes))

source("R/1-TensorDefinition/fromNumeric2Logic.R")
boolMutation <- fromNumeric2Logic(boolMutation) #takes time



#boolMutation[1:5,1:5]





# RTGCA -----------
loadls("RTCGA.mutations dplyr reshape2")
library(RTCGA.mutations)
library(dplyr)
mutationsTCGA(OV.mutations) %>%
  filter(Variant_Classification != "Silent") %>% # cancer tissue
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
  mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> TCGOV.mutations

patients <- unique(TCGOV.mutations$bcr_patient_barcode)
genes <- unique(TCGOV.mutations$Hugo_Symbol)

binaMutation <- data.frame(matrix(0,length(patients),length(genes)),
                           row.names = patients)
names(binaMutation) <- genes

for(i in seq_len(dim(TCGOV.mutations)[1])){
  Apatient <- as.character(TCGOV.mutations$bcr_patient_barcode[i])
  APgene <- as.character(TCGOV.mutations$Hugo_Symbol[i])
  binaMutation[Apatient,APgene] <- 1
}

sum(binaMutation)



#### clinical -----

load("loadls.RData")
loadls("RTCGA")
data(OV.clinical)
names(OV.clinical)


names(OV.clinical)[10:20]
clinicalOV <- OV.clinical[which(OV.clinical[,"patient.bcr_patient_barcode"] %in% as.character(unlist(patients)) ) , ]

#remove all NAs
delCols <- apply(clinicalOV,2, function(x){ all(is.na(x)) } )
delC <- which(delCols)
clinicalOV <- clinicalOV[,-delC]

#Num different values
numDiffCols <- apply(clinicalOV, 2, function(x){ length(unique(x)) } )
delC <- which(numDiffCols==1)      #only one value
clinicalOV <- clinicalOV[,-delC]

numDiffCols <- apply(clinicalOV, 2, function(x){ length(unique(x)) } )
delC <- which(numDiffCols==355)      #all values different
clinicalOV <- clinicalOV[,-delC]

delC <- which( regexpr( "file", names(clinicalOV) ) != -1 )
clinicalOV <- clinicalOV[,-delC]

numDiffCols <- apply(clinicalOV, 2, function(x){ length(unique(x)) } )
max(numDiffCols)
delC <- which(numDiffCols==354)      #all values different
clinicalOV <- clinicalOV[,-delC]

# numDiffCols <- apply(clinicalOV, 2, function(x){ length(unique(x)) } )
# max(numDiffCols)
# delC <- which(numDiffCols==354)      #all values different
# clinicalOV <- clinicalOV[,-delC]


#COMPLETE TO KEEP
keepCols <- apply(clinicalOV,2, function(x){ all(!is.na(x)) } )
keepC <- which(keepCols)
clinicalOVcomplete <- clinicalOV[,keepC]


usefulCols <- c("patient.age_at_initial_pathologic_diagnosis", "patient.bcr_patient_barcode")


toKeepC <- which( regexpr( "race", names(clinicalOV) ) != -1 )
names(clinicalOV)[toKeepC]

usefulCols <- c(usefulCols, "patient.race")

toKeepC <- which( regexpr( "histo", names(clinicalOV) ) != -1 )
names(clinicalOV)[toKeepC]
usefulCols <- c(usefulCols,"patient.neoplasm_histologic_grade" ,
                "patient.biospecimen_cqcf.histological_type",
                "patient.biospecimen_cqcf.history_of_neoadjuvant_treatment" )

toKeepC <- which( regexpr( "stage", names(clinicalOV) ) != -1 )
names(clinicalOV)[toKeepC]
usefulCols <- c(usefulCols,"patient.stage_event.clinical_stage")


## get useful cols

clinicalOV <- OV.clinical[which(OV.clinical[,"patient.bcr_patient_barcode"] %in% as.character(unlist(patients)) ) , usefulCols]
keepCols <- apply(clinicalOV,2, anyNA )
keepCols

table(clinicalOV$patient.race,useNA =  "always" )
#NA = white

clinicalOV$patient.race[which(is.na(clinicalOV$patient.race))] = "white"

table(clinicalOV$patient.neoplasm_histologic_grade, useNA =  "always" )
#NA = g3
clinicalOV$patient.neoplasm_histologic_grade[which(is.na(clinicalOV$patient.neoplasm_histologic_grade))]="g3"

table(clinicalOV$patient.biospecimen_cqcf.histological_type , useNA =  "always" )
#delete only one value
which(names(clinicalOV) == "patient.biospecimen_cqcf.histological_type")
clinicalOV <- clinicalOV[,-5]

table(clinicalOV$patient.biospecimen_cqcf.history_of_neoadjuvant_treatment , useNA =  "always" )
#delete only one value
which(names(clinicalOV) == "patient.biospecimen_cqcf.history_of_neoadjuvant_treatment")
clinicalOV <- clinicalOV[,-5]


table(clinicalOV$patient.stage_event.clinical_stage , useNA =  "always" )
#NA= stage iiic
clinicalOV$patient.stage_event.clinical_stage[which(is.na(clinicalOV$patient.stage_event.clinical_stage))]="stage iiic"

row.names(clinicalOV) <- clinicalOV$patient.bcr_patient_barcode
clinicalOV <- clinicalOV[,-2]

patients <- as.character(row.names(clinicalOV))

save(clinicalOV,file="temp/clinicalOV.RData")
load(file="temp/clinicalOV.RData")
#### dicomotization of clinical vars -------

source('~/CLIGEN_tgit/R/1-TensorDefinition/dicotomFunctions.R')

clinical <- clinicalOV

NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
NumVals



#4 level columns-----------------
cols4 <- which(NumVals==4)
ncols <- names(cols4)
table(clinical[,ncols[1]], useNA =  "always" )
table(clinical[,ncols[2]], useNA =  "always" )

matDich4 <- dichoSeveralCols(cols=cols4)


#Stage
table(clinical$patient.stage_event.clinical_stage , useNA =  "always" )
clinical$patient.stage_event.clinical_stage[which(regexpr("stage iii",
                                                          clinical$patient.stage_event.clinical_stage) != -1)] = "STAGE III"
clinical$patient.stage_event.clinical_stage[which(regexpr("stage ii",
                                                            clinical$patient.stage_event.clinical_stage) != -1)] = "STAGE I and II"
clinical$patient.stage_event.clinical_stage[which(regexpr("stage iv",
                                                          clinical$patient.stage_event.clinical_stage) != -1)] = "STAGE IV"
clinical$patient.stage_event.clinical_stage[which(regexpr("stage i",
                                                            clinical$patient.stage_event.clinical_stage) != -1)] = "STAGE I and II"


NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
NumVals
cols3 <- which(NumVals==3)
matDich3 <- dichoSeveralCols(cols=cols3)

#using Sturges
##2	Diagnosis.Age
C <- "patient.age_at_initial_pathologic_diagnosis"
aCOl <- as.numeric(as.character(clinical[patients,C]))
range(aCOl) #34 89
aMat <- classesSturges(C)


#to Logical
aMati <- as.data.frame(lapply(aMat, as.logical))
names(aMati) <- names(aMat)
row.names(aMati) <- patients

conc1 <- merge(matDich4, matDich3, by='row.names' )
row.names(conc1) <- conc1$Row.names
conc1 <- conc1[,-1]

conc2 <- merge(aMati, conc1, by='row.names' )
row.names(conc2) <- conc2$Row.names
conc2 <- conc2[,-1]


boolClinical <- conc2
anyNA(boolClinical)
dim(boolClinical)
row.names(boolClinical) <- row.names(clinical)

patients <- row.names(boolClinical)
save(boolClinical, file="data/boolClinical.RData")
save(patients, file="data/patients.RData")

#write csv files for rubik matlab ------------
rm(list = ls())
load(file="data/binaMutation.Rd")
load(file="data/boolClinical.RData")
load(file="data/patients.RData")
source( '~/CLIGEN_tgit/R/1-TensorDefinition/fromLogic2Numeric.R' )
binaClinical <- fromLogic2Numeric(boolClinical[patients,])
tMut <- t(binaMutation[patients,])
#write.csv2(boolMutation, file = "temp/10-boolM.csv",row.names = FALSE)
write.csv2(tMut, file = "output4rubik/1-binaM_OVA.csv",row.names = FALSE)
write.csv2(binaClinical, file = "output4rubik/1-binaC_OVA.csv",row.names = FALSE)


