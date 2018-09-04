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
NBSmatrixMutation <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_sampleXgenes_smoothed.csv",
                       sep = ",", quote = "", header = FALSE, stringsAsFactors = FALSE)
NBSpatients <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_sample_raw.csv",
                       sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)
NBSgenes <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_geneKeys_smoothed.csv",
                       sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)
NBSpatients <- gsub(',','',as.character(unlist(NBSpatients)))
NBSpatients <- tolower(NBSpatients)
sum(NBSmatrixMutation)
max(NBSmatrixMutation)
min(NBSmatrixMutation)
binaMutation <- NBSmatrixMutation
row.names(binaMutation)  <- as.character(unlist(NBSpatients))
names(binaMutation) <- as.character(unlist(NBSgenes))

#Number of patients with mutation x
pat4genes <- apply(binaMutation,2,sum)
#Genes to remove
if ( length(which(pat4genes<=0)) > 0 ){
  binaMutation <- binaMutation[,-which(pat4genes<=0)]
}

x<-apply(binaMutation,2,sum)
min(x)
save(binaMutation, file="data/matrixMutationSmooth.Rd")

patients <- NBSpatients
save(patients, file="data/patientsSmooth.Rd")


# data from Matlab raw ----
NBSboolMutation <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_sampleXgenes_raw.csv",
                              sep = ",", quote = "",
                              header = FALSE, stringsAsFactors = FALSE)
NBSpatients <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_sample_raw.csv",
                          sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)

NBSgenes <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_genes_raw.csv",
                       sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)

NBSpatients <- gsub(',','',as.character(unlist(NBSpatients)))
NBSpatients <- tolower(NBSpatients)
sum(NBSboolMutation)
max(NBSboolMutation)
boolMutation <- NBSboolMutation
row.names(boolMutation)  <- NBSpatients
names(boolMutation) <- as.character(unlist(NBSgenes))


#Number of patients with mutation x
pat4genes <- apply(boolMutation,2,sum)
#Genes to remove
if ( length(which(pat4genes<=0)) > 0 ){
  boolMutation <- boolMutation[,-which(pat4genes<=0)]
}


x<-apply(boolMutation,2,sum)
min(x)
# source("R/1-TensorDefinition/fromLogic2Numeric.R")
# binaMutation <- fromLogic2Numeric(boolMutation)
binaMutation <- boolMutation
anyNA(binaMutation)
save(binaMutation, file="data/binaMutation.RData")

patients <- NBSpatients
save(patients, file="data/patients.Rd")

source("R/1-TensorDefinition/fromNumeric2Logic.R")
boolMutation <- fromNumeric2LogicFast(boolMutation)
anyNA(boolMutation)
save(boolMutation, file="data/boolMutation.Rd")
#boolMutation[1:5,1:5]


# RTGCA -----------
load("data/loadls.RData")
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

load("data/loadls.RData")
loadls("RTCGA RTCGA.clinical",F)
data(OV.clinical)
names(OV.clinical)
load(file="data/patients.Rd")

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
#NA = gx which means that cannot be assessed
clinicalOV$patient.neoplasm_histologic_grade[which(is.na(clinicalOV$patient.neoplasm_histologic_grade))]="gx"

table(clinicalOV$patient.biospecimen_cqcf.histological_type , useNA =  "always" )
#delete column, it has only one value
which(names(clinicalOV) == "patient.biospecimen_cqcf.histological_type")
clinicalOV <- clinicalOV[,-5]

table(clinicalOV$patient.biospecimen_cqcf.history_of_neoadjuvant_treatment , useNA =  "always" )
#delete column, it has only one value
which(names(clinicalOV) == "patient.biospecimen_cqcf.history_of_neoadjuvant_treatment")
clinicalOV <- clinicalOV[,-5]


table(clinicalOV$patient.stage_event.clinical_stage , useNA =  "always" )
#NA= stage iiic
clinicalOV$patient.stage_event.clinical_stage[which(is.na(clinicalOV$patient.stage_event.clinical_stage))]="stage iiic"

row.names(clinicalOV) <- clinicalOV$patient.bcr_patient_barcode
clinicalOV <- clinicalOV[,-2]

patients <- as.character(row.names(clinicalOV))

save.image("temp/1-2v3-before_dico")
#### dicomotization of clinical vars -------

source('~/CLIGEN_tgit/R/1-TensorDefinition/dicotomFunctions.R')

clinical <- clinicalOV

NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
NumVals

#4 level columns-----------------
cols4 <- which(NumVals==4)
ncols <- names(cols4)
table(clinical[,ncols[1]], useNA =  "always" )

matDich4 <- dichoSeveralCols(cols=cols4)

#OPTIONAL delete one column
#I'm afraid we are getting this error ('MATLAB:sylvester:solutionNotUnique')
#because there is always a redundant column so Im getting rid of them
table(clinical[,ncols[1]], useNA =  "always" )
#In this case delete the less frequent aka alaskan
matDich4 <- matDich4[ , -1]


#5 level columns-----------------
#patient.neoplasm_histologic_grade
cols5 <- which(NumVals==5)
ncols <- names(cols5)
table(clinical[,ncols[1]], useNA =  "always" )

#based on SEER move gb to gx
clinical$patient.neoplasm_histologic_grade[which(clinical$patient.neoplasm_histologic_grade == "gb" )]="gx"
#and too little g1 move to g1 and g2
clinical$patient.neoplasm_histologic_grade[which(clinical$patient.neoplasm_histologic_grade == "g1" )]="g1_2"
clinical$patient.neoplasm_histologic_grade[which(clinical$patient.neoplasm_histologic_grade == "g2" )]="g1_2"

#after merging values it hase three
table(clinical[,ncols[1]], useNA =  "always" )


##Stage
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


#3 level columns-----------------
cols3 <- which(NumVals==3)
ncols <- names(cols3)
table(clinical[,ncols[1]], useNA =  "always" )
table(clinical[,ncols[2]], useNA =  "always" )

matDich3 <- dichoSeveralCols(cols=cols3)

#OPTIONAL delete one column
#I'm afraid we are getting this error ('MATLAB:sylvester:solutionNotUnique')
#because there is always a redundant column so Im getting rid of them
table(clinical[,ncols[1]], useNA =  "always" )
table(clinical[,ncols[2]], useNA =  "always" )
#In this case delete the less frequent aka gx and STAGE I and II
matDich3 <- matDich3[ , -c(3,4)]



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

#OPTIONAL delete one column
#I'm afraid we are getting this error ('MATLAB:sylvester:solutionNotUnique')
#because there is always a redundant column so Im getting rid of them
apply(aMati,2,sum)
#In this case delete the less frequent aka NA(79,89]
aMati <- aMati[ , -10]


#merge all levels ----

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

#WHEN OPTIONAL
save(boolClinical,file="temp/boolCSmooth_beforeSave.RData")


#save data again -----
#clinical
load(file="data/boolClinical.RData")
patients <- row.names(boolClinical)
save(boolClinical, file="data/boolClinical.RData")
save(patients, file="data/patients.Rd")
save(clinicalOV,file="temp/clinicalOV_raw.RData")
source( '~/CLIGEN_tgit/R/1-TensorDefinition/fromLogic2Numeric.R' )
binaClinical <- fromLogic2Numeric(boolClinical[patients,])
save(binaClinical,file="temp/binaClinical.Rd")


#smooth
source( '~/CLIGEN_tgit/R/1-TensorDefinition/fromLogic2Numeric.R' )
binaClinical <- fromLogic2Numeric(boolClinical[patients,])
load("data/matrixMutationSmooth.Rd")
ppp <- intersect(patients, row.names(binaMutation))
binaClinical <- binaClinical[ppp,]
anyNA(binaClinical)
save(binaClinical,file="data/binaClinicalSmooth.Rd")

patients <- ppp
save(patients,file="data/patientsSmooth.Rd")


#mutation
load(file="temp/clinicalOV.RData")
load("data/boolMutation.RData")
boolMutation <- boolMutation[patients,]
anyNA(boolMutation)
save(boolMutation, file="data/boolMutation.RData")
load("data/binaMutation.Rd")
binaMutation <- binaMutation[patients,]
anyNA(binaMutation)
save(binaMutation, file="data/binaMutation.RData")

#smooth
binaMutation <- binaMutation[ppp,]
save(binaMutation, file="data/matrixMutationSmooth.Rd")



#write csv files for rubik matlab ------------
rm(list = ls())
load(file="data/boolClinical.RData")
load(file="data/patients.Rd")
load(file="temp/binaClinical.Rd")
load(file="data/binaMutation.Rd")

#smooth
load(file="data/matrixMutationSmooth.Rd")
load(file="data/patientsSmooth.Rd")
load(file="temp/binaClinicalSmooth.Rd")



tMut <- t(binaMutation[patients,])
anyNA(tMut)
#write.csv2(boolMutation, file = "temp/10-boolM.csv",row.names = FALSE)
write.csv2(tMut, file = "output4rubik/1-binaM_OVA.csv",row.names = FALSE)
write.csv2(binaClinical, file = "output4rubik/1-binaC_OVA.csv",row.names = FALSE)
