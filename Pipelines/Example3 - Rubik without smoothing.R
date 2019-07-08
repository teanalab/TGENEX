#Example of loading manually downloaded files
rm(list = ls())

source("R/0-ConfigurationVariables.R")
minNumOfMutationsPerPatient = 5
minNumOfPatientsPerGene = 5

source("R/1-ReadRawFiles.R")
# 1- read the raw files from the 'rawData' folder and store them in 'dataFile'
readBRCAFiles (dataFile = "temp/readBRCAFiles.RData",
               rawMutationFile="rawData/data_mutations_extended.txt",
               rawClinicalFile="rawData/BRCA.clin.merged.txt")

# 2- filter the clinical variables specified in "rawData/clinical_useful.csv" and "rawData/merge_useful.csv"
filterClinicalVariables()
load("data/tensorClinical.RData")

source("R/1-createBoolMutation.R")
createBoolMutation(boolOrBin="bina")
load("data/patients.RData")
load("data/genes.RData")
load("data/binaMutation.RData")

#Filter tensor clinical
load("data/tensorClinical.RData")
row.names(tensorClinical) <- tensorClinical$patient.bcr_patient_barcode
tensorClinical <- tensorClinical[patients, ]
save(tensorClinical, file="data/tensorClinical.RData")


#Dichotomize Clinical vars


# Step 2 - Construct the tensor
source("R/2-BoolMatrices.R")
#source("R/2-Smooth_MutationM_SVD.R")
#source("R/2-Construct_tensor.R")

# Step 3 - Tensor Factorization
#source("R/3-RunTF.R")
