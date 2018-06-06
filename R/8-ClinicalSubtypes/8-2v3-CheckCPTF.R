#' version 3 (Jun/6/2018)
#' Description:
#' Check if the subtypes obtained from 8.1 match anything that we got from CP TF

rm(list = ls())
load("data/myLib.RData")
load("temp/6-4v1.RData")
library("caret")
library("e1071")
CompToBC = read.csv(file="R/8-ClinicalSubtypes/CompVsBCtype.csv")

patiAfi_NTF <- patientAfiliation_NTF
row.names(survClinical) <- survClinical$patient.bcr_patient_barcode
PAM50_surv_c <- survClinical[row.names(patiAfi_NTF), c("patient.bcr_patient_barcode", "PAM50.Subtype")]
patiAfi_NTF <- merge(as.data.frame(PAM50_surv_c), as.data.frame(patientAfiliation_NTF), by='row.names', all=TRUE)
patiAfi_NTF <- patiAfi_NTF[,-1]
row.names(patiAfi_NTF) <- patiAfi_NTF$patient.bcr_patient_barcode
patiAfi_NTF <- patiAfi_NTF[,-1]
# delete subtypes that are not predicted
toDel <- row.names(patiAfi_NTF[patiAfi_NTF$PAM50.Subtype == "Normal-like",])
toDel1 <- row.names(patiAfi_NTF[patiAfi_NTF$PAM50.Subtype == "HER2-enriched",])
patiAfi_NTF <- patiAfi_NTF[-which(row.names(patiAfi_NTF) %in% c(toDel,toDel1)),]


realpSubtype <- factor(patiAfi_NTF$PAM50.Subtype)
predictedpSubtype <- patiAfi_NTF$component
predictedpSubtype <- CompToBC[predictedpSubtype,2]
confusionMatrix(realpSubtype,predictedpSubtype)

