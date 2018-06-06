#' version 3 (Jun/6/2018)
#' Description:
#' With 2 luminal A groups: Luminal A PR+ and Luminal A PR-



rm(list = ls())
load("data/myLib.RData")
load("temp/6-4v1.RData")
library("caret")
library("e1071")
CompToBC = read.csv(file="R/8-ClinicalSubtypes/CompVsBCtype_2Luminals.csv")

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

#diferentiate Luminal A PR+ from Luminal A PR-

toAddPRstatus <- tensorClinical[row.names(patiAfi_NTF),]
patiAfi_NTF <- merge(as.data.frame(patiAfi_NTF), as.data.frame(toAddPRstatus), by='row.names', all=TRUE)
row.names(patiAfi_NTF) <- patiAfi_NTF$patient.bcr_patient_barcode
patiAfi_NTF <- patiAfi_NTF[,-1]
table(patiAfi_NTF[which(patiAfi_NTF$PAM50.Subtype == "Luminal A"),"patient.breast_carcinoma_progesterone_receptor_status"])
#delete indeterminate
toDel <- row.names(patiAfi_NTF[patiAfi_NTF$patient.breast_carcinoma_progesterone_receptor_status == "indeterminate",])
patiAfi_NTF <- patiAfi_NTF[-which(row.names(patiAfi_NTF) %in% toDel),]
#PR +
newLuminalA <- apply(patiAfi_NTF,1, function(x){
    if(x["PAM50.Subtype"] == "Luminal A")
      if(x["patient.breast_carcinoma_progesterone_receptor_status"] == "negative")
      {
        return ("Luminal A PR-")
      } else {
          return ("Luminal A PR+")
      }
    return (x["PAM50.Subtype"])
  } )
patiAfi_NTF$PAM50.Subtype <- newLuminalA


realpSubtype <- factor(patiAfi_NTF$PAM50.Subtype)
predictedpSubtype <- patiAfi_NTF$component
predictedpSubtype <- CompToBC[predictedpSubtype,2]
confusionMatrix(realpSubtype,predictedpSubtype)

