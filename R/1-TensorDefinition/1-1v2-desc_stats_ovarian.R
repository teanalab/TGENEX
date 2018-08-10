# version 2 (8/8/2018)
# MERGE FIREBROWSE AND cBiolite data
# Save new clinical matrix
# ovarian

rm(list = ls())

readFiles <- function(){
  mutation<-read.table(file="/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_mutations_extended.txt", sep = "\t", quote = "",
                       header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)

  #mutationBarcodes <- substring(mutation$Tumor_Sample_Barcode, 1, 12)
  #mutationWithBC <- cbind(mutation,mutationBarcodes)

  #DATASETS DESCRIPTION
  dim(mutation)
  mutationRange <- t(sapply(mutation, range))
  #write.csv(mutationRange,file="descriptionData/brca/4cBiomutationRange.csv")

  #THIS FILE WAS TESTED AND IT HAS 2400 COLUMNS BUT ZERO INFORMATION
  # clinicalLong<-read.table(file="../../12Datasets/BRCA/gdacLevel_1/BRCA.merged_only_biospecimen_clin_format.txt", sep = "\t",
  #                          header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  # [16] "patient.patient_id"

  clinicalLongNom<-read.table(file="/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_clinical_patient.txt", sep = "\t",
                           header = FALSE, stringsAsFactors = FALSE, na.strings = "nasrting", fill = TRUE)

  clinicalLong<-read.table(file="/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_clinical_sample.txt", sep = "\t",
                           header = FALSE, stringsAsFactors = FALSE,  fill = TRUE)


  row.names(clinicalLong) <- clinicalLong[,1]
  VariablesLong <- row.names(clinicalLong)
  colsWithID <- which(regexpr("barcode",VariablesLong, ignore.case = TRUE) > 0)
  VariablesLong[colsWithID[1:10]]
  # [16] patient.bcr_patient_barcode
  clinicalLong[16,2:5]
  names(clinicalLong) <- clinicalLong[1,]
  clinicalLong <- clinicalLong[-1,]

  mutation$Tumor_Sample_Barcode[1:5] #"Tumor_Sample_Barcode"
  #unique(mutation$Tumor_Sample_UUID) #Tumor_Sample_UUID
  substring(mutation$Tumor_Sample_Barcode[1:5],1,12)

  patientsIDinMut <- tolower(substring(mutation$Tumor_Sample_Barcode,1,12)) #get just patient ID
  length(patientsIDinMut)
  length(unique(patientsIDinMut))
#
#   patientsIDinClin <- clinicalLong[16,] #get just patient ID
#   length(patientsIDinClin)
#   length(unique(patientsIDinClin))
#
#   #number of patients
#   length(intersect(patientsIDinClin, patientsIDinMut))
#   patientsID <- intersect(patientsIDinClin, patientsIDinMut)
#
#   n <- which(patientsIDinClin %in% patientsID)
#   clinicalLong <- clinicalLong[,n]
#   clinicalLong <- as.data.frame(t(clinicalLong))
#
#   # newpatients <- clinicalLong[,16]
#   # length(intersect(patientsIDinClin, newpatients))
#
#   #clinicalShort <- clinical1[which(clinical1$Patient.ID %in% patients),]
#
#
#   rm(list=ls()[-which(ls() %in% c("mutation", "clinicalLong"))])
#
#
#
#   libs<-c("Packages.R")
#   libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
#   sapply(libs, function(u) {source(u)})
#
#
#
#   write.csv(clinicalLong, "temp/clinicalLong.csv")
#
#   #save.image(file = "temp/1.2-data.RData")
#
# }
#
# load(file = "temp/1.2-data.RData")
# loadlib("reporttools",FALSE)
#
#
#
# ############ analyze variables in a Excel ------
# clini_use_cols<-read.table(file="descriptionData/clinical_useful.csv", sep = ",",
#                            header = TRUE, stringsAsFactors = FALSE, na.strings = "nasrting", fill = TRUE)
#
# useClini <- clinicalLong[, which(clini_use_cols[,2]=="Y")]
# row.names(useClini) <- useClini$patient.bcr_patient_barcode
#
#
#
# nominalCols <- clini_use_cols$column[which(clini_use_cols$type == "nominal")]
# cap1 <- "Patient characteristics: nominal variables. firebrowse"
# write(tableNominal(vars = useClini[,nominalCols], cap = cap1, vertical = TRUE, lab = "tab: nominal1",
#                    longtable = TRUE), file = "temp/1-1.2-Table1cNominalfb.tex")
#
#
# continousCols <- clini_use_cols$column[which(clini_use_cols$type == "continuos")]
# numClini <- useClini[,continousCols]
# indx <- sapply(numClini, is.factor)
# numClini[indx] <- lapply(numClini[indx], function(x) as.numeric(as.character(x)))
# cap1 <- "Patient characteristics: continuous variables.firebrowse"
# write(tableContinuous(vars = numClini, cap = cap1, lab = "tab: continuous2",
#                       longtable = TRUE), file = "temp/1-1.2-Table2cContifb.tex")
#
# ####### read clinical
# clinical1<-read.table(file="../../12Datasets/BRCA/brca4mCBioportal/data_clinical.txt", sep = "\t",
#                       header = TRUE, stringsAsFactors = FALSE)
# clinical1$PATIENT_ID <- tolower(substring(clinical1$PATIENT_ID,1,12)) #get just patient ID
# row.names(clinical1) <- clinical1$PATIENT_ID
#
# length(unique (intersect (useClini$patient.bcr_patient_barcode, clinical1$PATIENT_ID)))
# mergedClini <- merge(x = useClini, y = clinical1, by.x = "patient.bcr_patient_barcode", by.y = "PATIENT_ID")
#
# write.csv(as.data.frame(t(mergedClini)), "temp/mergedClini.csv")
#
# ####### verify repeated columns
# which(mergedClini$patient.gender != tolower(mergedClini$Gender) )
# #Do not use this patient
# mergedClini[471, c("patient.bcr_patient_barcode","SAMPLE_ID","patient.gender", "Gender") ]
#
# which(mergedClini$patient.age_at_initial_pathologic_diagnosis != mergedClini$AGE )
# #Do not use this patient
# mergedClini[296, c("patient.bcr_patient_barcode","SAMPLE_ID","patient.age_at_initial_pathologic_diagnosis",
#                    "AGE") ]
#
# which(mergedClini$patient.breast_carcinoma_estrogen_receptor_status != tolower(mergedClini$ER.Status) )
#
# which(mergedClini$patient.breast_carcinoma_progesterone_receptor_status != tolower(mergedClini$PR.Status) )
#
# tempStatus <-mergedClini$OS_STATUS
# tempStatus[which(tempStatus == "LIVING")] = "alive"
# tempStatus[which(tempStatus == "DECEASED")] = "dead"
# which(mergedClini$patient.vital_status != tempStatus )
#
#
# ###### remove patients ------
# mergedClini <- mergedClini[-c(296,471),]
#
# ##### remove repeated and non useful columns
# merge_use_cols<-read.table(file="descriptionData/merge_useful.csv", sep = ",",
#                            header = TRUE, stringsAsFactors = FALSE, na.strings = "nasrting", fill = TRUE)
# name_useful_mer_col <- merge_use_cols[which(merge_use_cols$repeated == '' ),1]
#
# mergedClini <- mergedClini[,name_useful_mer_col]
#
#
# #### describe remaning variables
# nominalCols <- merge_use_cols$columns[which(merge_use_cols$type == "nominal")]
# cap1 <- "Patient characteristics: nominal variables. Merged firebrowse and cBio"
# write(tableNominal(vars = mergedClini[,nominalCols], cap = cap1, vertical = TRUE, lab = "tab: nominal1",
#                    longtable = TRUE), file = "temp/1-1.2-NominalVarsFinal.tex")
#
#
# continousCols <- merge_use_cols$columns[which(merge_use_cols$type == "continuos")]
# numClini <- mergedClini[,continousCols]
# indx <- sapply(numClini, is.factor)
# numClini[indx] <- lapply(numClini[indx], function(x) as.numeric(as.character(x)))
# cap1 <- "Patient characteristics: continuous Merged firebrowse and cBio"
# write(tableContinuous(vars = numClini, cap = cap1, vertical = TRUE, lab = "tab: continuous2",
#                       longtable = TRUE), file = "temp/1-1.2-ContinousVarsFinal.tex")
#
# ##TODO: MOVE OUTPUT MANUALLY TO FILE!!!
#
# ######## construct clinical tables -----
# merge_use_cols[15,c(2:4)]  <- c("no_useful","","")
# merge_use_cols[16,c(2:4)]  <- c("no_useful","","")
# tensorCols <- merge_use_cols$columns[which(merge_use_cols$for.tensor == "tensor")]
# tensorClinical <- mergedClini[,c("patient.bcr_patient_barcode", tensorCols)]
# save(tensorClinical, file="data/tensorClinical.RData")
#
# survCols <- merge_use_cols$columns[which(merge_use_cols$for.tensor == "survival")]
# survClinical <- mergedClini[,c("patient.bcr_patient_barcode", survCols)]
# save(survClinical, file="data/survClinical.RData")
#
# treatCols <- merge_use_cols$columns[which(merge_use_cols$for.tensor == "treatment")]
# treatClinical <- mergedClini[,c("patient.bcr_patient_barcode", treatCols)]
# save(treatClinical, file="data/treatClinical.RData")
#
# save.image(file = "temp/1-1.2.RData")
# #load(file = "temp/1-1.2.RData")
