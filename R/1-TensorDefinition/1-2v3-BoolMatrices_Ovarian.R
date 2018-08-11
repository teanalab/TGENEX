# version 3 (8/8/2018)
#1) patients with fewer than 10 mutations were discarded initially, but after trimming, remaining patients have three or more mutations
#2) only important variables
#3) using cBio data
#4) using firebrowse and cBiolite data

##### OVARIAN ##########
rm(list = ls())


######### mutation ----
source('~/CLIGEN_tgit/R/1-TensorDefinition/createBoolMutation.R')

boolMutation <- createBoolMutation(fileName = "/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_mutations_extended.txt")
#description boolean table
genes <- names(boolMutation)
patients <- row.names(boolMutation)
length(genes)
length(patients)
dim(boolMutation)
save(boolMutation, patients, file="temp/boolMutation_OVA.RData")

#write.csv2(boolMutation, file = "temp/10-boolM.csv",row.names = FALSE)
binaMutation <- as.data.frame(lapply(boolMutation, as.numeric), stringsAsFactors = FALSE)
write.csv2(t(binaMutation), file = "temp/1-binaM_OVA.csv",row.names = FALSE)


######### clinical cBioportal ----
clinical<-read.table(file="/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_clinical_patient.txt", sep = "\t",
                     header = FALSE, stringsAsFactors = FALSE,  fill = TRUE)
names(clinical) <- clinical[1,]
clinical <- clinical[-1,]

clinicalSample<-read.table(file="/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_clinical_sample.txt", sep = "\t",
                           header = FALSE, stringsAsFactors = FALSE,  fill = TRUE)
names(clinicalSample) <- clinicalSample[1,]
clinicalSample <- clinicalSample[-1,]
clinical <- merge(clinical,clinicalSample)

clinicalPatients <- substring(tolower(clinical$PATIENT_ID), 1, 12)
row.names(clinical) <- clinicalPatients
patients <- intersect(patients, clinicalPatients)
clinical <- clinical[patients,]


write.csv2(clinical, file = "temp/1-clinical_OVA.csv",row.names = FALSE)
#write.csv2(boolClinical, file = "temp/10-boolC.csv",row.names = FALSE)
#binaClinical <- as.data.frame(lapply(boolClinical, as.numeric), stringsAsFactors = FALSE)
#write.csv2(binaClinical, file = "temp/10-binaC.csv",row.names = FALSE)



#source
readFiles <- function(){
  #dataFile <- "temp/1-2v3data.RData"
  #file.remove(dataFile)
  if( file.exists(dataFile) ){
    load(dataFile, env = globalenv())
  } else {
    #load my libraries
    libs<-c("Packages.R")
    libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
    sapply(libs, function(u) {source(u)})
  }
}



######### start here ----
readFiles()
createBoolMutation()

row.names(clinical) <- clinical$patient.bcr_patient_barcode
clinical <- clinical[patients,]
#DATASETS DIMENSIONS
dim(clinical)
dim(boolMutation)
length(patients)

###### source these functions ----

#plots for clinical variables
checkCol <- function(column, main=""){
  barplot(table(column), main=main,
          legend.text = c("> Num Vars: ", length(unique(column)),
                          "> VALUES: " ,unique(column) ),
          args.legend = list(x ='topright', bty='n'))
}

#plots for clinical variables
analyzeClinicVariables <- function(){
  pdf(file = "temp/10-Fig5_ClinicalVariablesRaw.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  NumVals <- c()
  for(i in seq_along(names(clinical))){
    checkCol(clinical[,i],i)
    NumVals <- c(NumVals, length(unique(clinical[,i])) )
  }
  dev.off()

  # output the number of values per each column
  # to see if the variable is dichotomous or not
  NumVals
}


#1.1.6.	Dichotomization of clinical variables with more than 2 categories
#and blank values
deleteBlanks <- function(columnId)
{
  matCols <- as.data.frame(clinical[, columnId])
  names(matCols) <- names(clinical)[columnId]
  matCols[is.na(matCols)]<-''
  matCols <- lapply(matCols, factor)

  ##Recode categories to columns
  mati <- model.matrix(~ . + 0, data=matCols,
                       contrasts.arg = lapply(matCols, contrasts, contrasts=FALSE))
  mati <- as.data.frame(mati)

  if(dim(mati)[1] != dim(clinical)[1])
    stop("ERROR. rows lost")

  #delete blank cols = col with shorter name
  nchars <- sapply(names(mati),nchar)
  ToDel <- which(nchars == min(nchars))
  if(length(ToDel) == 0){
    stop("No col to delete")
  } else if(length(ToDel) > 1){
    stop("more than one col to delete")
  } else {
    mati <- mati[,-ToDel]
  }
  mati
}


#Dichomotization of columns with more than 2 categories
#and [Not available] values
dichomoNotAvail <- function(columnId)
{
  matCols <- as.data.frame(clinical[, columnId])
  names(matCols) <- names(clinical)[columnId]
  matCols <- lapply(matCols, factor)

  ##Recode categories to columns
  mati <- model.matrix(~ . + 0, data=matCols,
                       contrasts.arg = lapply(matCols, contrasts, contrasts=FALSE))
  mati <- as.data.frame(mati)

  #delete [Not Available] cols
  ToDel <- grep("\\[", names(mati))
  if(length(ToDel) !=0){
    mati <- mati[,-ToDel]
  }
  mati
}

#Dichomotization of columns with more than 2 categories
dichoSeveralCols <-function(cols)
{
  dichoCols <- as.data.frame(clinical[patients, cols])
  dichoCols <- lapply(dichoCols, factor)

  ##Recode categories to columns
  matDich <- model.matrix(~ . + 0, data=dichoCols,
                          contrasts.arg = lapply(dichoCols, contrasts, contrasts=FALSE))
  matDich <- as.data.frame(matDich)
  matDich <- as.data.frame(lapply(matDich, as.logical))
  row.names(matDich) <- patients
  matDich
}


#Dichomotization of discrete columns with several categories
classesSturges <- function(colId)
{
  require(classInt)
  aInt <- classIntervals(as.numeric(as.character(clinical[patients,colId])), style = "jenks")
  aCategory <- findCols(aInt)

  ##Recode categories to columns
  boolClinical <- data.frame(aCategory)
  boolClinical <- lapply(boolClinical, factor)
  matA <- model.matrix(~ . + 0, data=boolClinical,
                       contrasts.arg = lapply(boolClinical, contrasts, contrasts=FALSE))
  matA <- as.data.frame(matA)

  interNames <-  print(aInt)
  interNames <- names(interNames)

  names(matA) <- paste(names(clinical)[colId],interNames, sep='')

  checkCol(aCategory, main = names(clinical)[colId])
  matA
}


###### run clinical line by line -----
createBoolClinical <-function(){
  #0- get only the patients that will be part of the analysis
  AllColsClinical <- clinical
  #remove cols after meeting (5/25/2018)
  clinical <- clinical[,c("patient.bcr_patient_barcode","patient.age_at_initial_pathologic_diagnosis",
                          "patient.breast_carcinoma_estrogen_receptor_status", "patient.breast_carcinoma_progesterone_receptor_status",
                          "patient.ethnicity", "patient.gender",
                          "patient.histological_type", "HER2.Final.Status",
                          "Tumor",
                          "Node.Coded", "Metastasis.Coded",
                          "Converted.Stage")]
  #transform tumor to tumor.coded after meeting (5/25/2018)
  Tumor.Coded <- clinical$Tumor
  Tumor.Coded[which(Tumor.Coded != "T1")] <- "T_other"
  Tumor.Coded -> clinical$Tumor
  names(clinical) = c("patient.bcr_patient_barcode","patient.age_at_initial_pathologic_diagnosis",
                      "patient.breast_carcinoma_estrogen_receptor_status", "patient.breast_carcinoma_progesterone_receptor_status",
                      "patient.ethnicity", "patient.gender",
                      "patient.histological_type", "HER2.Final.Status",
                      "Tumor.Coded",
                      "Node.Coded", "Metastasis.Coded",
                      "Converted.Stage")

  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  NumVals <- analyzeClinicVariables()
  NumVals

  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  #2- get dichotomous columns
  #3-level columns-----------------
  cols3 <- which(NumVals==3)  #patient.ethnicity
  #Node.Coded  Metastasis.Coded

  #patient.ethnicity
  #142 NA values mapped as non-hispanic because of the country of sample (germany, NA, poland, united states, vietnam)
  patient.ethnicity = clinical[,"patient.ethnicity"]
  patient.ethnicity [which(is.na(patient.ethnicity)) ] = "not hispanic or latino"
  table(patient.ethnicity, useNA = "ifany" )
  anyNA(patient.ethnicity)
  clinical[,"patient.ethnicity"] = patient.ethnicity

  #Node.Coded
  #1 NA, got rid of patient
  Node.Coded = clinical[,"Node.Coded"]
  toDelete <- which(is.na(Node.Coded))

  deletepatients <- function(toDelete){
    if(length(toDelete) != 0)
    {
      clinical <- clinical[-toDelete,]
      patients <- patients[-toDelete]
      boolMutation <- boolMutation[-toDelete,]
    }
  }

  deletepatients(toDelete)
  table(clinical$Node.Coded,useNA = "ifany" )


  #Metastasis.Coded
  #1 NA, got rid of patient
  Metastasis.Coded = clinical[,"Metastasis.Coded"]
  toDelete <- which(is.na(Metastasis.Coded))
  deletepatients(toDelete)
  table(Metastasis.Coded)


  #Tumor.Coded
  Tumor.Coded = clinical[,"Tumor.Coded"]
  table(Tumor.Coded)



  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  #2 level columns -----------------
  matDich2 <- dichoSeveralCols(cols=which(NumVals==2))
  #View(matDich2)
  #delete redundant info, leave only one column
  #reduCols <- seq(1,(2*length(which(NumVals==2))),by=2)
  #select manually for true being higher risk
  reduCols <- c(2,4,5,7,9,11)
  matDich2 <- matDich2[,-reduCols]

  #4 level columns-----------------
  cols4 <- which(NumVals==4)
  # ncols <- names(cols4)
  # # [1] "patient.breast_carcinoma_estrogen_receptor_status"
  # table(clinical[,ncols[1]])
  # # [2] "patient.breast_carcinoma_progesterone_receptor_status"
  # table(clinical[,ncols[2]])
  # # [3] "HER2.Final.Status"
  # table(clinical[,ncols[3]])


  # [1] "patient.breast_carcinoma_estrogen_receptor_status"
  #5 NA, got rid of patients
  patient.breast_carcinoma_estrogen_receptor_status = clinical[,"patient.breast_carcinoma_estrogen_receptor_status"]
  toDelete <- which(is.na(patient.breast_carcinoma_estrogen_receptor_status))
  deletepatients(toDelete)
  table(patient.breast_carcinoma_estrogen_receptor_status)
  #rename indeterminate to Equivocal
  patient.breast_carcinoma_estrogen_receptor_status = clinical[,"patient.breast_carcinoma_estrogen_receptor_status"]
  levels(patient.breast_carcinoma_estrogen_receptor_status)[levels(patient.breast_carcinoma_estrogen_receptor_status)=="indeterminate"] <- "Equivocal"
  table(patient.breast_carcinoma_estrogen_receptor_status)
  clinical[,"patient.breast_carcinoma_estrogen_receptor_status"] = patient.breast_carcinoma_estrogen_receptor_status

  # [2] "patient.breast_carcinoma_progesterone_receptor_status"
  #5 NA, got rid of patients
  patient.breast_carcinoma_progesterone_receptor_status = clinical[,"patient.breast_carcinoma_progesterone_receptor_status"]
  toDelete <- which(is.na(patient.breast_carcinoma_progesterone_receptor_status))

  # [3] "HER2.Final.Status"
  #9 NA, got rid of patients
  HER2.Final.Status = clinical[,"HER2.Final.Status"]
  toDelete <- which(HER2.Final.Status == "Not Available")
  deletepatients(toDelete)
  HER2.Final.Status = clinical[,"HER2.Final.Status"]
  table(HER2.Final.Status)
  #apply(clinical[,cols3], 2, unique)

  #Do like it has 3 levels
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  matDich3 <- dichoSeveralCols(which(NumVals==3))


  #7 levels -----------------
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==7)
  # ncols <- names(colsN)
  # # [1] patient.histological_type"
  # anyNA(clinical[,ncols[1]])
  # names(table(clinical[,ncols[1]]))


  #1 NA, got rid of patient
  patient.histological_type = clinical[,"patient.histological_type"]
  toDelete <- which(is.na(patient.histological_type))
  # deletepatients(toDelete)
  # anyNA(clinical[,ncols[1]])
  # names(table(clinical[,ncols[1]]))


  #Do like it has 6 levels
 # matDich6 <- dichoSeveralCols(ncols[1])

  #8 levels -----------------
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==8)
  ncols <- names(colsN)
  names(table(clinical[,ncols[1]], useNA = "ifany"))
  anyNA(clinical[,ncols[1]])

  #Do like it has N levels
  matDich8 <- dichoSeveralCols(ncols[1])

  ### CONCATENATE ALL nominal
  matDich2 <- matDich2[patients,]
  matDich3 <- matDich3[patients,]
  #matDich6 <- matDich6[patients,]
  matDich8 <- matDich8[patients,]

  conc1 <- cbind(matDich2, matDich3, #matDich6,
                 matDich8)

  #5- add other columns, one by one --------

  #using Sturges
  ##2	Diagnosis.Age
  C <- "patient.age_at_initial_pathologic_diagnosis"
  aCOl <- as.numeric(as.character(clinical[patients,C]))
  range(aCOl) #29 89
  aMat <- classesSturges(C)

#   style: jenks
#   one of 8,996,462,475 possible partitions of this variable into 10 classes
#   [3,10] (10,16] (16,22] (22,27] (27,32] (32,37] (37,42] (42,48] (48,54] (54,60]
#     25      46      52      63      56      74      59      32      32      17
#   [29,37] (37,43] (43,49] (49,54] (54,59] (59,64] (64,70] (70,76] (76,82] (82,89]
#     25      46      52      63      56      74      65      33      29      13

  aMatNew  <- classesSturges(C)

  #to Logical
  aMati <- as.data.frame(lapply(aMat, as.logical))
  names(aMati) <- names(aMat)
  row.names(aMati) <- patients

  aMati <- aMati[patients, ]

  #6- concatenate all cols
  boolClinical <- cbind(conc1, aMati)
  anyNA(boolClinical)
  dim(boolClinical)
  row.names(boolClinical) <- row.names(clinical)

  # plots ##########

  #description boolean table
  clinicalColunms <- names(boolClinical)
  dim(boolClinical)

  #plot distributions
  #1 Number of clinical variables per patient
  cli4patient <- apply(boolClinical,1,sum)
  range(cli4patient)

  pdf(file = "temp/10-Fig6.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  plot(cli4patient,main = "Number of clinical variables per patient",
       ylab = "number of true clinical variables",
       xlab = "patient")
  dev.off()

  pdf(file = "temp/10-Fig7.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  hist(cli4patient, main = "Histogram of clinical variables per patient (clinical variables count)",
       ylab = "frequency of true clinical variables",
       xlab = "patient")
  dev.off()

  #NOTE: Patients to remove (If remove patients, remove from mutations also)
  which(cli4patient<=0)

  #2 Number of patients with clinical variable x
  pat4cli <- apply(boolClinical,2,sum)

  #NOTE: Clinical to remove zero patients with the variable
  length(which(pat4cli<=0))
  #boolClinical <- boolClinical[,-which(pat4cli<=0)]

  #NOTE: Clinical to remove all patients with the variable
  #length(which(pat4cli==dim(boolClinical)[1]))

  pdf(file = "temp/10-Fig8.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  plot(pat4cli,main = "Number of patients per clinical variable",
       ylab = "number of patients with true clinical variable x",
       xlab = "clinical variable")
  dev.off()

  pdf(file = "temp/10-Fig9.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  hist(pat4cli, main = "Histogram patient per clinical variable (patient count)",
       ylab = "frequency of number of patients with true clinical variable x",
       xlab = "clinical variable", breaks = 40)
  dev.off()

  # save ##########
  save(boolClinical,file= "data/boolClinical.RData")
}


# continue here --------

save(boolClinical,file= "data/boolClinical.RData")
save(boolMutation,patients, file="data/boolMutation.RData")

#load(file= "data/boolClinical.RData")
#load(file="data/boolMutation.RData")

#write.csv2(boolClinical, file = "temp/10-boolC.csv",row.names = FALSE)
#write.csv2(boolMutation, file = "temp/10-boolM.csv",row.names = FALSE)

binaClinical <- as.data.frame(lapply(boolClinical, as.numeric), stringsAsFactors = FALSE)
binaMutation <- as.data.frame(lapply(boolMutation, as.numeric), stringsAsFactors = FALSE)
write.csv2(binaClinical, file = "temp/10-binaC.csv",row.names = FALSE)
write.csv2(t(binaMutation), file = "temp/10-binaM.csv",row.names = FALSE)

save(patients,file="data/patients.RData")

save.image("temp/1-2v3.RData")
# load("temp/10-BoolMatrices_v2.RData")
