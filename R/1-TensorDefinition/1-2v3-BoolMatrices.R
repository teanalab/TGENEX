# version 3 (5/20/2018)
#1) patients with fewer than 10 mutations were discarded initially, but after trimming, remaining patients have three or more mutations
#2) only important variables
#3) using cBio data
#4) using firebrowse and cBiolite data

######### start here ----
rm(list = ls())

#source
readFiles <- function(){
  dataFile <- "temp/1-2v3data.RData"
  #file.remove(dataFile)
  if( file.exists(dataFile) ){
    load(dataFile, env = globalenv())
  } else {
    load("data/tensorClinical.RData")
    clinical <- tensorClinical
    mutation <-read.table(file="data/brca4mCBioportal/data_mutations_extended.txt", sep = "\t", quote = "",
                           header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)

    #uuid = NO A GOOD IDENTIFIER
    #Tumor_Sample_Barcode  = GOOD IDENTIFIER, BUT WE NEED TO CUT
    #TCGA-06-5411-01A-01D-1696-08
    mutationBarcodes <- substring(tolower(mutation$Tumor_Sample_Barcode), 1, 12)
    mutation$mutationBarcodes <- mutationBarcodes

    patients <<- intersect(mutationBarcodes, clinical$patient.bcr_patient_barcode)

    #load my libraries
    libs<-c("Packages.R")
    libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
    sapply(libs, function(u) {source(u)})
    clinical <<- clinical
    mutation <<- mutation
    save.image(dataFile)
  }
}

#source
createBoolMutation <- function(boolOrBin="bool"){
  dataFile = "data/boolMutation.RData"
  if( file.exists(dataFile) ){
    load(dataFile, env = globalenv())
  } else {
    genes <- unique(as.character(mutation$Hugo_Symbol))

    if(boolOrBin == "bool")
    {
      boolMutation <- data.frame(matrix(FALSE,length(patients),length(genes)),
                                 row.names = patients)

      trueVal = TRUE
    } else{
      boolMutation <- data.frame(matrix(0,length(patients),length(genes)),
                                 row.names = patients)
      trueVal = 1
    }
    names(boolMutation) <- genes



    for(i in seq_len(dim(mutation)[1])){
      Apatient <- as.character(mutation$mutationBarcodes[i])
      if(Apatient %in% patients){
        #only non-silent mutations
        TMuta <- unique(mutation$Variant_Classification[i])
        if(TMuta != "Silent"){
          APgene <- mutation$Hugo_Symbol[i]
          boolMutation[Apatient,APgene] <- trueVal
        }
      }
    }

    ####
    #plot distributions
    #1 Number of mutations per patient
    mut4patient <- apply(boolMutation,1,sum)
    range(mut4patient)

    # 8
    # 5 = six patients
    # 2 no patients
    #NOTE: Patients to remove
    length(which(mut4patient<=10))
    boolMutation <- boolMutation[-which(mut4patient<=10),]
    #4 patients were removed with less than 10 mutations
    #tcga-a2-a0es tcga-a8-a08c tcga-bh-a0dv tcga-bh-a0h5
    #38          113          355          372

    #after deletion
    mut4patient <- apply(boolMutation,1,sum)
    range(mut4patient)

    pdf(file = "temp/1-2v3-histomut4pati.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
    hist(mut4patient, main = "Histogram mutations per patient (mutation count)" ,
         ylab = "number of patients" ,
         xlab = "number of genes with mutations",
         breaks = 200)
    dev.off()

    #2 Number of patients with mutation x
    pat4genes <- apply(boolMutation,2,sum)
    range(pat4genes)

    #NOTE: Genes to remove
    length(genes)
    length(which(pat4genes<=0))
    boolMutation <- boolMutation[,-which(pat4genes<=0)]
    genes <- unique(names(boolMutation))
    length(genes)

    #again after delete 0s
    pat4genes <- apply(boolMutation,2,sum)
    range(pat4genes)


    #NOTE: Genes to remove with 1
    length(genes)
    length(which(pat4genes==1))
    boolMutation <- boolMutation[,-which(pat4genes==1)]
    genes <- unique(names(boolMutation))
    length(genes)

    #again after delete 0s
    pat4genes <- apply(boolMutation,2,sum)
    range(pat4genes)

    #NOTE: Genes to remove with 2
    length(genes)
    length(which(pat4genes==2))
    boolMutation <- boolMutation[,-which(pat4genes==2)]
    genes <- unique(names(boolMutation))
    length(genes)

    #again after delete 0s
    pat4genes <- apply(boolMutation,2,sum)
    range(pat4genes)

    # #NOTE: Genes to remove with 3
    # length(genes)
    # length(which(pat4genes==3))
    # boolMutation <- boolMutation[,-which(pat4genes==3)]
    # genes <- unique(names(boolMutation))
    # length(genes)
    #
    # #again after delete 0s
    # pat4genes <- apply(boolMutation,2,sum)
    # range(pat4genes)
    #
    # #NOTE: Genes to remove with 4
    # length(genes)
    # length(which(pat4genes==4))
    # boolMutation <- boolMutation[,-which(pat4genes==4)]
    # genes <- unique(names(boolMutation))
    # length(genes)
    #
    # #again after delete 4s
    # pat4genes <- apply(boolMutation,2,sum)
    # range(pat4genes)


    pdf(file = "temp/1-2v3-Fig4.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
    hist(pat4genes, main = "Histogram patient per gene (patient count)",
         ylab = "Number of genes",
         xlab = "Number of mutations", breaks=70)
    dev.off()


    #description boolean table
    genes <- names(boolMutation)
    patients <- row.names(boolMutation)

    length(genes)
    length(patients)
    dim(boolMutation)

    save(boolMutation,patients, file="data/boolMutation.RData")
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
