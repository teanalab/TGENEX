# version 3 (5/20/2018)

#1) patients with fewer than 10 mutations were discarded
#2) only important variables
#3) using cBio data

#4) using firebrowse and cBiolite data




#################
# >   length(genes)
# [1] 11996
# >   length(patients)
# [1] 503
# >   dim(boolMutation)
# [1]   503 11996

#original table with 29 variables
#2 are ID variables
#2 are cancer type = all bc
#9 are clusters
#3 vars are repeated:
## "Tumor" and "Tumor..T1.Coded" are almost same
## "Node" and "Node.coded" are almost the same
## "Metastasis" and "Metastasis.Coded" are exactly the same
#
#14 remaining variables (Total included)
# > names(NumVals)
# [1] "Gender"             "AGE"                "ER.Status"          "PR.Status"
# [5] "HER2.Final.Status"  "Tumor"              "Node"               "Metastasis"
# [9] "AJCC.Stage"         "Converted.Stage"    "Survival.Data.Form" "OS_STATUS"
# [13] "OS_MONTHS"          "PAM50.Subtype"
#12 categorical
#2 vars are continuous (Age and OS_Months)
#dicotomized using  Sturges


######### start here ----

rm(list = ls())
#Data adjustment
readFiles <- function(){
  dataFile <- "temp/1-2v3data.RData"
  #file.remove(dataFile)
  if( file.exists(dataFile) ){
    load(dataFile, env = globalenv())
  } else {
    load("data/tensorClinical.RData")
    clinical <- tensorClinical
    mutation <-read.table(file="../../12Datasets/BRCA/brca4mCBioportal/data_mutations_extended.txt", sep = "\t", quote = "",
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

#run one by one
createBoolMutation <- function(){
  dataFile = "data/boolMutation.RData"
  if( file.exists(dataFile) ){
    load(dataFile, env = globalenv())
  } else {
    genes <- unique(as.character(mutation$Hugo_Symbol))
    boolMutation <- data.frame(matrix(FALSE,length(patients),length(genes)),
                               row.names = patients)
    names(boolMutation) <- genes

    for(i in seq_len(dim(mutation)[1])){
      Apatient <- as.character(mutation$mutationBarcodes[i])
      if(Apatient %in% patients){
        #only non-silent mutations
        TMuta <- unique(mutation$Variant_Classification[i])
        if(TMuta != "Silent"){
          APgene <- mutation$Hugo_Symbol[i]
          boolMutation[Apatient,APgene] <- TRUE
        }
      }
    }

    ####
    #plot distributions
    #1 Number of mutations per patient
    mut4patient <- apply(boolMutation,1,sum)
    range(mut4patient)

    # 8
    #
    # 5 = six patients
    #
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
#DATASETS DESCRIPTION
dim(clinical)
dim(boolMutation)
#number of patients
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
  dichoCols <- as.data.frame(clinical[, cols])
  dichoCols <- lapply(dichoCols, factor)

  ##Recode categories to columns
  matDich <- model.matrix(~ . + 0, data=dichoCols,
                          contrasts.arg = lapply(dichoCols, contrasts, contrasts=FALSE))
  matDich <- as.data.frame(matDich)
  matDich <- as.data.frame(lapply(matDich, as.logical))
  row.names(matDich) <- row.names(clinical)
  matDich
}


#Dichomotization of discrete columns with several categories
classesSturges <- function(colId)
{
  require(classInt)
  aInt <- classIntervals(as.numeric(clinical[,colId]), style = "jenks")
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

  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  #NumVals

  #remove columns with only one value
  #done manually with excel
  # OnlyOneValsC <- which(NumVals==1)
  # clinical <- clinical [, -OnlyOneValsC]

  #remove columns with IDS
  #done manually with excel
  # OnlyOneValsC <- which(NumVals==503)
  # clinical <- clinical [, -OnlyOneValsC]

  #remove all cluter cols
  #done manually with excel
  # clustNIndex <- which(names(clinical) %in% clustersColsNames)
  # clinical <- clinical [, -clustNIndex]


  #"Tumor" and "Tumor..T1.Coded" are almost same
  #"Node" and "Node.coded" are almost the same
  #"Metastasis" and "Metastasis.Coded" are exactly the same
  #done manually with excel
  # clinical <- clinical [, -c(7,9,11)]

  #Node and Node.Details should be merged (dependent)

  #1- run analysis - Plots
  NumVals <- analyzeClinicVariables()
  NumVals

  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  #2- get dichotomous columns
  #3 levels columns-----------------
  cols3 <- which(NumVals==3)  #patient.ethnicity
  #Node.Coded  Metastasis.Coded

  #patient.ethnicity
  #142 NA values mapped as non-hispanic because of the country of sample (germany, NA, poland, united states, vietnam)
  patient.ethnicity = clinical[,"patient.ethnicity"]
  patient.ethnicity [which(is.na(patient.ethnicity)) ] = "not hispanic or latino"
  table(patient.ethnicity)
  clinical[,"patient.ethnicity"] = patient.ethnicity

  #Node.Coded
  #1 NA, got rid of patient
  Node.Coded = clinical[,"Node.Coded"]
  toDelete <- which(is.na(Node.Coded))
  clinical <- clinical[-toDelete,]
  patients <- patients[-toDelete]
  boolMutation <- boolMutation[-toDelete,]
  table(Node.Coded)


  #Metastasis.Coded
  #1 NA, got rid of patient
  Metastasis.Coded = clinical[,"Metastasis.Coded"]
  toDelete <- which(is.na(Metastasis.Coded))
  clinical <- clinical[-toDelete,]
  patients <- patients[-toDelete]
  boolMutation <- boolMutation[-toDelete,]
  table(Metastasis.Coded)


  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  #2 level columns -----------------
  matDich2 <- dichoSeveralCols(which(NumVals==2))
  #View(matDich2)
  #delete redundant info, leave only one column
  #reduCols <- seq(1,(2*length(which(NumVals==2))),by=2)
  #select manually for true being higher risk
  reduCols <- c(1,4,5,7,9)
  matDich2 <- matDich2[,-reduCols]

  #4 level columns-----------------
  cols4 <- which(NumVals==4)
  ncols <- names(cols4)
  # [1] "patient.breast_carcinoma_estrogen_receptor_status"
  table(clinical[,ncols[1]])
  # [2] "patient.breast_carcinoma_progesterone_receptor_status"
  table(clinical[,ncols[2]])
  # [3] "HER2.Final.Status"
  table(clinical[,ncols[3]])
  # [4] "Node"
  table(clinical[,ncols[4]])


  # [1] "patient.breast_carcinoma_estrogen_receptor_status"
  #5 NA, got rid of patients
  patient.breast_carcinoma_estrogen_receptor_status = clinical[,"patient.breast_carcinoma_estrogen_receptor_status"]
  toDelete <- which(is.na(patient.breast_carcinoma_estrogen_receptor_status))
  clinical <- clinical[-toDelete,]
  patients <- patients[-toDelete]
  boolMutation <- boolMutation[-toDelete,]
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
  clinical <- clinical[-toDelete,]
  patients <- patients[-toDelete]
  boolMutation <- boolMutation[-toDelete,]
  HER2.Final.Status = clinical[,"HER2.Final.Status"]
  table(HER2.Final.Status)
  #apply(clinical[,cols3], 2, unique)

  #Do like it has 3 levels
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  matDich3 <- dichoSeveralCols(which(NumVals==3))

  #Do like it has 4 levels
  matDich4 <- dichoSeveralCols(which(NumVals==4))

  ## OTHER --
  #na to character
  #clinical[which(is.na(clinical[,cols3])),cols3] = "NA"

  #assing nice names
  #names(matD3) <- c("Metastasis.M0", "Metastasis.M1")


  #5 levels -----------------
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  cols5 <- which(NumVals==5)
  ncols <- names(cols5)
  # [1] "Tumor"
  table(clinical[,ncols[1]])

  #Do like it has 5 levels
  matDich5 <- dichoSeveralCols(ncols[1])

  #7 levels -----------------
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==7)
  ncols <- names(colsN)
  # [1] patient.histological_type"
  anyNA(clinical[,ncols[1]])
  names(table(clinical[,ncols[1]]))


  #1 NA, got rid of patient
  patient.histological_type = clinical[,"patient.histological_type"]
  toDelete <- which(is.na(patient.histological_type))
  clinical <- clinical[-toDelete,]
  patients <- patients[-toDelete]
  boolMutation <- boolMutation[-toDelete,]
  anyNA(clinical[,ncols[1]])
  names(table(clinical[,ncols[1]]))


  #Do like it has 6 levels
  matDich6 <- dichoSeveralCols(ncols[1])

  #8 levels -----------------
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==8)
  ncols <- names(colsN)
  names(table(clinical[,ncols[1]]))
  anyNA(clinical[,ncols[1]])

  #Do like it has N levels
  matDich8 <- dichoSeveralCols(ncols[1])

  #13 levels -----------------
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==13)
  ncols <- names(colsN)
  ncols
  # [1] "patient.tissue_source_site"
  names(table(clinical[,ncols[1]]))
  anyNA(clinical[,ncols[1]])
  # [2] "AJCC.Stage"
  table(clinical[,ncols[2]])
  anyNA(clinical[,ncols[2]])

  # [2] "AJCC.Stage"
  AJCC.Stage = clinical[,"AJCC.Stage"]
  toDelete <- which(AJCC.Stage == "[Not Available]")
  clinical <- clinical[-toDelete,]
  patients <- patients[-toDelete]
  boolMutation <- boolMutation[-toDelete,]
  AJCC.Stage = clinical[,"AJCC.Stage"]
  table(AJCC.Stage)
  #apply(clinical[,cols3], 2, unique)

  #Do like it has N levels
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==12)
  ncols <- names(colsN)
  matDich12 <- dichoSeveralCols(ncols[1])

  #Do like it has N levels
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==13)
  ncols <- names(colsN)
  matDich13 <- dichoSeveralCols(ncols[1])




  ## OTHER --
  #na to character
  #clinical[which(is.na(clinical[,cols3])),cols3] = "NA"

  #assing nice names
  #names(matD3) <- c("Metastasis.M0", "Metastasis.M1")


  ### CONCATENATE ALL nominal
  matDich2 <- matDich2[patients,]
  matDich3 <- matDich3[patients,]
  matDich4 <- matDich4[patients,]
  matDich5 <- matDich5[patients,]
  matDich6 <- matDich6[patients,]
  matDich8 <- matDich8[patients,]
  matDich12 <- matDich12[patients,]
  matDich13 <- matDich13[patients,]

  conc1 <- cbind(matDich2, matDich3, matDich4, matDich5, matDich6,
                 matDich8, matDich12, matDich13)

  #5- add other columns, one by one --------

  #using Sturges
  ##2	Diagnosis.Age
  C <- "patient.age_at_initial_pathologic_diagnosis"
  aCOl <- as.numeric(clinical[,C])
  range(aCOl)
  aMat <- classesSturges(C)
  View(aMat)

  #to Logical
  aMati <- as.data.frame(lapply(aMat, as.logical))
  names(aMati) <- names(aMat)
  row.names(aMati) <- row.names(clinical)


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
  boolClinical <- boolClinical[,-which(pat4cli<=0)]

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

  save(boolClinical,file= "data/boolClinical.RData")
}


# continue here --------

save(boolClinical,file= "data/boolClinical.RData")
save(boolMutation,patients, file="data/boolMutation.RData")

#load(file= "data/boolClinical.RData")
#load(file="data/boolMutation.RData")

write.csv2(boolClinical, file = "temp/10-boolC.csv",row.names = FALSE)
write.csv2(boolMutation, file = "temp/10-boolM.csv",row.names = FALSE)

binaClinical <- as.data.frame(lapply(boolClinical, as.numeric), stringsAsFactors = FALSE)
binaMutation <- as.data.frame(lapply(boolMutation, as.numeric), stringsAsFactors = FALSE)
write.csv2(binaClinical, file = "temp/10-binaC.csv",row.names = FALSE,
           col.names = FALSE)
write.csv2(t(binaMutation), file = "temp/10-binaM.csv",row.names = FALSE)



save.image("temp/1-2v3.RData")
# load("temp/10-BoolMatrices_v2.RData")
