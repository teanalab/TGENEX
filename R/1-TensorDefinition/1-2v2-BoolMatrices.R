# version 2 (feb/3/2018)
#1) patients with fewer than 10 mutations were discarded
#2) only important variables
#3) using cBio data


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




##################

rm(list = ls())
load("temp/10-BoolMatrices_v2.RData")

#load my libraries
libs<-c("Packages.R")
libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
sapply(libs, function(u) {source(u)})

#required for the cox proportional hazard model
loadls("classInt")


#Data adjustment
readFiles <- function(){
  clinical1<-read.table(file="../../12Datasets/BRCA/brca4mCBioportal/data_clinical.txt", sep = "\t",
                       header = TRUE, stringsAsFactors = FALSE)

  mutation<-read.table(file="../../12Datasets/BRCA/brca4mCBioportal/data_mutations_extended.txt", sep = "\t", quote = "",
                       header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)

  #uuid = NO A GOOD IDENTIFIER
  #Tumor_Sample_Barcode  = GOOD IDENTIFIER, BUT WE NEED TO CUT
  #TCGA-06-5411-01A-01D-1696-08



  mutationBarcodes <- substring(mutation$Tumor_Sample_Barcode, 1, 12)
  mutationWithBC <- cbind(mutation,mutationBarcodes)

  #DATASETS DESCRIPTION

  dim(clinical)
  dim(mutation)
  #names(clinical)
  #names(mutation)
  mutationRange <- t(sapply(mutation, range))
  #write.csv(mutationRange,file="descriptionData/brca/4cBiomutationRange.csv")
  clinicalRange <- t(sapply(clinical, range))
  #write.csv(clinicalRange,file="descriptionData/brca/4cBioclinicalRange.csv")

  #to map stuff 4 cBio
  names(clinical)[2] = "Patient.ID"
  clinical$Patient.ID <- substring(clinical$Patient.ID,1,12)

  #number of patients
  length(intersect(mutationBarcodes, clinical$Patient.ID))
  patients <- intersect(mutationBarcodes, clinical$Patient.ID)
}

createBoolMutation <- function(){
  genes <- unique(as.character(mutation$Hugo_Symbol))

  boolMutation <- data.frame(matrix(FALSE,length(patients),length(genes)),
                             row.names = patients)
  names(boolMutation) <- genes

  for(i in seq_len(dim(mutation)[1])){
    Apatient <- as.character(mutationWithBC$mutationBarcodes[i])
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

  #NOTE: Patients to remove
  boolMutation <- boolMutation[-which(mut4patient<=10),]

  #after deletion
  mut4patient <- apply(boolMutation,1,sum)
  range(mut4patient)

  pdf(file = "temp/10-Fig1_Mutation4pati_v2.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  plot(mut4patient,main = "Number of mutations per patient",
       ylab = "number of genes with mutations",
       xlab = "patient")
  dev.off()

  pdf(file = "temp/10-Fig2_histomut4pati_v2.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  hist(mut4patient, main = "Histogram mutations per patient (mutation count)" ,
       ylab = "number of patients" ,
       xlab = "number of genes with mutations",
       breaks = 200)
  dev.off()

  #2 Number of patients with mutation x
  pat4genes <- apply(boolMutation,2,sum)
  range(pat4genes)

  #NOTE: Genes to remove
  length(which(pat4genes<=0))
  boolMutation <- boolMutation[,-which(pat4genes<=0)]

  #again after delete 0s
  pat4genes <- apply(boolMutation,2,sum)
  range(pat4genes)

  pdf(file = "temp/10-Fig3.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  plot(pat4genes,main = "Number of patients per mutation",
       ylab = "number of patients with mutation in gene x",
       xlab = "gene")
  dev.off()

  pdf(file = "temp/10-Fig4.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  hist(pat4genes, main = "Histogram patient per gene (patient count)",
       ylab = "frequency of number of patients with mutation in gene x",
       xlab = "gene", breaks=70)
  dev.off()


  #description boolean table
  genes <- names(boolMutation)

  patients <- row.names(boolMutation)

  length(genes)
  length(patients)
  dim(boolMutation)


}

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


#Dichomotization of columns with more than 2 categories
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


createBoolClinical <-function(){
  #0- get only the patients that will be part of the analysis
  clinical <- clinical[which(clinical$Patient.ID %in% patients),]
  row.names(clinical) <- clinical$Patient.ID
  AllColsClinical <- clinical

  ########
  # names(clinical) 30
  # NO INFO
  # [1] "SAMPLE_ID"                        "Patient.ID"
  # [29] "CANCER_TYPE" ("Cancer.Studies")                    "CANCER_TYPE_DETAILED"
  #
  # Important Info
  # [3] "Gender"                           "AGE" ("Diagnosis.Age")
  # [5] "ER.Status"                        "PR.Status"
  # [7] "HER2.Final.Status"                "Tumor"
  # [9] "Tumor..T1.Coded"                  "Node"
  # [11] "Node.Coded"                       "Metastasis"
  # [13] "Metastasis.Coded"                 "AJCC.Stage"
  # [15] "Converted.Stage"                  "Survival.Data.Form"
  # [17] "OS_STATUS"                        "OS_MONTHS"        ("Overall.Survival.Status")
  # [19] "PAM50.Subtype"
  #
  # Clusters
  # [20] "SigClust.Unsupervised.mRNA"
  # [21] "SigClust.Intrinsic.mRNA"          "miRNA.Clusters"
  # [23] "methylation.Clusters"             "RPPA.Clusters"
  # [25] "CN.Clusters"                      "Integrated.Clusters..with.PAM50."
  # [27] "Integrated.Clusters..no.exp."     "Integrated.Clusters..unsup.exp."

  #########

  clustersColsNames <- names(clinical)[c(20:28)]

  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  #remove columns with only one value
  OnlyOneValsC <- which(NumVals==1)
  clinical <- clinical [, -OnlyOneValsC]

  #remove columns with IDS
  OnlyOneValsC <- which(NumVals==503)
  clinical <- clinical [, -OnlyOneValsC]

  #remove all cluter cols
  clustNIndex <- which(names(clinical) %in% clustersColsNames)
  clinical <- clinical [, -clustNIndex]

  #17 variables remaning
  # > names(clinical)
  # [1] "Gender"             "AGE"                "ER.Status"          "PR.Status"
  # [5] "HER2.Final.Status"  "Tumor"              "Tumor..T1.Coded"    "Node"
  # [9] "Node.Coded"         "Metastasis"         "Metastasis.Coded"   "AJCC.Stage"
  # [13] "Converted.Stage"    "Survival.Data.Form" "OS_STATUS"          "OS_MONTHS"
  # [17] "PAM50.Subtype"

  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  #"Tumor" and "Tumor..T1.Coded" are almost same
  #"Node" and "Node.coded" are almost the same
  #"Metastasis" and "Metastasis.Coded" are exactly the same
  clinical <- clinical [, -c(7,9,11)]




  #Node and Node.Details should be merged (dependent)

  #1- run analysis - Plots
  NumVals <- analyzeClinicVariables()
  NumVals

  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals

  #2- get dichotomous columns

  #-----------------
  #lenght 2
  matDich <- dichoSeveralCols(which(NumVals==2))  #Gender Survival.Data.Form          OS_STATUS
  View(matDich)
  #delete redundant info, leave only one column
  matDich <- matDich[,-c(2,4,6)]
  #to Bool Matrix
  matDich <- as.data.frame(lapply(matDich, as.logical))

  #-----------------
  #columns with 3
  cols3 <- which(NumVals==3)  #Metastasis
  #apply(clinical[,cols3], 2, unique)
  unique(clinical[,cols3])
  #na to character
  clinical[which(is.na(clinical[,cols3])),cols3] = "NA"

  matD3 <- dichoSeveralCols(cols3)
  View(matD3)

  #delete last columns coz are redundant, but not this time
  names(matD3)
  ToDel <- seq(3,length(names(matD3)),by = 3) #c(1,4,7,10)
  matD3 <- matD3[,-ToDel]

  #assing nice names
  names(matD3) <- c("Metastasis.M0", "Metastasis.M1")

  #to Bool Matrix
  matD3 <- as.data.frame(lapply(matD3, as.logical))


  # conc1 <- cbind(matDich,matD3,matD4)
  conc1 <- cbind(matDich,matD3)


  #4- add columns with clear categories, but some cells with [not available]
  #one by one
  ListOfMat <- list()


  ####
  # [3] ER.Status
  C<-3
  i<-1
  unique(clinical[,C])
  #na to character
  clinical[which(clinical[,C]=="Performed but Not Available"),C] = "NA"
  mat9 <- deleteBlanks(C)
  #View(mat9)
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))

  #[4] PR.Status
  C<-4
  i <- i+1
  unique(clinical[,C])
  #na to character
  clinical[which(clinical[,C]=="Performed but Not Available"),C] = "NA"
  clinical[which(clinical[,C]=="Not Performed"),C] = "Not_Performed"
  mat9 <- deleteBlanks(C)
  #View(mat9)
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))

  ###
  #[5] HER2.Final.Status
  C<-5
  #i <- i+1
  unique(clinical[,C])
  #na to character
  clinical[which(is.na(clinical[,C])),C] = "NA"
  clinical[which(clinical[,C]=="Not Available"),C] = "NA"
  mat9 <- deleteBlanks(C)
  #View(mat9)
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))

  ##
  # [6] Tumor
  C<-C+1
  i <- i+1
  unique(clinical[,C])
  #na to character
  clinical[which(is.na(clinical[,C])),C] = "NA"
  #clinical[which(clinical[,C]=="Not Available"),C] = "NA"
  mat9 <- dichoSeveralCols(C)
  View(mat9)
  names(mat9) <- c("Tumor.NA","Tumor.T1","Tumor.T2","Tumor.T3","Tumor.T4","Tumor.TX")
  mat9 <- mat9[,-1]
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))


  ##
  # [7] Node
  C<-C+1
  i <- i+1
  unique(clinical[,C])
  #na to character
  clinical[which(is.na(clinical[,C])),C] = "NA"
  #clinical[which(clinical[,C]=="Not Available"),C] = "NA"
  mat9 <- dichoSeveralCols(C)
  View(mat9)
  names(mat9) <- paste("Node.N",c(0,1,2,3,'A'),sep='')
  mat9 <- mat9[,-5]
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))

  ##
  #[9] AJCC.Stage
  C<-9
  i <- i+1
  unique(clinical[,C])
  #na to character
  clinical[which(is.na(clinical[,C])),C] = "NA"
  clinical[which(clinical[,C]=="[Not Available]"),C] = "NA"
  mat9 <- deleteBlanks(C)
  View(mat9)
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))


  ##
  #[10] Converted.Stage
  C<-C+1
  i <- i+1
  unique(clinical[,C])
  #na to character
  clinical[which(is.na(clinical[,C])),C] = "NA"
  #clinical[which(clinical[,C]=="[Not Available]"),C] = "NA"
  mat9 <- deleteBlanks(C)
  View(mat9)
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))


  ##
  #[14] PAM50.Subtype
  C<- 14
  i <- i+1
  unique(clinical[,C])
  mat9 <- dichoSeveralCols(C)
  View(mat9)
  names(mat9) <- paste("PAM50.",c("Basal-like","HER2-enriched","Luminal-A","Luminal-B","Normal-like"),sep='')
  ListOfMat[[i]] <- as.data.frame(lapply(mat9, as.logical))


  #check all have the same dim
  sapply(ListOfMat,dim)

  #concatenate
  conc2 <- ListOfMat[[1]]
  for(i in seq_len(length(ListOfMat)-1))
  {
    conc2 <- data.frame(conc2,ListOfMat[[i+1]], check.names = TRUE, stringsAsFactors = FALSE)
  }
  #verify number of columns
  sum(sapply(ListOfMat,dim)[2,])
  dim(conc2)


  ############
  #5- add other columns, one by one

  #---------
  #using Sturges
  ##2	Diagnosis.Age
  C<-2
  aCOl <- clinical[,C]
  range(aCOl)
  aMat <- classesSturges(C)
  View(aMat)

  #to Logical
  aMati <- as.data.frame(lapply(aMat, as.logical))
  names(aMati) <- names(aMat)

  ##13	OS_MONTHS
  C<-13
  aCOl <- clinical[,C]
  range(aCOl)
  aMat2 <- classesSturges(C)
  View(aMat2)
  #to Logical
  aMati2 <- as.data.frame(lapply(aMat2, as.logical))
  names(aMati2) <- names(aMat2)

  # #concatenate step 5
  conc3 <- cbind(aMati,aMati2)

  ############

  #6- concatenate all cols
  boolClinical <- cbind(conc1, conc2, conc3)
  dim(boolClinical)
  row.names(boolClinical) <- row.names(clinical)

  ##############
  #PLOTS

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


}


boolClinicalN <- data.frame(patientID = row.names(boolClinical), boolClinical)
write.csv2(boolClinicalN, file = "temp/10-boolC.csv",row.names = FALSE)

boolMutationN <- data.frame(patientID = row.names(boolMutation), boolMutation)
write.csv2(boolMutationN, file = "temp/10-boolM.csv",row.names = FALSE)

#save.image("temp/10-BoolMatrices_v2.RData")
load("temp/10-BoolMatrices_v2.RData")
