# # version 1 (5/20/2017)
#
# #Data adjustment
# readFiles <- function(){
#   clinical<-read.table(file="../../12Datasets/BRCA/clinical.txt", sep = "\t",
#                        header = TRUE, stringsAsFactors = FALSE)
#
#   mutation<-read.table(file="../../12Datasets/BRCA/data_mutations_extended.txt", sep = "\t", quote = "",
#                        header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)
#
#   #uuid = NO A GOOD IDENTIFIER
#   #Tumor_Sample_Barcode  = GOOD IDENTIFIER, BUT WE NEED TO CUT
#   #TCGA-06-5411-01A-01D-1696-08
#
#   mutationBarcodes <- substring(mutation$Tumor_Sample_Barcode, 1, 12)
#   mutationWithBC <- cbind(mutation,mutationBarcodes)
#
#   #DATASETS DESCRIPTION
#
#   dim(clinical)
#   dim(mutation)
#   #names(clinical)
#   #names(mutation)
#   mutationRange <- t(sapply(mutation, range))
#   write.csv(mutationRange,file="descriptionData/brca/mutationRange.csv")
#   clinicalRange <- t(sapply(clinical, range))
#   write.csv(clinicalRange,file="descriptionData/brca/clinicalRange.csv")
#
#   #number of patients
#   length(intersect(mutationBarcodes, clinical$Patient.ID))
#   patients <- intersect(mutationBarcodes, clinical$Patient.ID)
# }
#
# createBoolMutation <- function(){
#   genes <- unique(as.character(mutation$Hugo_Symbol))
#
#   boolMutation <- data.frame(matrix(FALSE,length(patients),length(genes)),
#                              row.names = patients)
#   names(boolMutation) <- genes
#
#   for(i in seq_len(dim(mutation)[1])){
#     Apatient <- as.character(mutationWithBC$mutationBarcodes[i])
#     if(Apatient %in% patients){
#       #only non-silent mutations
#       TMuta <- unique(mutation$Variant_Classification[i])
#       if(TMuta != "Silent"){
#         APgene <- mutation$Hugo_Symbol[i]
#         boolMutation[Apatient,APgene] <- TRUE
#       }
#     }
#   }
#
#   ####
#   #plot distributions
#   #1 Number of mutations per patient
#   mut4patient <- apply(boolMutation,1,sum)
#   range(mut4patient)
#
#   #NOTE: Patients to remove
#   which(mut4patient<=0)
#   #boolMutation <- boolMutation[-which(mut4patient<=0),]
#
#   plot(mut4patient,main = "Number of mutations per patient",
#        ylab = "number of genes with mutations",
#        xlab = "patient")
#
#   hist(mut4patient, main = "Histogram mutations per patient (mutation count)",
#        ylab = "frequency of number of genes with mutations",
#        xlab = "patient", breaks = 70)
#
#
#
#   #2 Number of patients with mutation x
#   pat4genes <- apply(boolMutation,2,sum)
#
#   #NOTE: Genes to remove
#   length(which(pat4genes<=0))
#   boolMutation <- boolMutation[,-which(pat4genes<=0)]
#
#
#   plot(pat4genes,main = "Number of patients per mutation",
#        ylab = "number of patients with mutation in gene x",
#        xlab = "gene")
#
#   hist(pat4genes, main = "Histogram patient per gene (patient count)",
#        ylab = "frequency of number of patients with mutation in gene x",
#        xlab = "gene", breaks=70)
#
#   #description boolean table
#   genes <- names(boolMutation)
#
#   length(genes)
#   length(patients)
#   dim(boolMutation)
#
# }
#
# #plots for clinical variables
# checkCol <- function(column, main=""){
#   barplot(table(column), main=main,
#           legend.text = c("> Num Vars: ", length(unique(column)),
#                           "> VALUES: " ,unique(column) ),
#           args.legend = list(x ='topright', bty='n'))
# }
#
# #plots for clinical variables
# analyzeClinicVariables <- function(){
#   pdf( file="ClinicalVariablesRaw.pdf", onefile=TRUE)
#   opar <- par()
#   par(mfrow=c(1,1))
#   NumVals <- c()
#   for(i in seq_along(names(clinical))){
#     checkCol(clinical[,i],i)
#     NumVals <- c(NumVals, length(unique(clinical[,i])) )
#   }
#   par(opar)
#   graphics.off()
#
#   # output the number of values per each column
#   # to see if the variable is dichotomous or not
#   NumVals
# }
#
#
# #Dichomotization of columns with more than 2 categories
# #and blank values
# deleteBlanks <- function(columnId)
# {
#   matCols <- as.data.frame(clinical[, columnId])
#   names(matCols) <- names(clinical)[columnId]
#   matCols[is.na(matCols)]<-''
#   matCols <- lapply(matCols, factor)
#
#
#
#   ##Recode categories to columns
#   mati <- model.matrix(~ . + 0, data=matCols,
#                        contrasts.arg = lapply(matCols, contrasts, contrasts=FALSE))
#   mati <- as.data.frame(mati)
#
#   if(dim(mati)[1] != dim(clinical)[1])
#     stop("ERROR. rows lost")
#
#   #delete blank cols = col with shorter name
#   nchars <- sapply(names(mati),nchar)
#   ToDel <- which(nchars == min(nchars))
#   if(length(ToDel) == 0){
#     stop("No col to delete")
#   } else if(length(ToDel) > 1){
#     stop("more than one col to delete")
#   } else {
#     mati <- mati[,-ToDel]
#   }
#   mati
# }
#
#
#
# #Dichomotization of columns with more than 2 categories
# #and [Not available] values
# dichomoNotAvail <- function(columnId)
# {
#   matCols <- as.data.frame(clinical[, columnId])
#   names(matCols) <- names(clinical)[columnId]
#   matCols <- lapply(matCols, factor)
#
#   ##Recode categories to columns
#   mati <- model.matrix(~ . + 0, data=matCols,
#                        contrasts.arg = lapply(matCols, contrasts, contrasts=FALSE))
#   mati <- as.data.frame(mati)
#
#   #delete [Not Available] cols
#   ToDel <- grep("\\[", names(mati))
#   if(length(ToDel) !=0){
#     mati <- mati[,-ToDel]
#   }
#   mati
# }
#
# #Dichomotization of columns with more than 2 categories
# dichoSeveralCols <-function(cols)
# {
#   dichoCols <- as.data.frame(clinical[, cols])
#   dichoCols <- lapply(dichoCols, factor)
#
#   ##Recode categories to columns
#   matDich <- model.matrix(~ . + 0, data=dichoCols,
#                           contrasts.arg = lapply(dichoCols, contrasts, contrasts=FALSE))
#   matDich <- as.data.frame(matDich)
#   matDich
# }
#
# #Dichomotization of discrete columns with several categories
# classesSturges <- function(colId)
# {
#   require(classInt)
#   aInt <- classIntervals(as.numeric(clinical[,colId]), style = "jenks")
#   aCategory <- findCols(aInt)
#
#   ##Recode categories to columns
#   boolClinical <- data.frame(aCategory)
#   boolClinical <- lapply(boolClinical, factor)
#   matA <- model.matrix(~ . + 0, data=boolClinical,
#                        contrasts.arg = lapply(boolClinical, contrasts, contrasts=FALSE))
#   matA <- as.data.frame(matA)
#
#   interNames <-  print(aInt)
#   interNames <- names(interNames)
#
#   names(matA) <- paste(names(clinical)[colId],interNames, sep='')
#
#   checkCol(aCategory, main = names(clinical)[colId])
#   matA
# }
#
#
# createBoolClinical <-function(){
#   #0- get only the patients that will be part of the analysis
#   clinical <- clinical[which(clinical$Patient.ID %in% patients),]
#   row.names(clinical) <- clinical$Patient.ID
#
#   #1- run analysis
#   NumVals <- analyzeClinicVariables()
#
#   #2- Check pdf for categories that can be useful
#   #and document in excel file
#
#   #3- get dichotomous columns
#
#   #-----------------
#   #lenght 2
#   matDich <- dichoSeveralCols(c(20,23))
#   View(matDich)
#   #delete redundant info, leave only one column
#   matDich <- matDich[,-c(2,4)]
#   #to Bool Matrix
#   matDich <- as.data.frame(lapply(matDich, as.logical))
#
#
#
#
#   #-----------------
#   #columns with 3
#   cols3 <- which(NumVals==3)
#   apply(clinical[,cols3], 2, unique)
#
#   #after checking that these columns are not dichotomous
#   #cols3 <- cols3[-9] # icd_o_3_site
#   #cols3 <- cols3[-8] # repited with icd_o_3_site
#   #cols3 <- cols3[-7] # histological_type
#
#   matD3 <- dichoSeveralCols(cols3)
#
#
#   #delete [Not Available] cols
#   #ToDel <- grep("\\[", names(matD3))
#   #matD3 <- matD3[,-ToDel]
#
#   #delete blank cols
#   names(matD3)
#   ToDel <- seq(1,length(names(matD3)),by = 3) #c(1,4,7,10)
#   matD3 <- matD3[,-ToDel]
#
#   #to Bool Matrix
#   matD3 <- as.data.frame(lapply(matD3, as.logical))
#
#
#   #-----------------
#   #columns with 4
#   # cols4 <- which(NumVals==4)
#   # apply(clinical[,cols4], 2, unique)
#
#   #if only one column skip to next step
#
#   # #after checking that these columns are not dichotomous
#   # cols4 <- cols4[-3] #race
#   #
#   # matD4 <- dichoSeveralCols(cols4)
#   #
#   # #delete [Not Available] cols
#   # ToDel <- grep("\\[", names(matD4))
#   # matD4 <- matD4[,-ToDel]
#   #
#   # #to Bool Matrix
#   # matD4 <- as.data.frame(lapply(matD4, as.logical))
#
#   # #concatenate
#   # conc1 <- cbind(matDich,matD3,matD4)
#   conc1 <- cbind(matDich,matD3)
#
#   ############
#
#   #4- add columns with clear categories, but some cells with [not available]
#   #one by one
#   ListOfMat <- list()
#   ##9- HER2.Status
#   mat9 <- deleteBlanks(9)
#   ListOfMat[[1]] <- as.data.frame(lapply(mat9, as.logical))
#   ##17	Node
#   x <- deleteBlanks(17)
#   #names(x)
#   ListOfMat[[2]] <- as.data.frame(lapply(x, as.logical))
#   ##21	PAM50.subtype
#   x <- deleteBlanks(21)
#   ListOfMat[[3]] <- as.data.frame(lapply(x, as.logical))
#   ##24	RPPA.Cluster
#   x <- deleteBlanks(24)
#   names(x)
#   ListOfMat[[4]] <- as.data.frame(lapply(x, as.logical))
#   ##25	SigClust.Intrinsic.mRNA
#   x <- deleteBlanks(25)
#   names(x)
#   dim(x)
#   ListOfMat[[5]] <- as.data.frame(lapply(x, as.logical))
#   ##26	SigClust.Unsupervised.mRNA
#   x <- deleteBlanks(26)
#   dim(x)
#   ListOfMat[[6]] <- as.data.frame(lapply(x, as.logical))
#   ##28	Tumor
#   x <- deleteBlanks(28)
#   ListOfMat[[7]] <- as.data.frame(lapply(x, as.logical))
#   ##6	Converted.Stage
#   i <- 8
#   x <- deleteBlanks(6)
#   ListOfMat[[i]] <- as.data.frame(lapply(x, as.logical))
#   ##8	ER.Status
#   ##22 22	PR.Status
#   i <- i+1
#   x <- dichoSeveralCols(c(8,22))
#   ListOfMat[[i]] <- as.data.frame(lapply(x, as.logical))
#
#   #check all have the same dim
#   sapply(ListOfMat,dim)
#
#   #concatenate
#   conc2 <- ListOfMat[[1]]
#   for(i in seq_len(length(ListOfMat)-1))
#   {
#     conc2 <- data.frame(conc2,ListOfMat[[i+1]], check.names = TRUE, stringsAsFactors = FALSE)
#   }
#   #verify number of columns
#   sum(sapply(ListOfMat,dim)[2,])
#   dim(conc2)
#
#
#
#   ############
#   #5- add other columns, one by one
#
#   #---------
#   # ##3- form_completion_date
#   # yearCompletionF <- as.numeric(substr(clinical[,3],1,4))
#   #
#   # ##Recode categories to columns
#   # boolClinical <- data.frame(yearCompletionF)
#   # boolClinical <- lapply(boolClinical, factor)
#   # matYear <- model.matrix(~ . + 0, data=boolClinical,
#   #                          contrasts.arg = lapply(boolClinical, contrasts, contrasts=FALSE))
#   # matYear <- as.data.frame(matYear)
#
#
#   #---------
#   #using Sturges
#   ##7	Diagnosis.Age
#   aCOl <- clinical[,7]
#   range(aCOl)
#
#   aMat <- classesSturges(7)
#
#
#   # #concatenate step 5
#   # conc3 <- cbind(matYear,matAge,matDxYear)
#
#   conc3 <- aMat
#
#   ############
#
#   #6- concatenate all cols
#   boolClinical <- cbind(conc1, conc2, conc3)
#   dim(boolClinical)
#   row.names(boolClinical) <- row.names(clinical)
#
#   ##############
#   #PLOTS
#
#   #description boolean table
#   clinicalColunms <- names(boolClinical)
#   dim(boolClinical)
#
#   #plot distributions
#   #1 Number of clinical variables per patient
#   cli4patient <- apply(boolClinical,1,sum)
#   range(cli4patient)
#
#   plot(cli4patient,main = "Number of clinical variables per patient",
#        ylab = "number of true clinical variables",
#        xlab = "patient")
#
#   hist(cli4patient, main = "Histogram of clinical variables per patient (clinical variables count)",
#        ylab = "frequency of true clinical variables",
#        xlab = "patient")
#
#   #NOTE: Patients to remove (If remove patients, remove from mutations also)
#   which(cli4patient<=0)
#
#   #2 Number of patients with clinical variable x
#   pat4cli <- apply(boolClinical,2,sum)
#
#   #NOTE: Clinical to remove zero patients with the variable
#   length(which(pat4cli<=0))
#   #boolClinical <- boolClinical[,-which(pat4cli<=0)]
#
#   #NOTE: Clinical to remove all patients with the variable
#   #length(which(pat4cli==dim(boolClinical)[1]))
#
#   plot(pat4cli,main = "Number of patients per clinical variable",
#        ylab = "number of patients with true clinical variable x",
#        xlab = "clinical variable")
#
#   hist(pat4cli, main = "Histogram patient per clinical variable (patient count)",
#        ylab = "frequency of number of patients with true clinical variable x",
#        xlab = "clinical variable", breaks = 40)
#
#
# }
#
#
# #save.image("temp/10-BoolMatrices.RData")
# load("temp/10-BoolMatrices.RData")
#
