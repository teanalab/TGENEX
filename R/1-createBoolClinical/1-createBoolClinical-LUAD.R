source("R/1-createBoolClinical/1-dicotomFunctions.R")

###### run clinical line by line -----
createBoolClinical <-function(clinical){
  ## start analysis - create a pdf in temp folder with the variables distribution
  NumVals <- analyzeClinicVariables()
  NumVals
  
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals
  
  #  table(NumVals)
  
  #1- get columns with only one value
  cols1 <- which(NumVals==1)
  matDich1 <- clinical[,cols1]
  
  #2- get dichotomous columns
  cols2 <- which(NumVals==2)
  matDich2 <- clinical[,cols2]
  
  #3-level columns-----------------
  cols3 <- which(NumVals==3)  
  #RESIDUAL_TUMOR
  
  patients <<- row.names(clinical)
  matDich3 <- dichoSeveralCols(cols=cols3, logicVal = FALSE)
  diffvals <- names(table(clinical[,cols3]))
  names(matDich3) <-  paste0(names(cols3),"_",diffvals)
  
  #4 level columns-----------------
  cols4 <- which(NumVals==4)
  #AJCC_NODES_PATHOLOGIC_PN AJCC_TUMOR_PATHOLOGIC_PT 
  
  matDich4 <- dichoSeveralCols(cols=cols4, logicVal = FALSE)
  names(matDich4)
  # diffvals1 <- names(table(clinical[,cols4[1]]))
  # diffvals2 <- names(table(clinical[,cols4[2]]))
  # names(matDich4) <-  paste0(names(cols4),"_",diffvals)
 
  #6 level columns-----------------
  cols6 <- which(NumVals==6)
  #AJCC_PATHOLOGIC_TUMOR_STAGE 
  
  matDich6 <- dichoSeveralCols(cols=cols6, logicVal = FALSE)
  names(matDich6)
  diffvals <- names(table(clinical[,cols6]))
  names(matDich6) <-  paste0(names(cols6),"_",diffvals)
  
  ### CONCATENATE ALL nominal
  matDich1 <- matDich1[patients,]
  matDich2 <- matDich2[patients,]
  matDich3 <- matDich3[patients,]
  matDich4 <- matDich4[patients,]
  matDich6 <- matDich6[patients,]

  conc1 <- cbind(matDich1, matDich2, matDich3, matDich4, matDich6)
  names(conc1) <- c(names(cols1), names(conc1)[2:51])
  
  #5- add other columns, one by one --------
  #using Sturges and jenks
  
  #9 level columns-----------------
  (cols9 <- which(NumVals==9))
  # SHORTEST_DIMENSION
  aCOl <- as.numeric(as.character(clinical[patients,cols9]))
  range(aCOl) 
  matDich9 <- classesSturges(cols9, 5) #fixed number of classes
  row.names(matDich9) <- patients
  
  #19 level columns-----------------
  (cols19 <- which(NumVals==19))
  # SPECIMEN_SECOND_LONGEST_DIMENSION
  aCOl <- as.numeric(as.character(clinical[patients,cols19]))
  range(aCOl) 
  matDich19 <- classesSturges(cols19) #classes defined using Sturges
  row.names(matDich19) <- patients
  
  #29 level columns-----------------
  (cols29 <- which(NumVals==29))
  # LONGEST_DIMENSION
  aCOl <- as.numeric(as.character(clinical[patients,cols29]))
  range(aCOl) 
  matDich29 <- classesSturges(cols29) #classes defined using Sturges
  row.names(matDich29) <- patients
  
  #47 level columns-----------------
  (cols47 <- which(NumVals==47))
  # SMOKING_YEARS
  aCOl <- as.numeric(as.character(clinical[patients,cols47]))
  range(aCOl) 
  matDich47 <- classesSturges(cols47) #classes defined using Sturges
  row.names(matDich47) <- patients
  
  #50 level columns-----------------
  (cols50 <- which(NumVals==50))
  # AGE
  aCOl <- as.numeric(as.character(clinical[patients,cols50]))
  range(aCOl) 
  matDich50 <- classesSturges(cols50) #classes defined using Sturges
  row.names(matDich50) <- patients
  
  
  #89 level columns-----------------
  (cols89 <- which(NumVals==89))
  # SMOKING_PACK_YEARS
  aCOl <- as.numeric(as.character(clinical[patients,cols89]))
  range(aCOl) 
  matDich89 <- classesSturges(cols89) #classes defined using Sturges
  row.names(matDich89) <- patients
  
  table(NumVals)
  ### CONCATENATE ALL CONTINUOUS
  conc2 <- cbind(matDich9, matDich19, matDich29, matDich47, matDich50, matDich89)
  
  #6- concatenate all cols
  boolClinical <- cbind(conc1, conc2)
  anyNA(boolClinical)
  dim(boolClinical)
  row.names(boolClinical) <- row.names(clinical)
  
  # plots ##########
  
  #description boolean table
  (clinicalColunms <- names(boolClinical))
  
  #plot distributions
  #1 Number of clinical variables per patient
  cli4patient <- apply(boolClinical,1,sum)
  range(cli4patient)
  
  pa4cli <- apply(boolClinical,2,sum)
  length(pa4cli)
  range(pa4cli)
  
  
  pdf(file = "temp/10-Fig6.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
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
  
  pdf(file = "temp/10-Fig7.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  plot(pat4cli,main = "Number of patients per clinical variable",
       ylab = "number of patients with true clinical variable x",
       xlab = "clinical variable")
  dev.off()
  
  pdf(file = "temp/10-Fig8.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  hist(pat4cli, main = "Histogram patient per clinical variable (patient count)",
       ylab = "frequency of number of patients with true clinical variable x",
       xlab = "clinical variable", breaks = 40)
  dev.off()
  
  # save ##########
  save(boolClinical,file= "data/boolClinical.RData")
}

