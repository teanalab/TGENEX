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
      clinical <<- clinical[-toDelete,]
      patients <<- patients[-toDelete]
      boolMutation <<- boolMutation[-toDelete,]
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
  #head(matDich2)
  #View(matDich2)
  #delete redundant info, leave only one column
  #reduCols <- seq(1,(2*length(which(NumVals==2))),by=2)
  #select manually for true being higher risk
  #names(matDich2)
  #reduCols <- c(2,4,5,7,9,11)
  reduCols <- c(2,4,5,7,9)
  matDich2 <- matDich2[,-reduCols]
  
  
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  NumVals
  
  
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
  #9 NA, got rid of patients - 8 patients
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
  
  
  NumVals <- apply(clinical, 2 ,function(x) { length(unique(x))})
  colsN <- which(NumVals==6)
  ncols <- names(colsN)
  #matDich7 <- dichoSeveralCols(ncols[1])
  
  
  #Do like it has 6 levels
  matDich6 <- dichoSeveralCols(ncols[1])
  
  
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
  matDich6 <- matDich6[patients,]
  matDich8 <- matDich8[patients,]
  
  conc1 <- cbind(matDich2, matDich3, matDich6, matDich8)
  
  #5- add other columns, one by one --------
  
  #using Sturges
  ##2	Diagnosis.Age
  C <- "patient.age_at_initial_pathologic_diagnosis"
  aCOl <- as.numeric(as.character(clinical[patients,C]))
  range(aCOl) #29 89
  
  # here -------------
  load("temp/imhere.RData")
  aMat <- classesSturges(C)
  
  #   style: jenks
  #   one of 8,996,462,475 possible partitions of this variable into 10 classes
  #   [3,10] (10,16] (16,22] (22,27] (27,32] (32,37] (37,42] (42,48] (48,54] (54,60]
  #     25      46      52      63      56      74      59      32      32      17
  #   [29,37] (37,43] (43,49] (49,54] (54,59] (59,64] (64,70] (70,76] (76,82] (82,89]
  #     25      46      52      63      56      74      65      33      29      13
  
  # aMatNew  <- classesSturges(C)
  
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
  
  pa4cli <- apply(boolClinical,2,sum)
  length(pa4cli)
  range(pa4cli)
  
  boolClinical <- boolClinical[,-which(pa4cli<=5)]
  
  
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




getUsefulColumns <- function(disease = "BRCA"){
  if( disease == "BRCA"){
    clinicalData <<- BRCA.clinical[which(BRCA.clinical[,"patient.bcr_patient_barcode"] %in% tolower(as.character(unlist(patients)) ) ) , ]
  } else if (disease == "OV"){
    clinicalData <<- OV.clinical[which(OV.clinical[,"patient.bcr_patient_barcode"] %in% tolower(as.character(unlist(patients)) ) ) , ]
  } else {
    stop("disease not included")
  }
  if(disease == "BRCA"){
    #patient.gender 
    toKeepC <- which( regexpr( "gender", names(clinicalData) ) != -1 )
    usefulCols <- c(names(clinicalData)[toKeepC])
    
    #patient.ethnicity 
    toKeepC <- which( regexpr( "ethnicity", names(clinicalData) ) != -1 )
    # names(clinicalData)[toKeepC]
    # View(clinicalData[,toKeepC])
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC[1]] )
    
    #age_at_initial_pathologic_diagnosis  
    toKeepC <- which( regexpr( "age_at_initial_pathologic_diagnosis", names(clinicalData) ) != -1 )
    #names(clinicalData)[toKeepC]
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC])
    
    # Stage
    # toKeepC <- which( regexpr( "stage", names(clinicalData) ) != -1 )
    # names(clinicalData)[toKeepC]
    # View(clinicalData[,toKeepC])
    toKeepC <- which( regexpr( "pathologic_stage", names(clinicalData) ) != -1 )
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC])
    
    # Tumor
    # Node
    # Metastasis    
    toKeepC <- which( regexpr( "pathologic_categories.pathologic_", names(clinicalData) ) != -1 )
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC])
                    
    # Histology
    toKeepC <- which( regexpr( "histology", names(clinicalData) ) != -1 )
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC])
    
    # HER2.Final.Status
    toKeepC <- which( regexpr( "patient.lab_proc_her2", names(clinicalData) ) != -1 )
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC])
    
    #progesterone_receptor_status
    toKeepC <- which( regexpr( "progesterone_receptor_status", names(clinicalData) ) != -1 )
    # names(clinicalData)[toKeepC]
    # View(clinicalData[,toKeepC])
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC[1]] )
    
    # estrogen_receptor_status
    toKeepC <- which( regexpr( "estrogen_receptor_status", names(clinicalData) ) != -1 )
    # names(clinicalData)[toKeepC]
    # View(clinicalData[,toKeepC])
    usefulCols <- c(usefulCols, names(clinicalData)[toKeepC[1]] )
    return(usefulCols)
  } else if(disease == "OV"){
    usefulCols <- c("patient.age_at_initial_pathologic_diagnosis", "patient.bcr_patient_barcode")
    toKeepC <- which( regexpr( "race", names(clinicalData) ) != -1 )
    names(clinicalData)[toKeepC]
    usefulCols <- c(usefulCols, "patient.race")
    toKeepC <- which( regexpr( "histo", names(clinicalData) ) != -1 )
    names(clinicalData)[toKeepC]
    usefulCols <- c(usefulCols,"patient.neoplasm_histologic_grade" ,
                    "patient.biospecimen_cqcf.histological_type",
                    "patient.biospecimen_cqcf.history_of_neoadjuvant_treatment" )
    toKeepC <- which( regexpr( "stage", names(clinicalData) ) != -1 )
    names(clinicalData)[toKeepC]
    usefulCols <- c(usefulCols,"patient.stage_event.clinical_stage")
  } else {
    #remove all NAs
    delCols <- apply(clinicalData,2, function(x){ all(is.na(x)) } )
    delC <- which(delCols)
    clinicalData <- clinicalData[,-delC]

    #Num different values
    numDiffCols <- apply(clinicalData, 2, function(x){ length(unique(x)) } )
    delC <- which(numDiffCols==1)      #only one value
    clinicalData <- clinicalData[,-delC]

    numDiffCols <- apply(clinicalData, 2, function(x){ length(unique(x)) } )
    delC <- which(numDiffCols==dim(clinicalData)[1])      #all values different
    clinicalData <- clinicalData[,-delC]

    # delC <- which( regexpr( "file", names(clinicalData) ) != -1 )
    # clinicalData <- clinicalData[,-delC]

    #COMPLETE TO KEEP
    keepCols <- apply(clinicalData,2, function(x){ all(!is.na(x)) } )
    keepC <- which(keepCols)
    return(keepC)
  }
  
}


#' Get clinical data From RTCGA
#'
#' @return
#' @export
#'
#' @examples
#' load("data/loadls.RData")
#' clinicalFromRTCGA()
clinicalFromRTCGA <- function(disease ="BRCA"){
  loadls("RTCGA RTCGA.clinical",F)
  load(file="data/patients.RData")
  #data(OV.clinical)
  #names(OV.clinical)
    
  usefulCols <- getUsefulColumns("BRCA")  # get useful cols
  
  clinicalData <- clinicalData[ , usefulCols]
  # keepCols <- apply(clinicalData,2, anyNA )
  # keepCols
  
  
  #TODO function clean data for each cancer
  table(clinicalData$patient.race,useNA =  "always" )
  #NA = white
  clinicalData$patient.race[which(is.na(clinicalData$patient.race))] = "white"
  
  table(clinicalData$patient.neoplasm_histologic_grade, useNA =  "always" )
  #NA = gx which means that cannot be assessed
  clinicalData$patient.neoplasm_histologic_grade[which(is.na(clinicalData$patient.neoplasm_histologic_grade))]="gx"
  
  table(clinicalData$patient.biospecimen_cqcf.histological_type , useNA =  "always" )
  #delete column, it has only one value
  which(names(clinicalData) == "patient.biospecimen_cqcf.histological_type")
  clinicalData <- clinicalData[,-5]
  
  table(clinicalData$patient.biospecimen_cqcf.history_of_neoadjuvant_treatment , useNA =  "always" )
  #delete column, it has only one value
  which(names(clinicalData) == "patient.biospecimen_cqcf.history_of_neoadjuvant_treatment")
  clinicalData <- clinicalData[,-5]
  
  
  table(clinicalData$patient.stage_event.clinical_stage , useNA =  "always" )
  #NA= stage iiic
  clinicalData$patient.stage_event.clinical_stage[which(is.na(clinicalData$patient.stage_event.clinical_stage))]="stage iiic"
  
  row.names(clinicalData) <- clinicalData$patient.bcr_patient_barcode
  clinicalData <- clinicalData[,-2]
  
  patients <- as.character(row.names(clinicalData))
  
  save.image("temp/clinicalFromRTCGA.RData")
  
}
