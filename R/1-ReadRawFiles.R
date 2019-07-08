#what about this?
#rawClinicalFile="/Users/diamac/10-Research/Current/CLINGEN/9source/CLIGEN/publicPackage/CLIGEN/rawData/data_clinical.txt"


#' Read raw files from public datasets. reads the raw files from the 'rawData' folder and store them in 'dataFile'
#'
#' @param dataFile the RData file that will store the R objects with the mutation data
#' @param rawMutationFile the path to the raw data mutation file downloaded from CBioportal
#' @param rawClinicalFile the path to the raw clinical data file downloaded from firebrowse
#'
#' @examples
#' readBRCAFiles (dataFile = "temp/readBRCAFiles.RData",
#' rawMutationFile="rawData/data_mutations_extended.txt",
#' rawClinicalFile="rawData/BRCA.clin.merged.txt")
readBRCAFiles <- function(dataFile = "temp/old_moved_10-10/1-2v3data.RData",
                      rawMutationFile="/Users/diamac/DDfileSystem/10-Research/Current/CLINGEN/12Datasets/BRCA/brca4mCBioportal/data_mutations_extended.txt",
                      rawClinicalFile="/Users/diamac/DDfileSystem/10-Research/Current/CLINGEN/12Datasets/BRCA/gdacLevel_1/BRCA.clin.merged.txt"){
  #file.remove(dataFile)
  if( file.exists(dataFile) ){
    load(dataFile, env = globalenv())
  } else {
    clinicalLong<-read.table(file=rawClinicalFile, sep = "\t", quote = '"',
                             header = FALSE, stringsAsFactors = FALSE,  fill = TRUE)
    mutation <- read.table(file=rawMutationFile, sep = "\t", quote = "",
                          header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)

    #uuid = NO A GOOD IDENTIFIER
    #Tumor_Sample_Barcode  = GOOD IDENTIFIER, BUT WE NEED TO CUT
    #TCGA-06-5411-01A-01D-1696-08
    # mutationBarcodes <- substring(tolower(mutation$Tumor_Sample_Barcode), 1, 12)
    # mutation$mutationBarcodes <- mutationBarcodes


    # DATASETS DESCRIPTION
    # cat("Dimensions of Mutation file: ")
    # cat(dim(mutation))
    # mutationRange <- t(sapply(mutation, range))
    # cat("\nRange of Mutation file: ")
    # cat(mutationRange)
    # write.csv(mutationRange,file="descriptionData/brca/4cBiomutationRange.csv")

    # mutation$Tumor_Sample_Barcode[1:5] #"Tumor_Sample_Barcode"
    # #unique(mutation$Tumor_Sample_UUID) #Tumor_Sample_UUID
    # substring(mutation$Tumor_Sample_Barcode[1:5],1,12)

    mutation$mutationBarcodes <- tolower(substring(mutation$Tumor_Sample_Barcode,1,12)) #get just patient ID
    patientsIDinMut <- as.character(mutation$mutationBarcodes) #get just patient ID
    cat("\nNumber of patients in mutation file: ")
    cat( length(unique(patientsIDinMut)) )

    # reading clinical variables -----
    #If not transposed
    row.names(clinicalLong) <- clinicalLong[,1]
    VariablesLong <- row.names(clinicalLong)
    #colsWithID <- which(regexpr("barcode",VariablesLong, ignore.case = TRUE) > 0)
    #VariablesLong[colsWithID[1:10]]
    # [16] patient.bcr_patient_barcode
    #clinicalLong[16,2:5]

    patientsIDinClin <- clinicalLong["patient.bcr_patient_barcode",] #get just patient ID
    cat("\nNumber of patients in clinical file: ")
    cat( length(unique(patientsIDinClin)) )


    #number of patients
    patients <- intersect( as.character(unlist(patientsIDinClin)), unlist(patientsIDinMut) )
    cat("\nNumber of patients in both files: ")
    cat(length(patients))

    n <- which(patientsIDinClin %in% patients)
    clinicalLong <- clinicalLong[,n]
    clinicalLong <- as.data.frame(t(clinicalLong))

    # newpatients <- clinicalLong[,16]
    # length(intersect(patientsIDinClin, newpatients))
    # write.csv(clinicalLong, "temp/clinicalLong.csv")

    save(list=c("patients", "mutation", "clinicalLong"), file=dataFile)
    load(dataFile, env = globalenv())
  }
}



#' Filter only useful clinical variables. This requires that readBRCAFiles() need
#' to be run first and the file "rawData/clinical_useful.csv" needs to be present.
#' Write a file named 'clinical_useful.csv' with the list of clinical variables that
#' you want to include in your analysis. Locate the file in the folder 'rawData'.
#' This csv file should have the following column names:
#' "column", "USEFUL", "type".
#' "column" refers to the name of the clinical variable in the raw file,
#' "USEFUL" has a 'Y' or 'N' value indicating if the clinical variable is useful
#' for the analysis or not,
#' and "type" can have the values 'nominal' or 'continuos' indicating the type of
#' useful variables.
#' This function generates latex tables with the statistics of the clinical variables
#' and store them in the 'temp' folder. The file names are:  "temp/ContinousVars.tex"
#' for statistics on the continuos variables
#' and "temp/NominalVars.tex" for statistics on the nominal variables.
#'
#' @param mergeClinical This indicates if you want to merge the clinical file from
#' CBioportal or not. If this value is true the file data_clinical.txt should exist in the rawData folder.
#' In adition, the file "merge_useful.csv" must exist and it should have the following columns
#' "columns", "repeated", "type", "for.tensor". See the provided file as an example.
#'
#' @requires reporttools
#'
#' @examples
#' #We provide an example with the clinical variables that are relevant to Breast Cancer.
#'
filterClinicalVariables <- function(mergeClinical=T, tensorClinical){
  require("reporttools")
  if(mergeClinical){
    clini_use_cols<-read.table(file="rawData/clinical_useful.csv", sep = ",",
                               header = TRUE, stringsAsFactors = FALSE, na.strings = "nasrting", fill = TRUE)

    useClini <- clinicalLong[, clini_use_cols[which(clini_use_cols[,2]=="Y") ,"column"]]
    row.names(useClini) <- useClini$patient.bcr_patient_barcode

    nominalCols <- clini_use_cols$column[which(clini_use_cols$type == "nominal")]
    cap1 <- "Patient characteristics: nominal variables. firebrowse"
    write(tableNominal(vars = useClini[,nominalCols], cap = cap1, vertical = TRUE, lab = "tab: nominal1",
                       longtable = TRUE), file = "temp/NominalVars.tex")

    continousCols <- clini_use_cols$column[which(clini_use_cols$type == "continuos")]
    numClini <- useClini[,continousCols]
    indx <- sapply(numClini, is.factor)
    numClini[indx] <- lapply(numClini[indx], function(x) as.numeric(as.character(x)))
    cap1 <- "Patient characteristics: continuous variables.firebrowse"
    write(tableContinuous(vars = numClini, cap = cap1, lab = "tab: continuous2",
                          longtable = TRUE), file = "temp/ContinousVars.tex")


    # read clinical to merge
    clinical1<-read.table(file="rawData/data_clinical.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    clinical1$PATIENT_ID <- tolower(substring(clinical1$PATIENT_ID,1,12)) #get just patient ID
    row.names(clinical1) <- clinical1$PATIENT_ID

    #length(unique (intersect (useClini$patient.bcr_patient_barcode, clinical1$PATIENT_ID)))
    mergedClini <- merge(x = useClini, y = clinical1, by.x = "patient.bcr_patient_barcode", by.y = "PATIENT_ID")

    # verify repeated columns
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
    # which(mergedClini$patient.breast_carcinoma_progesterone_receptor_status != tolower(mergedClini$PR.Status) )
    #
    # tempStatus <- mergedClini$OS_STATUS
    # tempStatus[which(tempStatus == "LIVING")] = "alive"
    # tempStatus[which(tempStatus == "DECEASED")] = "dead"
    # which(mergedClini$patient.vital_status != tempStatus )

    # remove patients ------
    mergedClini <- mergedClini[-c(296,471),]

    # remove repeated and non useful columns
    merge_use_cols<-read.table(file="rawData/merge_useful.csv", sep = ",",
                               header = TRUE, stringsAsFactors = FALSE, na.strings = "nasrting", fill = TRUE)

    name_useful_mer_col <- merge_use_cols[which(merge_use_cols$repeated == '' ),1]
    mergedClini <- mergedClini[,name_useful_mer_col]

    # describe remaning variables
    nominalCols <- merge_use_cols$columns[which(merge_use_cols$type == "nominal")]
    cap1 <- "Patient characteristics: nominal variables. Merged firebrowse and cBio"
    write(tableNominal(vars = mergedClini[,nominalCols], cap = cap1, vertical = TRUE, lab = "tab: nominal1",
                       longtable = TRUE), file = "temp/NominalVars.tex")

    continousCols <- merge_use_cols$columns[which(merge_use_cols$type == "continuos")]
    numClini <- mergedClini[,continousCols]
    indx <- sapply(numClini, is.factor)
    numClini[indx] <- lapply(numClini[indx], function(x) as.numeric(as.character(x)))
    cap1 <- "Patient characteristics: continuous Merged firebrowse and cBio"
    write(tableContinuous(vars = numClini, cap = cap1, vertical = TRUE, lab = "tab: continuous2",
                          longtable = TRUE), file = "temp/ContinousVars.tex")
    cat("\n\nFor some reason the last file is not writen, so MOVE OUTPUT MANUALLY TO FILE!!! temp/ContinousVars.tex")

    #maybe try to generate other statistics with https://cran.r-project.org/web/packages/Gmisc/vignettes/Descriptives.html


    # construct clinical tables -----
    merge_use_cols[15,c(2:4)]  <- c("no_useful","","")
    merge_use_cols[16,c(2:4)]  <- c("no_useful","","")
    tensorCols <- merge_use_cols$columns[which(merge_use_cols$for.tensor == "tensor")]
    tensorClinical <- mergedClini[,c("patient.bcr_patient_barcode", tensorCols)]
    save(tensorClinical, file="data/tensorClinical.RData")

    survCols <- merge_use_cols$columns[which(merge_use_cols$for.tensor == "survival")]
    survClinical <- mergedClini[,c("patient.bcr_patient_barcode", survCols)]
    save(survClinical, file="data/survClinical.RData")

    treatCols <- merge_use_cols$columns[which(merge_use_cols$for.tensor == "treatment")]
    treatClinical <- mergedClini[,c("patient.bcr_patient_barcode", treatCols)]
    save(treatClinical, file="data/treatClinical.RData")

    patients <- as.character(mergedClini$patient.bcr_patient_barcode)
    save(patients,file="data/patients.RData")

    write.csv(as.data.frame(t(mergedClini)), "temp/mergedClini.csv")
  }

}


#' get Survival data directly from RTCGA
#' more at https://rtcga.github.io/RTCGA/reference/survivalTCGA.html
#'
#' @return
#' @export
#'
#' @examples
#' load("data/OV/loadls.RData")
#' source('~/10-Research/Current/CLINGEN/9source/CLIGEN/publicPackageCLIGEN/R/1-createBoolMutation.R', echo=F)
#' mutationFromRTCGA()
survivalFromRTCGA <- function(disease ="LUAD")
{
  require(RTCGA.clinical)
  require(dplyr)
  if(disease == "LUAD") {
      survivalTCGA(LUAD.clinical) -> survivalVars
  } else {
    stop("disease not included")
  }
  return(survivalVars)
}
