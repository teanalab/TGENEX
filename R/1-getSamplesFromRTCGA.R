# Construct Binary Matrix of Breast Cancer dataset
# using the RTCGA package to download the data and
# write a csv file to be read by python for training

library(RTCGA)
print(checkTCGA('Dates'))
#From my gist get entrez from hugo symbol
library(org.Hs.eg.db)
library(dplyr)


(cohorts <- infoTCGA() %>%
    rownames() %>%
    sub("-counts", "", x=.))


#options (move all the resulting folders before getting the next set of data)
includeSilentMutations = FALSE
allGenes = TRUE
entrez = FALSE #Turn this true to save entrez IDs instead of the gene symbols

#' Added parameter file name
#' get mutation data directly from RTCGA
#' this excludes silent mutations
#'
#' @return
#' @export
#'
#' @examples
#' load("data/OV/loadls.RData")
#' source('~/10-Research/Current/CLINGEN/9source/CLIGEN/publicPackageCLIGEN/R/1-createBoolMutation.R', echo=F)
#' mutationFromRTCGA()
mutationFromRTCGA <- function(disease = "BRCA", includeSilentMutations,
                              fileNameMutation = paste0("RDataFiles/binaMutation",disease,".RData"),
                              fileNamePatients = paste0("RDataFiles/patients",disease, ".RData") )
{
  library(RTCGA.mutations)
  library(dplyr)
  library(reshape2)

  if (!includeSilentMutations){
    if(disease == "BRCA") {
      mutationsTCGA(BRCA.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "OV") {
      mutationsTCGA(OV.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "LUAD") {
      mutationsTCGA(LUAD.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "ACC") {
      mutationsTCGA(ACC.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "BLCA") {
      mutationsTCGA(BLCA.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "CESC") {
      mutationsTCGA(CESC.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "CHOL") {
      mutationsTCGA(CHOL.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "COAD") {
      mutationsTCGA(COAD.mutations) %>%
        filter(Variant_Classification != "Silent") %>% # cancer tissue
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else {
      stop("disease not included")
    }
  } else {
    if(disease == "BRCA") {
      mutationsTCGA(BRCA.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "OV") {
      mutationsTCGA(OV.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "LUAD") {
      mutationsTCGA(LUAD.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "ACC") {
      mutationsTCGA(ACC.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "BLCA") {
      mutationsTCGA(BLCA.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "CESC") {
      mutationsTCGA(CESC.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "CHOL") {
      mutationsTCGA(CHOL.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else if(disease == "COAD") {
      mutationsTCGA(COAD.mutations) %>%
        filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% # cancer tissue
        mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12)) -> CANCER.mutations
    } else {
      stop("disease not included")
    }
  }

  patients <- unique(CANCER.mutations$bcr_patient_barcode)
  genes <- unique(CANCER.mutations$Hugo_Symbol)

  binaMutation <- data.frame(matrix(0,length(patients),length(genes)),
                             row.names = patients)
  names(binaMutation) <- genes

  for(i in seq_len(dim(CANCER.mutations)[1])){
    Apatient <- as.character(CANCER.mutations$bcr_patient_barcode[i])
    APgene <- as.character(CANCER.mutations$Hugo_Symbol[i])
    binaMutation[Apatient,APgene] <- 1
  }


  save(binaMutation, file=fileNameMutation)
  save(patients,  file=fileNamePatients)


  return (list(binaMutation,patients))
}


#' get gene expression data directly from RTCGA
geneExpressionFromRTCGA <- function( )
{
  library(RTCGA.mRNA)
  library(dplyr)
  library(reshape2)


  dir.create( "data2" )
  releaseDate <- "2016-01-28"
  sapply( cohorts, function(element){
    tryCatch({
      downloadTCGA( cancerTypes = element,
                    dataSet = "Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3",
                    destDir = "BRCA",
                    date = releaseDate )},
      error = function(cond){
        cat("Error: Maybe there weren't mutations data for ", element, " cancer.\n")
      }
    )
  })


  list.files( "data2") %>%
    file.path( "data2", .) %>%
    file.rename( to = substr(.,start=1,stop=50))


  list.files( "data2") %>%
    file.path( "data2", .) %>%
    sapply(function(x){
      file.path(x, list.files(x)) %>%
        grep(pattern = "MANIFEST.txt", x = ., value=TRUE) %>%
        file.remove()
    })


  list.files("data2") %>%
    file.path("data2", .) %>%
    sapply(function(y){
      file.path(y, list.files(y)) %>%
        assign( value = .,
                x = paste0(list.files(y) %>%
                             gsub(x = .,
                                  pattern = "\\..*",
                                  replacement = "") %>%
                             gsub(x=., pattern="-", replacement = "_"),
                           ".mRNA.path"),
                envir = .GlobalEnv)
    })

  ls() %>%
    grep("mRNA\\.path", x = ., value = TRUE) %>%
    sapply(function(element){
      tryCatch({
        readTCGA(get(element, envir = .GlobalEnv),
                 dataType = "mRNA") %>%
          assign(value = .,
                 x = sub("\\.path", "", x = element),
                 envir = .GlobalEnv )
      }, error = function(cond){
        cat(element)
      })
      invisible(NULL)
    }
    )


  cohorts <- ls()
  cohorts <- cohorts[grep("mRNA$", cohorts)]


  save.image("allmRNA.RData")
  save(list=cohorts, file="cohortsMRNA.RData")
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
survivalFromRTCGA <- function(disease)
{
  require(RTCGA.clinical)
  require(dplyr)
  if(disease == "LUAD") {
    survivalTCGA(LUAD.clinical) -> survivalVars
  } else if(disease == "BRCA") {
    survivalTCGA(BRCA.clinical) -> survivalVars
  } else if(disease == "OV") {
    survivalTCGA(OV.clinical) -> survivalVars
  } else if(disease == "ACC") {
    survivalTCGA(ACC.clinical) -> survivalVars
  } else if(disease == "BLCA") {
    survivalTCGA(BLCA.clinical) -> survivalVars
  } else if(disease == "CESC") {
    survivalTCGA(CESC.clinical) -> survivalVars
  } else if(disease == "CHOL") {
    survivalTCGA(CHOL.clinical) -> survivalVars
  } else if(disease == "COAD") {
    survivalTCGA(COAD.clinical) -> survivalVars
  } else {
    eval(parse(text= paste("survivalTCGA(",disease,".clinical) -> survivalVars", sep='')))
  }
  return(survivalVars)
}



# testing functions #############

#
# disease = cohorts[1]
# dir.create(disease)
# #X = mutationFromRTCGA(disease, includeSilentMutations)
# #binaMutation <- X[[1]]
# #patients <- X[[2]]
# #rm(X)
#
# if(!allGenes){
#   print( "***********\nTO FILTER RUN THE SCRIPTS IN THE scripts2filterData FOLDER.\n***************\n" )
# }
#
#
# #get gene expression
#
#
#
# #get clinical data
#
#
#
#
# #get survival
# survivalVars <- survivalFromRTCGA(disease)
# row.names(survivalVars) <- survivalVars$bcr_patient_barcode
# #patients <- intersect(patients, row.names(survivalVars))
# #survivalData <- survivalVars[patients,-2]
# survivalData <- survivalVars
#
# cat(dim(survivalVars))
#
# #Plot for paper - complete example in SurvivalAnalysis.R
# plotName<-paste0(disease,"/OS_Kaplan-Meier.pdf")
# pdf( file = plotName,  onefile = TRUE, width = 9, height = 7)
# ## Kaplan-Meier Survival Curves
# kmTCGA(survivalData)
# dev.off()
#
# #write mutation to csv file for python script
# binaMutation <- binaMutation[patients,]
# tBina = t(binaMutation)
# if(!entrez){
#   #remove the genes that have no entrez eventhough we will use only gene symbol?
#   #no, because, these information can be useful and if we discover something that the other method
#   #cannot because it is not in the protein networks, then our method is in advantage compare with the others
#   write.table(tBina, file = paste0(disease,"/mutation.csv"),
#               sep =';', dec = '.', row.names = T, col.names = T )
# } else {
#   stop("I have not tested these functionality, but I have the script fromSym2Entrez on the scripts2smoothData project")
#
# }
#
# write(c("Mutations table\nNumber of patients: ", dim(binaMutation)[1],
#         "Number of genes: ", dim(binaMutation)[2]),
#         file=paste0(disease, "/dimensions.txt"), append=TRUE )
#
#
# #write survival data for python script
# write.table(survivalData, file = paste0(disease, "/survival.csv"), sep =';',
#             dec = '.', row.names = T, col.names = T )
#
# write(c("\nClinical table\nNumber of patients: ", dim(survivalData)[1]),
#       file=paste0(disease, "/dimensions.txt"), append=TRUE )



# # ---- main script ------
# #1- Example how to get Breast Cancer mutation data using RTCGA
# diseases = c("BRCA","OV","BLCA","CESC", "CHOL", "COAD", "LUAD", "ACC")
#
# for (i in seq_along(diseases))
# {
#   disease = diseases[i]
#   dir.create(disease)
#   X = mutationFromRTCGA(disease, includeSilentMutations)
#   binaMutation <- X[[1]]
#   patients <- X[[2]]
#   rm(X)
#
#   if(!allGenes){
#     print("***********\nTO FILTER RUN THE SCRIPTS IN THE scripts2filterData FOLDER.\n***************\n")
#   }
#
#   #get survival
#   survivalVars <- survivalFromRTCGA(disease)
#
#   #Plot for paper - complete example in SurvivalAnalysis.R
#   plotName<-paste0(disease,"/OS_Kaplan-Meier.pdf")
#   pdf( file = plotName,  onefile = TRUE, width = 9, height = 7)
#   ## Kaplan-Meier Survival Curves
#   kmTCGA(survivalData)
#   dev.off()
#
#   #write mutation to csv file for python script
#   binaMutation <- binaMutation[patients,]
#   tBina = t(binaMutation)
#   if(!entrez){
#     #remove the genes that have no entrez eventhough we will use only gene symbol?
#     #no, because, these information can be useful and if we discover something that the other method
#     #cannot because it is not in the protein networks, then our method is in advantage compare with the others
#     write.table(tBina, file = paste0(disease,"/mutation.csv"),
#                 sep =';', dec = '.', row.names = T, col.names = T )
#   } else {
#     stop("I have not tested these functionality, but I have the script fromSym2Entrez on the scripts2smoothData project")
#
#   }
#
#   write(c("Mutations table\nNumber of patients: ", dim(binaMutation)[1],
#           "Number of genes: ", dim(binaMutation)[2]),
#         file=paste0(disease, "/dimensions.txt"), append=TRUE )
#
#
#   #write survival data for python script
#   write.table(survivalData, file = paste0(disease, "/survival.csv"), sep =';',
#               dec = '.', row.names = T, col.names = T )
#
#   write(c("\nClinical table\nNumber of patients: ", dim(survivalData)[1]),
#         file=paste0(disease, "/dimensions.txt"), append=TRUE )
# }





# # ---- main script get survival data------
# #1- Example how to get Breast Cancer mutation data using RTCGA
# (diseases = cohorts[11:38])
#
# for (i in seq_along(diseases))
# {
#   disease = diseases[i]
#   dir.create(disease)
#
#   #get survival
#   survivalVars <- survivalFromRTCGA(disease)
#   row.names(survivalVars) <- survivalVars$bcr_patient_barcode
#   survivalData <- survivalVars
#
#   #Plot for paper - complete example in SurvivalAnalysis.R
#   plotName<-paste0(disease,"/OS_Kaplan-Meier.pdf")
#   pdf( file = plotName,  onefile = TRUE, width = 9, height = 7)
#   ## Kaplan-Meier Survival Curves
#   kmTCGA(survivalData)
#   dev.off()
#
#
#   #write survival data object
#   save(survivalData, file = paste0(disease, "/survival.Rd") )
# }




#' Title
#'
#' @param disease
#'
#' @return
#' @export
#'
#' @examples
#' load(file="cohortsMRNA.RData")
#' disease = "BRCA"
getOnlyTumorFromRNAandFilter <- function(disease, folder_p){
  eval(parse(text= paste0("TCGA_DB <- ",disease,".mRNA")))
  patients <- unique(TCGA_DB$bcr_patient_barcode)
  sampletypes <- substr(patients, 14, 15)
  print(table(sampletypes))
  #https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
  # 01 primary tummor
  # 06 Metastatic
  # 11	Solid Tissue Normal
  #remove non primary tummor samples
  toDel = which(sampletypes != '01')
  if (length(toDel) > 0 ){
    TCGA_DB <- TCGA_DB[-toDel, ]
    patients <- TCGA_DB$bcr_patient_barcode
  }
  patients <- substr(patients, 1, 12)
  if( length(unique(patients)) != dim(TCGA_DB)[1])
    stop("ERROR: length(unique(patients)) != dim(TCGA_DB)[1]")
  TCGA_DB$bcr_patient_barcode <- patients
  #remove NA
  if(anyNA(TCGA_DB)){
    nCols = dim(TCGA_DB)[2]
    cat(sum(is.na(TCGA_DB)), " NA values found for ", disease, "\n")
    #delete genes because patients are more valuable
    delCols <- apply(TCGA_DB, 2, anyNA )
    cat(sum(delCols), " genes to be deleted\n")
    remNA = which(delCols)
    TCGA_DB <- TCGA_DB[,-remNA]
    cat( (nCols - dim(TCGA_DB)[2] ) , "genes were deleted", "\n")
  }
  if ( anyNA(TCGA_DB) ){
    stop("ERROR: Table still have NA values for ", disease)
  }
  row.names(TCGA_DB) <- TCGA_DB$bcr_patient_barcode
  TCGA_DB = TCGA_DB[, -1]
  eval(parse(text= paste0("TCGA_DB -> ",disease,".mRNA")))
  save(list=paste0(disease,".mRNA"), file = paste0(folder_p,disease,".RData"))
}
