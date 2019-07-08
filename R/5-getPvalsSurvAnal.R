rm(list = ls())
source("R/5-survivalAnalysis.R")

# #BRCA
# load("data/patients.RData")
# load("dataExample/survival/survivalClinical.RData")
# row.names(survClinical) <- as.character(survClinical$patient.bcr_patient_barcode)
# survClinical <- survClinical[patients,]
# survivalClinical <- survClinical[,c( "patient.vital_status", "OS_MONTHS")]
# survivalClinical$patient.vital_status = survivalClinical$patient.vital_status != "alive"
# names(survivalClinical) <- c( "Overall.Survival.Status", "Overall.Survival..Months.")

#LUAD
load("data/patients.RData")
load("data/survivalVars.RData")
survivalClinical <-survivalVars
survivalClinical <- survivalClinical[patients,]
survivalClinical$Vital_Status = survivalClinical$Vital_Status != "living"
names(survivalClinical) <- c( "Overall.Survival..Months.", "Overall.Survival.Status" )

threshold_groups = 10 #10 or 3

rRange = c(3:90)

clustersK <- c(2:15) 

results <- data.frame(matrix(NA, nrow = length(clustersK), ncol = length(rRange) ))
row.names(results) <- clustersK
names(results) <- rRange

resultsOK <- results
numGroups <- results
numSamples <- results

# inputFolder = "dataExample/affiliations" 
inputFolder = "temp"

#K means
prefix = paste0(inputFolder, "/affiListKmeans_r") #r
outputN = "kmeans_"
# #HC means
# prefix = paste0(inputFolder, "/affiListHierClu_r") #r
# outputN = "hierarc_"


# # #K means without OL:
# prefix = paste0(inputFolder, "/affiliationList_k") #r
# outputN = "kmeans_without_OL"
# # #K means with OL:
# prefix = paste0(inputFolder, "/affiliationListWithOutL_k") #r
# outputN = "kmeans_with_OL"
# # # H clustering without OL:
# prefix = paste0(inputFolder, "/affiliationHC_k") #r
# outputN = "HC_without_OL"
# # # H clustering with OL:
# prefix = paste0(inputFolder, "/affiliationHCwithOL_k") #r
# outputN = "HC_with_OL"
# # #sprerical K means without OL:
# prefix = paste0(inputFolder, "/affiliationListSp_r") #r
# outputN = "skmeans_without_OL"
# # #K means with OL:
# prefix = paste0(inputFolder, "/affiliationListWithOutL_sk") #r
# outputN = "skmeans_with_OL"

outputN = paste0(outputN,'_t',threshold_groups,'_')

# Survival Analysis -----------
for(r in rRange)
{
  for(k in clustersK)
  {
    load(file=paste(prefix, r,"_k",k,".RData",sep=''))
    list[pval, ok, patiF, nFactors] = survivalAnalysis(r,k, threshold_groups=threshold_groups)
    results[as.character(k),as.character(r)] <- as.numeric(pval)
    resultsOK[as.character(k),as.character(r)] <- ok
    numGroups[as.character(k),as.character(r)] <- nFactors
    numSamples[as.character(k),as.character(r)] <- dim(patiF)[[1]]
    rm(list = ls()[-which(ls() %in% c("k","r","numSamples","survivalClinical" ,"rRange","clustersK", "prefix", "threshold_groups", "numGroups", "outputN", "survivalAnalysis", "results", "resultsOK")  )])
  }
}

# write.csv(results,file=paste("finalResults/results_", outputN, ".csv",sep='') )
# write.csv(resultsOK,file=paste("finalResults/resultsOK_", outputN, ".csv",sep='') )
# write.csv(numGroups,file=paste("finalResults/numGroups_", outputN, ".csv",sep='') )
# write.csv(numSamples,file=paste("finalResults/numSamples_", outputN, ".csv",sep='') )


finalResults <- data.frame(matrix(NA, nrow = length(clustersK), ncol = 3*length(rRange) ))
row.names(finalResults) <- clustersK

j=1
for(i in  seq_along(rRange) ){
  finalResults[,j] <- results[,i]
  finalResults[,j+2] <- numGroups[,i]
  finalResults[,j+1] <- numSamples[,i]
  j=j+3
}


#write.csv(finalResults,file=paste("finalResults/finalTable_", outputN, ".csv",sep='') )

View(finalResults)

minK <- apply(results,1,min)
which.min(minK)
min(minK)
minR <- apply(results,2,min)
which.min(minR)
min(minR)
