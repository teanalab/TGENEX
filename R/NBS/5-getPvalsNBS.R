rm(list = ls())
source("R/1-ReadRawFiles.R")
source("R/5-survivalAnalysis.R")

#LUAD
#patients from NBS
patientsAll <- read.table("/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/3_LUAD/output_standardNMF_sampleID.csv",
                       header = FALSE)
patientsAll <-as.character(unlist(patientsAll))
#survival from RTCGA
survivalVars <- survivalFromRTCGA("LUAD")
survivalClinical <-survivalVars
row.names(survivalClinical) <- survivalClinical[,2]
survivalClinical <- survivalClinical[patientsAll,]
survivalClinical$patient.vital_status  = survivalClinical$patient.vital_status == 1
survivalClinical <- survivalClinical[,-2]
names(survivalClinical) <- c( "Overall.Survival..Months.", "Overall.Survival.Status" )

#remove NA
if(anyNA(survivalClinical)) {
  remNA = which(is.na(survivalClinical$Overall.Survival..Months.))
  survivalClinical <- survivalClinical[-remNA,]
}

patients <- intersect(row.names(survivalClinical),patientsAll) #only 310 in survival file


threshold_groups = 3 #10 or 3


clustersK <- c(3:80) #80

results <- data.frame(matrix(NA, nrow = length(clustersK), ncol = 1 ))
row.names(results) <- clustersK

resultsOK <- results
numGroups <- results
numSamples <- results

# inputFolder = "dataExample/affiliations" 
outputN = "NBS_kmeans_"
outputN = paste0(outputN,'_t',threshold_groups,'_')

# Survival Analysis -----------
for(k in clustersK)
{
  clustersNBS <-read.table(file=paste0("/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/3_LUAD/output_clustersNBS_k",k,".csv"),
                           sep = "\t", quote = "",
                           header = FALSE, stringsAsFactors = FALSE)
  row.names(clustersNBS) <- patientsAll
  

  clustersNBS<-clustersNBS[patients,]
  
  affiliationList = clustersNBS
  affiliationList = unlist(affiliationList)
  names(affiliationList) = patients
  save(affiliationList, file = paste0("temp/affiliationNBSk",k,".RData") )
  
  
  
  list[pval, ok, patiF, nFactors] = survivalAnalysis('-',k, threshold_groups=threshold_groups)
  
  results[k-2] <- as.numeric(pval)
  resultsOK[k-2] <- ok
  numGroups[k-2] <- nFactors
  numSamples[k-2]<- dim(patiF)[[1]]
  rm(list = ls()[-which(ls() %in% c("k", "numSamples","survivalClinical" ,"rRange","clustersK", 
                                    "prefix", "threshold_groups", "numGroups", "outputN", 
                                    "survivalAnalysis", "results", "resultsOK",
                                    "patientsAll")  )])
}


# write.csv(results,file=paste("finalResults/results_", outputN, ".csv",sep='') )
# write.csv(resultsOK,file=paste("finalResults/resultsOK_", outputN, ".csv",sep='') )
# write.csv(numGroups,file=paste("finalResults/numGroups_", outputN, ".csv",sep='') )
# write.csv(numSamples,file=paste("finalResults/numSamples_", outputN, ".csv",sep='') )


finalResults <- data.frame(matrix(NA, nrow = length(clustersK), ncol = 3 ))
row.names(finalResults) <- clustersK


finalResults[,1] <- as.numeric(results[1,])
finalResults[,2] <-  as.numeric(numGroups[1,])
finalResults[,3] <-  as.numeric(numSamples[1,])


write.csv(finalResults,file=paste("finalResults/finalTable_NBS", outputN, ".csv",sep='') )

View(finalResults)

minK <- min(finalResults[,1])
which.min(finalResults[,1])

