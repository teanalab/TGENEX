rm(list = ls())
load("data/loadls.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)
load(file="dataExample/Filtered/survivalClinical.RData")

results <- rep(0,7)
names(results) <- c(4:10)

# resultsOK <- data.frame(matrix(NA, nrow = 14, ncol = 9))
# row.names(resultsOK) <- c(2:15)
# names(resultsOK) <- c(2:10)

# binaMredu <-read.table(file="dataExample/fromNBS/10-binaM.csv",sep = ";", quote = '"',
#                        header = TRUE, stringsAsFactors = FALSE)
# standardNMF <-read.table(file="dataExample/fromNBS/standardNMF_k4.csv", sep = ",", quote = "",
#                          header = FALSE, stringsAsFactors = FALSE)
# genes_redu <-read.table(file="dataExample/fromNBS/genes.csv", sep = "\t", quote = '"',
#                         header = TRUE, stringsAsFactors = FALSE)
# clinical <-read.table(file="dataExample/fromNBS/10-binac.csv",sep = ";", quote = '"',
#                       header = TRUE, stringsAsFactors = FALSE)
load("dataExample/Filtered/patients.RData")

for(k in names(results) )
{
  standardNMF <-read.table(file=paste0("dataExample/fromNBS/standardNMF_k",k,".csv"), sep = ",", quote = "",
                           header = FALSE, stringsAsFactors = FALSE)
  
  affiliationList = standardNMF
  affiliationList = unlist(affiliationList)
  names(affiliationList) = as.character(unlist(patients))
  save(affiliationList, file = paste0("dataExample/afiliations/affiliationNMFk",k,".RData") ) 
  
  patiF <- data.frame(patients,t(standardNMF) )
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2)]
  
  #remove NA
  #patiF <- patiF[-c(1,2),]
  if(anyNA(patiF)){
    remNA = which(is.na(patiF$Overall.Survival.Status))
    patiF <- patiF[-remNA,]
  }
  
  # #Kaplan Meier all samples
  # plot(survfit(coxFit_NBS), main = "NBS")
  
  #NBS components
  SurFunction <- Surv(time = Overall.Survival..Months.,
                      event = Overall.Survival.Status)~.
  coxFit_NBS <- coxph(SurFunction,
                      x=T, y=T, model=T, method = "breslow",
                      data=patiF)
  summary(coxFit_NBS)
  test <- summary(coxFit_NBS)
  results[k] <- test$sctest[3]
  rm(list = ls()[-which(ls() %in% c("k" ,"results", "survivalClinical", "patients")  )])
}

