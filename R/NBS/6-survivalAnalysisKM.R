rm(list = ls())


k = 5
#load(file=paste("dataExample/afiliations/affiliationListWithOutL_k",k,"_i",i,".RData",sep='' ) )
load(file=paste0("temp/affiliationNBSk",k,".RData") )

plotName = paste0("temp/kaplan-nbs_k",k,".pdf")


# #Figure 4.b) NBS
# k = 6
# load(file=paste("dataExample/afiliations/affiliationNBSk",k,".RData",sep='' ) )
# plotName = paste0("temp/kaplan-NBS_k",k,".pdf")


## different configurations ----
# rm(list = ls())
# k=3
# i=6
# plotName = paste0("temp/kaplan-cligen_k",k,"i",main.i,".pdf")
#k-means no OL
#load(file=paste("dataExample/afiliations/affiliationList_k",k,"_i",i,".RData",sep='' ) )
#k-means with OL
#load(file=paste("dataExample/afiliations/affiliationListWithOutL_k",k,"_i",i,".RData",sep='' ) )
#hierarchical
#load(file=paste("dataExample/afiliations/affiliationHC_k",k,"_i",i,".RData",sep='' ) )
#hierarchical
#load(file=paste("dataExample/afiliations/affiliationHCwithOL_k",k,"_i",i,".RData",sep='' ) )
########


#load(file="dataExample/Filtered/survivalClinical.RData")

#loat plotting functions
source("R/6-plotSurvivalUpdate.R")
#source('R/5-getPvalsSurvAnal.R') #get minimum p-val
source('R/5-survivalAnalysis.R')

#LUAD
#patients from NBS
patients <- names(affiliationList) #only 310 in survival file

source("R/1-ReadRawFiles.R")
#survival from RTCGA
survivalVars <- survivalFromRTCGA("LUAD")
survivalClinical <-survivalVars
row.names(survivalClinical) <- survivalClinical[,2]
survivalClinical <- survivalClinical[patients,]
survivalClinical$patient.vital_status  = survivalClinical$patient.vital_status == 1
survivalClinical <- survivalClinical[,-2]
names(survivalClinical) <- c( "Overall.Survival..Months.", "Overall.Survival.Status" )
# FOR LUAD
survivalClinical$Overall.Survival..Months. <- (survivalClinical$Overall.Survival..Months.*12/356)

#remove NA
if(anyNA(survivalClinical)) {
  remNA = which(is.na(survivalClinical$Overall.Survival..Months.))
  survivalClinical <- survivalClinical[-remNA,]
}



threshold_groups = 10 #10 or 3

list[pval, ok, patiF, nFactors] =  survivalAnalysis('-', k, threshold_groups, plots=1, plotName, location = "topright")

#survAna(k,0)
