rm(list = ls())

# > min(minR)
# [1] 3.686312e-07
# > which.min(minK)
# 15 
# 14 
# > which.min(minR)
# 62 
# 60 

#Figure 4.a) best subtypes obtained with cligen: Kmeans with outliers
k = 9
r = 27
#load(file=paste("dataExample/afiliations/affiliationListWithOutL_k",k,"_i",i,".RData",sep='' ) )
load(file=paste("/Users/diamac/CLIGEN_server/temp/temp/affiListKmeans_r",r,"_k",k,".RData",sep=''))

plotName = paste0("temp/kaplan-nmf_k",k,"r", r,".pdf")


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
load("data/survivalVars.RData")
load("data/patients.RData")
survivalClinical <-survivalVars
survivalClinical <- survivalClinical[patients,]
survivalClinical$Vital_Status = survivalClinical$Vital_Status != "living"
names(survivalClinical) <- c( "Overall.Survival..Months.", "Overall.Survival.Status" )
survivalClinical$Overall.Survival..Months. <- (survivalClinical$Overall.Survival..Months.*12/356)


threshold_groups = 3 #10 or 3

list[pval, ok, patiF, nFactors] =  survivalAnalysis(r, k, threshold_groups, plots=1, plotName, location = "bottomright")

#survAna(k,0)
