# version 2 (8/15/2018)
# run all algorithms
# with rank = 30



#######
rm(list = ls())

#load("data/binaMutation.Rd") #mutation matrix as integer
load("data/boolMutation.RData") #mutation matrix as bool


load("loadls.RData")  #load my libraries
loadls("plyr survival lmtest", T) #survival analysis
loadls("stats registry", T)  #random numbers
loadls("methods utils pkgmaker registry rngtools cluster", T) #NMF prerequisites
loadlib("NMF", T) #For non-negative matrix factorization

rankNMF <- 4
res <- nmf(boolMutation, rankNMF, "snmf/l", seed = 123456)


# get matrix W
w <- basis(res)
dim(w)
# get matrix H
h <- coef(res)
dim(h)
## [1] 3 38

########
for(i in c(1:rankNMF))
{
  w.sub <- w[,1:i]  #get the first i-th elements
  patientVsK <- as.data.frame(w.sub)
  patientVsK[, "max"] <- apply(patientVsK, 1, max)
  patientVsK[, "whichmax"] <- apply(patientVsK, 1, which.max)

  patientAfiliationNMF_PxM <- data.frame(patient=row.names(patientVsK),
                                  component=patientVsK$whichmax)
}

patientNMF_PxM <- w
mutationNMF_PxM <- h

save(patientAfiliationNMF_PxM, patientNMF_PxM, mutationNMF_PxM,
     file="data/NMF_k4_PxM.RData")

# clinical NMF -----------

load("data/boolClinical.RData")
res <- nmf(boolClinical, rankNMF, "snmf/l", seed = 123456)

# get matrix W
w <- basis(res)
dim(w)
# get matrix H
h <- coef(res)
dim(h) # [1] 15 21

for(i in c(1:rankNMF))
{
  w.sub <- w[,1:i]  #get the first i-th elements
  patientVsK <- as.data.frame(w.sub)
  patientVsK[, "max"] <- apply(patientVsK, 1, max)
  patientVsK[, "whichmax"] <- apply(patientVsK, 1, which.max)
  patientAfiliationNMF_PxC <- data.frame(patient=row.names(patientVsK),
                                     component=patientVsK$whichmax)
}
patientNMF_PxC <- w
mutationNMF_PxC <- h

save(patientAfiliationNMF_PxC, patientNMF_PxC, mutationNMF_PxC,
     file="data/NMF_4_PxC.RData")


sess <- sessionInfo() #save session on variable
save.image(file="temp/2-1v2_30_NMFs.RData")



###### old stuff -----

