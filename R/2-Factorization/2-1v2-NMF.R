# version 2 (5/26/2018)
# run all algorithms
# with rank = 30



#######
rm(list = ls())


#load my libraries
libs<-c("Packages.R", "survival_analysis.R")
libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
sapply(libs, function(u) {source(u)})

save.image("data/BC/myLib.RData")

load("data/BC/boolMutation.RData")

loadls("plyr survival lmtest", T)
loadls("stats registry", T)
loadls("methods utils pkgmaker registry rngtools cluster", T)
loadlib("NMF", T) #For non-negative matrix factorization
#The paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-367
#The vignette: https://cran.r-project.org/web/packages/NMF/index.html

##List of algorithms on the package
# nmfAlgorithm()
# ##  [1] "brunet"    "KL"        "lee"
# ##  [6] "nsNMF"     "ls-nmf"    "pe-nmf"
# ## [11] "snmf/l"

# meth <- nmfAlgorithm(version = "R")
# meth <- c(names(meth), meth)
# meth

rankNMF <- 10

res <- nmf(boolMutation, rankNMF, "snmf/l", seed = 123456)




# get matrix W
w <- basis(res)
dim(w)

# get matrix H
h <- coef(res)
dim(h)
## [1] 3 38

table(survivalClinical$Overall.Survival.Status)


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
     file="data/NMF_30_PxM.RData")

# clinical NMF -----------

load("data/boolClinical.RData")

res <- nmf(boolClinical, rankNMF, "snmf/l", seed = 123456)

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

  patientAfiliationNMF_PxC <- data.frame(patient=row.names(patientVsK),
                                     component=patientVsK$whichmax)
}


patientNMF_PxC <- w
mutationNMF_PxC <- h

save(patientAfiliationNMF_PxC, patientNMF_PxC, mutationNMF_PxC,
     file="data/NMF_30_PxC.RData")


sess <- sessionInfo() #save session on variable
save.image(file="temp/2-1v2_30_NMFs.RData")




###### old stuff -----

