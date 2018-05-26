#version 1 (1/29/2018)
#run all algorithms


#######
rm(list = ls())
load("data/boolMutation.RData")

#save.image(file="temp/101_NMFs.RData")

# #load my libraries
libs<-c("Packages.R", "survival_analysis.R")
libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
sapply(libs, function(u) {source(u)})
##install.packages("digest", repos="http://R-Forge.R-project.org") #required

loadls("plyr survival lmtest", FALSE)
loadls("stats registry", FALSE)
loadls("methods utils pkgmaker registry rngtools cluster", FALSE)
loadlib("NMF", FALSE) #For non-negative matrix factorization
#The paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-367
#The vignette: https://cran.r-project.org/web/packages/NMF/index.html


##List of algorithms on the package
# nmfAlgorithm()
# ##  [1] "brunet"    "KL"        "lee"
# ##  [6] "nsNMF"     "ls-nmf"    "pe-nmf"
# ## [11] "snmf/l"



meth <- nmfAlgorithm(version = "R")
meth <- c(names(meth), meth)
meth


res <- nmf(boolMutation, 15, meth, seed = 123456)

res <- nmf(boolMutation, 15, "snmf/l", seed = 123456)

# get matrix W
w <- basis(res)
dim(w)

# get matrix H
h <- coef(res)
dim(h)
## [1] 3 38


########
#get survival data
clinical<-read.table(file="../../12Datasets/BRCA/clinical.txt", sep = "\t",
                     header = TRUE, stringsAsFactors = FALSE)
# #remove unwanted columns
# #[2] "Sample.ID"
# #[3] "CN.Cluster"                       "Cancer.Studies"
# #[5] "Cancer.Type.Detailed"
# #[10] "Integrated.Clusters..no.exp."
# #[11] "Integrated.Clusters..unsup.exp."  "Integrated.Clusters..with.PAM50."
# #[13] "MIRNA.Cluster"
# #[16] "Methylation.Cluster"
# #[21] "PAM50.subtype"
# #[24] "RPPA.Cluster"
# #[25] "SigClust.Intrinsic.mRNA"          "SigClust.Unsupervised.mRNA"
#
#
# #include survival columns
# #[19] "Overall.Survival..Months."        "Overall.Survival.Status"
# #[27] "Survival.Data.Form"

survivalClinical = clinical[,c(1,19,20)]
table2 = count(survivalClinical, 'Overall.Survival.Status')
survivalClinical$Overall.Survival.Status <- (survivalClinical$Overall.Survival.Status=="DECEASED")
table3 = count(survivalClinical, 'Overall.Survival.Status')

######


#get membership
#attach survival data
#do survival analysis

listSurv <- c()

for(i in c(1:15))
{
  w.sub <- w[,1:i]  #get the first i-th elements
  patientVsK <- as.data.frame(w.sub)
  patientVsK[, "max"] <- apply(patientVsK, 1, max)
  patientVsK[, "whichmax"] <- apply(patientVsK, 1, which.max)

  patientAfiliation <- data.frame(patient=row.names(patientVsK),
                                  component=patientVsK$whichmax)
  survivalCl.sub = merge(patientAfiliation, survivalClinical, by.x="patient", by.y= "Patient.ID")

  response = Surv(time = survivalCl.sub$Overall.Survival..Months.,
                  event = survivalCl.sub$Overall.Survival.Status)
  coxFit1 <- coxph(response~survivalCl.sub$component)
  listSurv <- c(listSurv, list(summary(coxFit1)))
}

for(i in c(1:15))
{
  cat('\n',listSurv[[i]]$sctest)
}


sess <- sessionInfo() #save session on variable
save.image(file="temp/101_NMFs.RData")
