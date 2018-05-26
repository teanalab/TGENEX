#version 1 (1/29/2018)
#matrix factorization at different k

#Cox p values for K
# 0 1 1
# 0.1028745 1 0.748407
# 0.8084114 1 0.3685903
# 0.556556 1 0.4556512
# 0.04287207 1 0.8359663
# 0.2239175 1 0.6360711
# 1.912268 1 0.1667115
# 4.242531 1 0.03942336
# 3.287555 1 0.06980691
# 7.852071 1 0.005076252   ---- best 10
# 1.673628 1 0.1957733
# 2.531754 1 0.1115761
# 2.845336 1 0.09163946
# 3.206856 1 0.07333025
# 0.7705026 1 0.3800616


#######
rm(list = ls())
load("results/10March2018/10-BoolMatrices_v2.RData")

x=ls()[-10] #keep boolClinical only
rm(list = x)
rm(x)


#save.image(file="temp/100_NMFsurv_cli.RData")

# #load my libraries
libs<-c("Packages.R", "survival_analysis.R")
libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
sapply(libs, function(u) {source(u)})
install.packages("digest", repos="http://R-Forge.R-project.org") #required

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

res <- nmf(boolClinical,rank=10,"snmf/l",seed = 123456)

# get matrix W
w <- basis(res)
dim(w)

save(w,file="temp/w_clinical.RData")
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

for(i in c(1:10))
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

for(i in c(1:10))
{
  cat('\n',listSurv[[i]]$sctest)
}


sess <- sessionInfo() #save session on variable
save.image(file="temp/100_NMFsurv_clinical.RData")
load(file="temp/100_NMFsurv_clinical.RData")
