# version 1 (1/24/2018)
# Likelihood test of survival curves

#########
#Here I compared two Cox models:
#1) with clinical variables only. Variables that did not gave me Inf betha : survivalClinical$Diagnosis.Age +
#       survivalClinical$HER2.Status +survivalClinical$Metastasis.Coded + survivalClinical$Tumor..T1.Coded
#2) Same variables from model 1 plus component number survivalClinical$component

# we compared the Cox models with a likelihood test and we concluded that these models
# are not significantly different  Pr(>Chisq) = 1 and Chisq = 0

# the log rank test of model 1 = p=5.09e-05
# the log rank test of model 2 = p=5.09e-05

# The clinical variable with lowest p-value is Tumor..T1.Coded equal to p=0.5128


# other tests I do
# 91_ testing with different samples due to stage
# 92_ comparing with other clusters
# 93_ reducing my components like I did in 8


rm(list = ls())
load(file="temp/9_likelihoodT.RData")


# #load my libraries
# libs<-c("Packages.R", "permutation.R", "clustering.R", "Plots.R", "survival_analysis.R")
# libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
# sapply(libs, function(u) {source(u)})
#

#loadlib("plyr", FALSE) #for function count
#loadlib("lmtest", FALSE) #for comparing two models - likelihood


#required for the cox proportional hazard model
loadls("plyr survival Rcpp survAUC", FALSE)


# #my functions
# source("R/plotSurvivalUpdate.R")

#loading and saving data
#####
# #load components by patient, generated with 41-AnalizeCP_BRCA
# load("results/BRCA15.RData")
#
# A <- Boolcp_3comp$A #patients vs components
# rm(MatrixOrderPxGxC)
# rm(Boolcp_3comp)
#
# patientVsComponents <- as.data.frame(A)
# patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
# patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)
#
# patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
#                                 component=patientVsComponents$whichmax)
#
# #frequency of components
# table1 = count(patientAfiliation, 'component')
#
# #attach survival data
# clinical<-read.table(file="../../12Datasets/BRCA/clinical.txt", sep = "\t",
#                      header = TRUE, stringsAsFactors = FALSE)
#
#
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
#
#
# survivalClinical = clinical[,-c(2,3,4,5,10,11,12,13,16,21,24,25,26)]
# survivalClinical = merge(patientAfiliation, survivalClinical, by.x="patient", by.y= "Patient.ID")
#
# table2 = count(survivalClinical, 'Overall.Survival.Status')
# table3 = count(survivalClinical, 'Person.Gender')
#
# groups = factor(survivalClinical$component)
#
#
# survivalClinical$Overall.Survival.Status <- (survivalClinical$Overall.Survival.Status=="DECEASED")
#




######

response = Surv(time = survivalClinical$Overall.Survival..Months.,
                event = survivalClinical$Overall.Survival.Status)
#with components
coxFit1 <- coxph(response~cluster(survivalClinical$component)+
                   survivalClinical$Diagnosis.Age +
                   survivalClinical$Converted.Stage + survivalClinical$ER.Status +
                   survivalClinical$HER2.Status +
                   survivalClinical$Metastasis.Coded + survivalClinical$PR.Status +
                    survivalClinical$Tumor..T1.Coded) #,
                #ties ="exact",method = "breslow")
summary(coxFit1)


#with no components
response2 = Surv(time = survivalClinical$Overall.Survival..Months.,
                 event = survivalClinical$Overall.Survival.Status)
#with components
coxFit2 <- coxph(response2~survivalClinical$Diagnosis.Age +
                   survivalClinical$Converted.Stage + survivalClinical$ER.Status +
                   survivalClinical$HER2.Status +
                   survivalClinical$Metastasis.Coded + survivalClinical$PR.Status +
                   survivalClinical$Tumor..T1.Coded) #,
#ties ="exact",method = "breslow")
summary(coxFit2)


######

response = Surv(time = survivalClinical$Overall.Survival..Months.,
                event = survivalClinical$Overall.Survival.Status)


#with components
coxFit1 <- coxph(response~ survivalClinical$Diagnosis.Age +
                 survivalClinical$HER2.Status +
                 survivalClinical$Metastasis.Coded +
                 survivalClinical$Tumor..T1.Coded,
                 ties ="exact",method = "breslow")
summary(coxFit1)


#with no components
coxFit2 <- coxph(response ~ survivalClinical$Diagnosis.Age +
                   survivalClinical$HER2.Status +
                   survivalClinical$Metastasis.Coded +
                   survivalClinical$Tumor..T1.Coded+
                   cluster(survivalClinical$component),
                 ties ="exact",method = "breslow")
summary(coxFit2)

#compare models
#https://stackoverflow.com/questions/43143757/likelihood-ratio-test-of-two-cox-models-with-multiple-readings-per-person

#log likelihood
lrtest(coxFit2,coxFit1)


######

plot(survfit(response~cluster(survivalClinical$Diagnosis.Age)))

plot(survfit(coxFit1))

plot(survfit(coxFit2))

#plotCompleteKplan(groups, survivalData, resultsFolder=getwd(), pdfName1="temp/75_kaplan.pdf")

sess <- sessionInfo() #save session on variable
save.image(file="temp/9_likelihoodT.RData")
