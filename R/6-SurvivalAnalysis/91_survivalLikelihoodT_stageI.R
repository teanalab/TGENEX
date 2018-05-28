# version 1 (1/25/2018)
# Likelihood test of survival curves for stage I patients only

#########
# From script 9:
# Model 1) with clinical variables only. Variables that did not gave me Inf betha : survivalClinical$Diagnosis.Age +
#       survivalClinical$HER2.Status +survivalClinical$Metastasis.Coded + survivalClinical$Tumor..T1.Coded
# Model 2) Same variables from model 1 plus component number survivalClinical$component

# the previous test was done with the function cluster() around 'component'
# if we test without cluster the results are:
# Model 1 (Without)= 5.09e-05
# Model 2 (With component)= 8.004e-05
# Chisq = NA 0.7646354
# Pr(>Chisq) = NA 0.3818818
# No significant difference

# So I m not going to use cluster for these tests

# These are the frequencies of stages:
# > table4
# Converted.Stage freq
# 1                    1
# 2   No_Conversion  166
# 3         Stage I   78
# 4       Stage IIA  145
# 5       Stage IIB   54
# 6      Stage IIIA   41
# 7      Stage IIIB    4
# 8      Stage IIIC   14
# 9        Stage IV    4

# Test 1) there is a patient with missing data. What if I remove that patient?
# Model 2 gets worst but not significantly yet
# Model 1 (Without)= 0.0007819
# Model 2 (With component)= 0.001209
# Chisq = NA 0.5564876
# Pr(>Chisq) = NA 0.4556789


# Test 2) there are 5 males in the dataset. If I remove those 5 males...
# little worst both models no biggy
# Model 1 (Without)= 5.664e-05
# Model 2 (With component)= 8.6e-05
# Chisq = NA 0.8349007
# Pr(>Chisq) = NA 0.3608593

# Test 3) Using only samples from stage I (most difficult to predict)
# a lot worst both models, mine still the worst
# Model 1 (Without)= 0.8642
# Model 2 (With component)= 0.9254
# Chisq = NA 0.1274886
# Pr(>Chisq) = NA 0.7210503

# Test 4) Using all samples except for stage 4 (Metastasis)
# same old
# Model 1 (Without)= 3.596e-05
# Model 2 (With component)= 5.737e-05
# Chisq = NA 0.7686998
# Pr(>Chisq) = NA 0.3806196

# Test 5) Using only samples from stage= 'No_Conversion'
# worst both models, mine still the worst
# Model 1 (Without)= 0.06021
# Model 2 (With component)= 0.1217
# Chisq = NA 0.6729997
# Pr(>Chisq) = NA 0.4120077





########

rm(list = ls())
# load(file="temp/91_likelihoodT.RData")

# # #load my libraries
# libs<-c("Packages.R", "permutation.R", "clustering.R", "Plots.R", "survival_analysis.R")
# libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
# sapply(libs, function(u) {source(u)})
# loadls("survival lmtest plyr")

#loading and saving data
#######
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
# survivalClinical$Overall.Survival.Status <- (survivalClinical$Overall.Survival.Status=="DECEASED")
# survivalClinicalBackUp = survivalClinical



#selecting and describing stage I samples
######
#define models and compare them
testCoxModel <- function()
{
  response = Surv(time = survivalClinical$Overall.Survival..Months.,
                  event = survivalClinical$Overall.Survival.Status)

  #with components
  coxFit1 <- coxph(response~ survivalClinical$Diagnosis.Age +
                     survivalClinical$HER2.Status +
                     survivalClinical$Metastasis.Coded +
                     survivalClinical$Tumor..T1.Coded,
                   ties ="exact",method = "breslow")
  #summary(coxFit1)
  resCox <- mySumm(coxFit1)
  pVal<-format.pval(resCox$sctest[3], digits = 4)
  cat("Model 1 (Without)= ")
  cat(pVal)

  #with no components
  coxFit2 <- coxph(response ~ survivalClinical$Diagnosis.Age +
                     survivalClinical$HER2.Status +
                     survivalClinical$Metastasis.Coded +
                     survivalClinical$Tumor..T1.Coded+
                     survivalClinical$component,
                   ties ="exact",method = "breslow")
  #summary(coxFit2)
  resCox <- mySumm(coxFit2)
  pVal<-format.pval(resCox$sctest[3], digits = 4)
  cat("\nModel 2 (With component)= ")
  cat(pVal)

  #compare models
  #https://stackoverflow.com/questions/43143757/likelihood-ratio-test-of-two-cox-models-with-multiple-readings-per-person
  #log likelihood
  #lrtest(coxFit1,coxFit2)
  ans1=lrtest(coxFit1,coxFit2)
  cat("\nChisq = ")
  cat(ans1$Chisq)
  cat("\nPr(>Chisq) = ")
  cat(ans1$`Pr(>Chisq)`)
}



#####
#define models and compare them
testCoxModelOnlyAge <- function()
{
  response = Surv(time = survivalClinical$Overall.Survival..Months.,
                  event = survivalClinical$Overall.Survival.Status)

  #with components
  coxFit1 <- coxph(response~ survivalClinical$Diagnosis.Age, # +
                   #survivalClinical$HER2.Status + Inf
                   #survivalClinical$Metastasis.Coded +  2 or more levels
                   #survivalClinical$Tumor..T1.Coded, 2 or more levels
                   ties ="exact",method = "breslow")
  #summary(coxFit1)
  resCox <- mySumm(coxFit1)
  pVal<-format.pval(resCox$sctest[3], digits = 4)
  cat("Model 1 (Without)= ")
  cat(pVal)

  #with no components
  coxFit2 <- coxph(response ~ survivalClinical$Diagnosis.Age +
                     survivalClinical$component,
                   ties ="exact",method = "breslow")
  #summary(coxFit2)
  resCox <- mySumm(coxFit2)
  pVal<-format.pval(resCox$sctest[3], digits = 4)
  cat("\nModel 2 (With component)= ")
  cat(pVal)

  #compare models
  #https://stackoverflow.com/questions/43143757/likelihood-ratio-test-of-two-cox-models-with-multiple-readings-per-person
  #log likelihood
  #lrtest(coxFit1,coxFit2)
  ans1=lrtest(coxFit1,coxFit2)
  cat("\nChisq = ")
  cat(ans1$Chisq)
  cat("\nPr(>Chisq) = ")
  cat(ans1$`Pr(>Chisq)`)
}





######

# Test 1
# table4 = count(survivalClinical, 'Converted.Stage')
# which(survivalClinical$Converted.Stage=='') #remove this patient with missing data
# survivalClinical = survivalClinical[-415,]
# testCoxModel()

# #Test 2
# survivalClinical = survivalClinicalBackUp
# table5 = count(survivalClinical, 'Person.Gender')
# #removed 5 male from table
# survivalClinical = survivalClinical[-which(survivalClinical$Person.Gender=='MALE'),]
# testCoxModel()

# #Test 3
# survivalClinical = survivalClinicalBackUp
# survivalClinical = survivalClinical[which(survivalClinical$Converted.Stage=='Stage I'),]
# testCoxModelOnlyAge()


# #Test 4
# survivalClinical = survivalClinicalBackUp
# survivalClinical = survivalClinical[-which(survivalClinical$Converted.Stage=='Stage IV'),]
# testCoxModel()


# #Test 5
# survivalClinical = survivalClinicalBackUp
# survivalClinical = survivalClinical[which(survivalClinical$Converted.Stage=='No_Conversion'),]
# testCoxModelOnlyAge()


##############
plot(survfit(response~cluster(survivalClinical$Diagnosis.Age)))

plot(survfit(coxFit1))

plot(survfit(coxFit2))

#plotCompleteKplan(groups, survivalData, resultsFolder=getwd(), pdfName1="temp/75_kaplan.pdf")

sess <- sessionInfo() #save session on variable
save.image(file="temp/91_likelihoodT.RData")
