# version 1 (1/22/2018)
# AUC comparing model without components versus components

#Using all the clinical variables

rm(list = ls())


#load my libraries
libs<-c("Packages.R", "permutation.R", "clustering.R", "Plots.R", "survival_analysis.R")
libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
sapply(libs, function(u) {source(u)})


loadls("plyr") #for function count

#required for the cox proportional hazard model
loadls("survival Rcpp missForest survAUC")


#my functions
source("R/plotSurvivalUpdate.R")

#loading and saving data
#####
# #load components by patient, generated with 41-AnalizeCP_BRCA
# load("results/BRCA15.RData")
#
#
# A <- Boolcp_3comp$A #patients vs components
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
# # survivalData = data.frame(PatientID=survivalClinical$patient,
# #                           Survival=survivalClinical$Overall.Survival..Months.,
# #                           Death=(survivalClinical$Overall.Survival.Status=="DECEASED"))
#
# #plot.SurvivalKN (groups, "mainTitle", survivalData)
#
# save.image(file="temp/81_components.RData")

#####
load("temp/81_components.RData")

######

# set training and testing subsets (https://stackoverflow.com/questions/17200114/how-to-split-data-into-training-testing-sets-using-sample-function)
## 75% of the sample size
smp_size <- floor(0.75 * nrow(survivalClinical))

## set the seed to make your partition reproductible
set.seed(123)
train_ind <- sample(seq_len(nrow(survivalClinical)), size = smp_size)

train <- survivalClinical[train_ind, ]
test <- survivalClinical[-train_ind, ]

groupsTrain <- groups[train_ind]
groupsTest <- groups[-train_ind]

#cox proportional hazard model

coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                     event = Overall.Survival.Status)~Nodes.Details,
                x=TRUE, y=TRUE, method = "breslow", data=train)
summary(coxFit)

# Concordance= 0.668  (se = 0.056 )
# Rsquare= 0.011   (max possible= 0.639 )
# Likelihood ratio test= 4.19  on 21 df,   p=1
# Wald test            = 5.59  on 21 df,   p=0.9997
# Score (logrank) test = 194.5  on 21 df,   p=0


coxFit <- coxph(Surv(time = Overall.Survival..Months., event = Overall.Survival.Status)~Converted.Stage+
                  Diagnosis.Age+ER.Status+HER2.Status+Metastasis.Coded+PR.Status+
                  Tumor..T1.Coded+groupsTrain,
                x=TRUE, y=TRUE, method = "breslow", data=train)
summary(coxFit)


# Concordance= 0.721  (se = 0.056 )
# Rsquare= 0.054   (max possible= 0.639 )
# Likelihood ratio test= 21.21  on 35 df,   p=0.9679
# Wald test            = 12.2  on 35 df,   p=0.9999
# Score (logrank) test = 225.5  on 35 df,   p=0


#AUC
lp <- predict(coxFit)
lpnew <- predict(coxFit, newdata=test)
Surv.rsp <- Surv(time = train$Overall.Survival..Months., event = train$Overall.Survival.Status)
Surv.rsp.new <- Surv(time = test$Overall.Survival..Months., event = test$Overall.Survival.Status)
times <- seq(0.74, 221, 0.74)
AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
print(AUC_CD$iauc)
plot(AUC_CD)
#write.csv(cbind(AUC_CD$times, AUC_CD$auc), file = "c:/cancerstudy-rprogram/data/extraction/tauc-exp111.csv")
#write.csv(cbind(AUC_CD$times, AUC_CD$auc), file = "G:/teanalab/cancer-study-rprogramming/data/extraction/tauc-exp11.csv")

# #c-statistics
# cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
# print(cstat)
#
# #plot Kaplan-Meier curve
# kmsurvival <- survfit(Surv(train$obsPeriod,train$label)~1)
# plot(kmsurvival, ylim = c(.7,1), xlab="Time (days)", ylab="Proportion of patients surviving")



