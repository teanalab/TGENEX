#version 1 (Jun/6/2018)

rm(list = ls())

load("temp/6-6v1.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)


#NMF clinical all components
patiAfi_NMF_MxP <- patientAfiliation_NMF_mut
coxFit_Afi_NMF_mut <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~factor(component),
                        x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NMF_MxP)
summary(coxFit_Afi_NMF_mut)
# Concordance= 0.59  (se = 0.044 )
# Rsquare= 0.013   (max possible= 0.606 )
# Likelihood ratio test= 6.08  on 5 df,   p=0.3
# Wald test            = 7.22  on 5 df,   p=0.2
# Score (logrank) test = 8.79  on 5 df,   p=0.1


#Top 5 phenotypes
patiAfi_NMF_MxP <- patientAfiliation_NMF_mut[which(as.numeric(patientAfiliation_NMF_mut$component) %in% bc4),]
coxFit_Afi_NMF_mut <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~factor(component),
                        x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NMF_MxP)
summary(coxFit_Afi_NMF_mut)
# exp(coef) exp(-coef) lower .95 upper .95
# factor(component)2 2.958e-08  3.381e+07   0.00000       Inf
# factor(component)3 6.779e-01  1.475e+00   0.05730     8.020
# factor(component)4 1.603e-01  6.238e+00   0.02070     1.241
# factor(component)5 3.579e-01  2.794e+00   0.02177     5.884
#
# Concordance= 0.537  (se = 0.029 )
# Rsquare= 0.016   (max possible= 0.542 )
# Likelihood ratio test= 5.31  on 4 df,   p=0.3
# Wald test            = 6.97  on 4 df,   p=0.1
# Score (logrank) test = 8.86  on 4 df,   p=0.06



#add clinical variables
patiAfi_NMF_MxP_and_cliniVars <- merge( as.data.frame(patiAfi_NMF_MxP), as.data.frame(tensorClinical[row.names(patiAfi_NMF_MxP),]), by='row.names', all=TRUE)
patiAfi_NMF_MxP_and_cliniVars <- patiAfi_NMF_MxP_and_cliniVars[,-1]
#Stage IIIB has zero TRUE status
#solution:
#group all stage III together
patiAfi_NMF_MxP_and_cliniVars[grep("III",patiAfi_NMF_MxP_and_cliniVars$Converted.Stage),
                          "Converted.Stage"] = "Stage Three"
patiAfi_NMF_MxP_and_cliniVars[grep("II",patiAfi_NMF_MxP_and_cliniVars$Converted.Stage),
                          "Converted.Stage"] = "Stage Two"
table(patiAfi_NMF_MxP_and_cliniVars$Converted.Stage)


##Rerun Analysis
coxFit_plusClinical <- coxph(Surv(time = Overall.Survival..Months.,
                                 event = Overall.Survival.Status)~
                              factor(component)+
                              #patient.age_at_initial_pathologic_diagnosis+
                              Node.Coded+
                              Tumor+
                              HER2.Final.Status+
                              patient.breast_carcinoma_estrogen_receptor_status+
                              patient.breast_carcinoma_progesterone_receptor_status,
                            x=T, y=T, model=T, method = "breslow",
                            data=patiAfi_NMF_MxP_and_cliniVars)
summary(coxFit_plusClinical)
# Warning message:
#   In fitter(X, Y, strats, offset, init, control, weights = weights,  :
#               Ran out of iterations and did not converge

with(patiAfi_NMF_MxP_and_cliniVars, table(patient.breast_carcinoma_progesterone_receptor_status,
                                          Overall.Survival.Status))
#indeterminate PR has zero TRUE status
#solution:
#remove indeterminate
patiAfi_NMF_MxP_and_cliniVars<-patiAfi_NMF_MxP_and_cliniVars[-which(patiAfi_NMF_MxP_and_cliniVars$patient.breast_carcinoma_progesterone_receptor_status == "indeterminate"), ]



##Rerun Analysis
coxFit_plusClinical <- coxph(Surv(time = Overall.Survival..Months.,
                                  event = Overall.Survival.Status)~
                               factor(component)+
                               #patient.age_at_initial_pathologic_diagnosis+
                               Node.Coded+
                               Tumor+
                               HER2.Final.Status+
                               patient.breast_carcinoma_estrogen_receptor_status+
                               patient.breast_carcinoma_progesterone_receptor_status,
                             x=T, y=T, model=T, method = "breslow",
                             data=patiAfi_NMF_MxP_and_cliniVars)
summary(coxFit_plusClinical)
