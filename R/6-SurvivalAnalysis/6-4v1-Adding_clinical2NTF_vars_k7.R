#version 1 (Jun/6/2018)

rm(list = ls())


load("temp/6-4v1.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)

patiAfi_NTF <- patientAfiliation_NTF
#NFT all components
coxFit_Afi_NTF <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~factor(component),
                        x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF)
summary(coxFit_Afi_NTF)
# Concordance= 0.632  (se = 0.052 )
# Rsquare= 0.022   (max possible= 0.606 )
# Likelihood ratio test= 10.06  on 6 df,   p=0.1
# Wald test            = 10.79  on 6 df,   p=0.1
# Score (logrank) test = 11.84  on 6 df,   p=0.07


#NFT top 5 components
patiAfi_NTF <- patientAfiliation_NTF[which(as.numeric(patientAfiliation_NTF$component) %in% bc4),]
coxFit_Afi_NTF <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~factor(component),
                        x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF)
summary(coxFit_Afi_NTF)
# Concordance= 0.63  (se = 0.055 )
# Rsquare= 0.027   (max possible= 0.613 )
# Likelihood ratio test= 9.16  on 4 df,   p=0.06
# Wald test            = 9.24  on 4 df,   p=0.06
# Score (logrank) test = 10.16  on 4 df,   p=0.04


patiAfi_NTF_and_cliniVars <- merge( as.data.frame(patiAfi_NTF),as.data.frame(tensorClinical[row.names(patiAfi_NTF),]), by='row.names', all=TRUE)
patiAfi_NTF_and_cliniVars <- patiAfi_NTF_and_cliniVars[,-1]
coxFit_tensorClini <- coxph(Surv(time = Overall.Survival..Months.,
                                 event = Overall.Survival.Status)~
                             factor(component)+
                              as.numeric(as.character(patient.age_at_initial_pathologic_diagnosis))+
                              Node.Coded+
                              Tumor+
                              HER2.Final.Status+
                              patient.breast_carcinoma_estrogen_receptor_status+
                              patient.breast_carcinoma_progesterone_receptor_status,
                            x=T, y=T, model=T, method = "breslow",
                            data=patiAfi_NTF_and_cliniVars)
summary(coxFit_tensorClini)
#WARNING: coxph ran out of iterations and did not converge because:
with(tensorClini, table(Overall.Survival.Status, Converted.Stage))
#Stage IIIB has zero TRUE status
#solution:
#group all stage III together

patiAfi_NTF_and_cliniVars[grep("III",patiAfi_NTF_and_cliniVars$Converted.Stage),
                          "Converted.Stage"] = "Stage Three"
patiAfi_NTF_and_cliniVars[grep("II",patiAfi_NTF_and_cliniVars$Converted.Stage),
                          "Converted.Stage"] = "Stage Two"
table(patiAfi_NTF_and_cliniVars$Converted.Stage)


##Rerun Analysis
coxFit_tensorClini <- coxph(Surv(time = Overall.Survival..Months.,
                                 event = Overall.Survival.Status)~
                              factor(component)+
                              patient.age_at_initial_pathologic_diagnosis+
                              Node.Coded+
                              Tumor+
                              HER2.Final.Status+
                              patient.breast_carcinoma_estrogen_receptor_status+
                              patient.breast_carcinoma_progesterone_receptor_status,
                            x=T, y=T, model=T, method = "breslow",
                            data=patiAfi_NTF_and_cliniVars)
summary(coxFit_tensorClini)


