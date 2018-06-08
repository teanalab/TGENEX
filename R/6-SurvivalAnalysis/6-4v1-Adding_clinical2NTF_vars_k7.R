# #version 1 (Jun/6/2018)
#
# rm(list = ls())
#
#
# load("temp/6-4v1.RData")
# loadls("plyr survival missForest survAUC prodlim survminer", F)
#
# patiAfi_NTF <- patientAfiliation_NTF
#
# # NFT all components
# # coxFit_Afi_NTF <- coxph(Surv(time = Overall.Survival..Months.,
# #                             event = Overall.Survival.Status)~factor(component),
# #                        x=T, y=T, model=T, method = "breslow",
# #                        data=patiAfi_NTF)
# # summary(coxFit_Afi_NTF)
# # Concordance= 0.632  (se = 0.052 )
# # Rsquare= 0.022   (max possible= 0.606 )
# # Likelihood ratio test= 10.06  on 6 df,   p=0.1
# # Wald test            = 10.79  on 6 df,   p=0.1
# # Score (logrank) test = 11.84  on 6 df,   p=0.07
#
#
# # top 5 components
# patiAfi_NTF <- patientAfiliation_NTF[which(as.numeric(patientAfiliation_NTF$component) %in% bc4),]
#
# # coxFit_Afi_NTF <- coxph(Surv(time = Overall.Survival..Months.,
# #                             event = Overall.Survival.Status)~factor(component),
# #                        x=T, y=T, model=T, method = "breslow",
# #                        data=patiAfi_NTF)
# # summary(coxFit_Afi_NTF)
# # Concordance= 0.63  (se = 0.055 )
# # Rsquare= 0.027   (max possible= 0.613 )
# # Likelihood ratio test= 9.16  on 4 df,   p=0.06
# # Wald test            = 9.24  on 4 df,   p=0.06
# # Score (logrank) test = 10.16  on 4 df,   p=0.04
#
#
# patiAfi_NTF_and_cliniVars <- merge( as.data.frame(patiAfi_NTF),as.data.frame(tensorClinical[row.names(patiAfi_NTF),]), by='row.names', all=TRUE)
# patiAfi_NTF_and_cliniVars <- patiAfi_NTF_and_cliniVars[,-1]
# patiAfi_NTF_and_cliniVars$age <- as.numeric(as.character(patiAfi_NTF_and_cliniVars$patient.age_at_initial_pathologic_diagnosis))
# patiAfi_NTF_and_cliniVars$age.f <- cut(patiAfi_NTF_and_cliniVars$age, c(28,50,70,85), ordered_result = TRUE)
# patiAfi_NTF_and_cliniVars$Converted.Tumor <- with(patiAfi_NTF_and_cliniVars, ifelse(patiAfi_NTF_and_cliniVars$Tumor == "T1", "T1", "T_other"))
#
# # Stage IIIB has zero TRUE status
# # solution: group all stage III together
# patiAfi_NTF_and_cliniVars[grep("I",patiAfi_NTF_and_cliniVars$AJCC.Stage),
#                           "Normalized.Stage"] = "I"
# patiAfi_NTF_and_cliniVars[grep("II",patiAfi_NTF_and_cliniVars$AJCC.Stage),
#                           "Normalized.Stage"] = "II"
# patiAfi_NTF_and_cliniVars[grep("III",patiAfi_NTF_and_cliniVars$AJCC.Stage),
#                           "Normalized.Stage"] = "III"
# patiAfi_NTF_and_cliniVars[grep("IV",patiAfi_NTF_and_cliniVars$AJCC.Stage),
#                           "Normalized.Stage"] = "IV"
#
# # variable stats
# table(patiAfi_NTF_and_cliniVars$Normalized.Stage)
# with(patiAfi_NTF_and_cliniVars, table(Overall.Survival.Status, patient.histological_type))
#
# coxFit_tensorClini <- coxph(Surv(time = Overall.Survival..Months.,
#                                  event = Overall.Survival.Status)~
#                               factor(component)+
#                               # age.f +
#                               factor(Normalized.Stage) +
#                               factor(Node.Coded) +
#                               #factor(Metastasis.Coded)
#                               factor(Converted.Tumor)#+
#                               # factor(patient.histological_type)+
#                               # factor(HER2.Final.Status)+
#                               # factor(patient.breast_carcinoma_estrogen_receptor_status)+
#                               # factor(patient.breast_carcinoma_progesterone_receptor_status)
#                               ,x=T, y=T, model=T, method = "breslow",
#                               data=patiAfi_NTF_and_cliniVars)
# summary(coxFit_tensorClini)
