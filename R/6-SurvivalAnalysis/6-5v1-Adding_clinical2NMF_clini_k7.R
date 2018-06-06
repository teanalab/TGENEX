#version 1 (Jun/6/2018)

rm(list = ls())

load("temp/6-5v1.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)


#NMF clinical all components
patiAfi_NMF_CxP <- patientAfiliation_NMF_clini
coxFit_Afi_NMF_clin <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~factor(component),
                        x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NMF_CxP)
summary(coxFit_Afi_NMF_clin)
# Concordance= 0.686  (se = 0.052 )
# Rsquare= 0.026   (max possible= 0.606 )
# Likelihood ratio test= 12.16  on 6 df,   p=0.06
# Wald test            = 11.51  on 6 df,   p=0.07
# Score (logrank) test = 12.73  on 6 df,   p=0.05



#Top 5 phenotypes
patiAfi_NMF_CxP <- patientAfiliation_NMF_clini[which(as.numeric(patientAfiliation_NMF_clini$component) %in% bc4),]
coxFit_Afi_NMF_clin <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~factor(component),
                        x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NMF_CxP)
summary(coxFit_Afi_NMF_clin)
# exp(coef) exp(-coef) lower .95 upper .95
# factor(component)2     1.463     0.6833    0.3159     6.779
# factor(component)3     1.521     0.6574    0.2784     8.311
# factor(component)5     2.990     0.3344    0.6024    14.840
# factor(component)7     4.287     0.2333    0.9893    18.576



#add clinical variables

patiAfi_NMF_CxP_and_cliniVars <- merge( as.data.frame(patiAfi_NMF_CxP), as.data.frame(tensorClinical[row.names(patiAfi_NMF_CxP),]), by='row.names', all=TRUE)
patiAfi_NMF_CxP_and_cliniVars <- patiAfi_NMF_CxP_and_cliniVars[,-1]
#Stage IIIB has zero TRUE status
#solution:
#group all stage III together
patiAfi_NMF_CxP_and_cliniVars[grep("III",patiAfi_NMF_CxP_and_cliniVars$Converted.Stage),
                          "Converted.Stage"] = "Stage Three"
patiAfi_NMF_CxP_and_cliniVars[grep("II",patiAfi_NMF_CxP_and_cliniVars$Converted.Stage),
                          "Converted.Stage"] = "Stage Two"
table(patiAfi_NMF_CxP_and_cliniVars$Converted.Stage)


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
                            data=patiAfi_NMF_CxP_and_cliniVars)
summary(coxFit_plusClinical)


