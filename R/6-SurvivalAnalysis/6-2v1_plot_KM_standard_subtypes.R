#version 1 (Jun/5/2018)

rm(list = ls())

load("data/survivalClinical-5-1_30.RData")
load("data/myLib.RData")
source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")
loadls("plyr survival missForest survAUC prodlim survminer", F)


#PAM50
row.names(survClinical) <- survClinical$patient.bcr_patient_barcode
PAM50_surv_c <- survClinical[patients,c("patient.bcr_patient_barcode", "PAM50.Subtype")]
PAM50_surv_c <- merge(as.data.frame(PAM50_surv_c), as.data.frame(survivalClinical), by='row.names', all=TRUE)
PAM50_surv_c <- PAM50_surv_c[,-1]
coxFit_PAM50_surv_c <- coxph(Surv(time = Overall.Survival..Months.,
                                  event = Overall.Survival.Status)~factor(PAM50.Subtype),
                             x=T, y=T, model=T, method = "breslow",
                             data=PAM50_surv_c)
summary(coxFit_PAM50_surv_c)
# Concordance= 0.585  (se = 0.05 )
# Rsquare= 0.019   (max possible= 0.606 )
# Likelihood ratio test= 8.67  on 4 df,   p=0.07
# Wald test            = 8.07  on 4 df,   p=0.09
# Score (logrank) test = 8.6  on 4 df,   p=0.07

# table(PAM50_surv_c$PAM50.Subtype)

mfit_PAM50_surv_c <- survfit(Surv(time = Overall.Survival..Months.,
                                  event = Overall.Survival.Status)~factor(PAM50.Subtype),
                             data=PAM50_surv_c)

test <- summary(coxFit_PAM50_surv_c)
ggsurvplot(mfit_PAM50_surv_c, conf.int = F, pval = test$sctest[3])

PAM50_surv_c_no_normal <- which(PAM50_surv_c$PAM50.Subtype == "Normal-like")
coxFit_PAM50_surv_c_no_normal <- coxph(Surv(time = Overall.Survival..Months.,
                                            event = Overall.Survival.Status)~factor(PAM50.Subtype),
                                       x=T, y=T, model=T, method = "breslow",
                                       data=PAM50_surv_c[-PAM50_surv_c_no_normal,])

summary(coxFit_PAM50_surv_c_no_normal)

mfit_PAM50_no_normal <- survfit(Surv(time = Overall.Survival..Months.,
                                     event = Overall.Survival.Status)~factor(PAM50.Subtype),
                                data=PAM50_surv_c[-PAM50_surv_c_no_normal,])

test <- summary(coxFit_PAM50_surv_c_no_normal)
ggsurvplot(mfit_PAM50_no_normal, conf.int = F, pval = test$sctest[3])

table(PAM50_surv_c$PAM50.Subtype)

setEPS()
postscript("output4paper/KM_standard_subtypes.eps", width = 12, height = 9)
#tiff(file = "output4paper/KM_standard_subtypes.tiff",  width = 12, height = 9, units = 'in', res=300)
#pdf( file = "output4paper/KM_standard_subtypes.pdf",  onefile = TRUE, width = 12, height = 9)
plot.Survival4paper(coxFit_PAM50_surv_c_no_normal, mfit_PAM50_no_normal, location = "topright",
                    colorsP = c("red","black","blue","#1E8449" ),
                    colorsL = c("#1E8449","blue","black","red" ),
                    labelClu = c("Luminal B (110)", "Luminal A (200)",  "HER2-enriched (54)", "Basal-like (85)") )
# font_size_times = 0.5)
dev.off()
