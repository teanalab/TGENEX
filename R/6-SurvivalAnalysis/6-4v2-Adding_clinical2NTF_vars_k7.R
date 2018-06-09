#version 2 (Jun/8/2018)
# Adding AUC

rm(list = ls())


load("temp/6-4v2.RData")
load("data/myLib.RData")
loadls("plyr survival missForest survAUC prodlim survminer perry", F)


# ## Adding clinical variables ----
# patiAfi_NTF <- patientAfiliation_NTF
# patiAfi_NTF_and_cliniVars <- merge( as.data.frame(patiAfi_NTF),as.data.frame(tensorClinical[row.names(patiAfi_NTF),]), by='row.names', all=TRUE)
# row.names(patiAfi_NTF_and_cliniVars) <- patiAfi_NTF_and_cliniVars$Row.names
# patiAfi_NTF_and_cliniVars <- patiAfi_NTF_and_cliniVars[,-1]
# patiAfi_NTF_and_cliniVars$age <- as.numeric(as.character(patiAfi_NTF_and_cliniVars$patient.age_at_initial_pathologic_diagnosis))
# patiAfi_NTF_and_cliniVars$age.f <- cut(patiAfi_NTF_and_cliniVars$age, c(28,50,70,89), ordered_result = TRUE)
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
# table(patiAfi_NTF_and_cliniVars$Normalized.Stage, useNA = "always")
# anyNA(patiAfi_NTF_and_cliniVars$Normalized.Stage)
# # remove patients without AJCC.Stage
# patiAfi_NTF_and_cliniVars<-patiAfi_NTF_and_cliniVars[-which(is.na(patiAfi_NTF_and_cliniVars$Normalized.Stage) ),]
#
#
# ## check for empty values
# variabl <- "patient.histological_type"
# colmnAsChar <- as.character(patiAfi_NTF_and_cliniVars[,variabl])
# anyNA(colmnAsChar)
# table(colmnAsChar, useNA = "always")
# sum(table(colmnAsChar, useNA = "always"))
# # check NAs
# # View(patiAfi_NTF_and_cliniVars[which(is.na(colmnAsChar)),])
# ## remove patients that belong to factors with only one patient
# disCompo <- table(colmnAsChar)
# patiAfi_NTF_and_cliniVars[,"patient"] <- row.names(patiAfi_NTF_and_cliniVars)
# compo2del <- names(disCompo[which(disCompo==1)])
# if( length(compo2del) != 0){
#   patients2keep <- patiAfi_NTF_and_cliniVars [-which(colmnAsChar %in% compo2del), "patient"]
#   patients2keep <- as.character(patients2keep)
#   patiAfi_NTF_and_cliniVars <- patiAfi_NTF_and_cliniVars[patients2keep,]
# }
# patiAfi_NTF_and_cliniVars[,variabl] <- colmnAsChar
# table(patiAfi_NTF_and_cliniVars[,variabl], useNA = "always")
#
# ## check for empty values
# variabl <- "patient.breast_carcinoma_progesterone_receptor_status"
# colmnAsChar <- as.character(patiAfi_NTF_and_cliniVars[,variabl])
# anyNA(colmnAsChar)
# table(colmnAsChar, useNA = "always")
#
#
# save.image("temp/6-4v2.RData")

# TEST 1 = JUST COMPONENTs ----
totalNumOfSamples <- dim(patiAfi_NTF_and_cliniVars)[[1]]  #340 samples
numSplits = 10
smp_size <- floor(0.2 * totalNumOfSamples ) #testing sample size
set.seed(123)
obsInTest <- perrySplits(totalNumOfSamples, splitControl(m = smp_size, R = numSplits))
AUC_CD_K <- rep(0,numSplits)
Cstat_K <- rep(0,numSplits)
LogrankP_K <- rep(0,numSplits)
#times	The vector of time points at which AUC is evaluated.
utimes <- unique( patiAfi_NTF_and_cliniVars[patiAfi_NTF_and_cliniVars$Overall.Survival.Status == T,"Overall.Survival..Months."] )
utimes <- utimes[ order(utimes) ]
i=1
while (i <= numSplits) {
  test_ind <- obsInTest$subsets[,i]
  train <- patiAfi_NTF_and_cliniVars[-test_ind, ]
  test <- patiAfi_NTF_and_cliniVars[test_ind, ]
  #model
  coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                       event = Overall.Survival.Status)~
                    factor(component),
                  model=T, data=train)
  lp <- predict(coxFit)
  lpnew <- predict(coxFit, newdata=test)
  if(anyNA(lpnew)) {
    obsInTest <- perrySplits(totalNumOfSamples, splitControl(m = smp_size, R = numSplits))
    i=i-1
    print("*******lpnew had NA values**********")
    next
  }
  Surv.rsp <- Surv(time = train$Overall.Survival..Months.,
                   event = train$Overall.Survival.Status)
  Surv.rsp.new <- Surv(time = test$Overall.Survival..Months.,
                       event = test$Overall.Survival.Status)
  Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
  try(if(Cstat < 0.01) stop("**** Cstat equal zero *****"))
  AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
  AUC_CD_K[i] <- AUC_CD$iauc
  Cstat_K[i] <- Cstat
  i=i+1
}
mean(AUC_CD_K)
0.5426689
#0.5389344 slicely different to the AUC obtained with all patients


# ADD one by one the 9 variables -----

foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component) )
#sortedAUCs[1]



