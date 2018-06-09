#version 2 (Jun/8/2018)
# Adding AUC

rm(list = ls())


load("temp/6-4v1.RData")
load("data/myLib.RData")
loadls("plyr survival missForest survAUC prodlim survminer perry", F)

# top 5 components alone p-value=0.03779645

#Adding clinical variables
patiAfi_NTF <- patientAfiliation_NTF[which(as.numeric(patientAfiliation_NTF$component) %in% bc4),]
patiAfi_NTF_and_cliniVars <- merge( as.data.frame(patiAfi_NTF),as.data.frame(tensorClinical[row.names(patiAfi_NTF),]), by='row.names', all=TRUE)
patiAfi_NTF_and_cliniVars <- patiAfi_NTF_and_cliniVars[,-1]
patiAfi_NTF_and_cliniVars$age <- as.numeric(as.character(patiAfi_NTF_and_cliniVars$patient.age_at_initial_pathologic_diagnosis))
patiAfi_NTF_and_cliniVars$age.f <- cut(patiAfi_NTF_and_cliniVars$age, c(28,50,70,85), ordered_result = TRUE)
patiAfi_NTF_and_cliniVars$Converted.Tumor <- with(patiAfi_NTF_and_cliniVars, ifelse(patiAfi_NTF_and_cliniVars$Tumor == "T1", "T1", "T_other"))

# Stage IIIB has zero TRUE status
# solution: group all stage III together
patiAfi_NTF_and_cliniVars[grep("I",patiAfi_NTF_and_cliniVars$AJCC.Stage),
                          "Normalized.Stage"] = "I"
patiAfi_NTF_and_cliniVars[grep("II",patiAfi_NTF_and_cliniVars$AJCC.Stage),
                          "Normalized.Stage"] = "II"
patiAfi_NTF_and_cliniVars[grep("III",patiAfi_NTF_and_cliniVars$AJCC.Stage),
                          "Normalized.Stage"] = "III"
patiAfi_NTF_and_cliniVars[grep("IV",patiAfi_NTF_and_cliniVars$AJCC.Stage),
                          "Normalized.Stage"] = "IV"

# variable stats
table(patiAfi_NTF_and_cliniVars$Normalized.Stage)
with(patiAfi_NTF_and_cliniVars, table(Overall.Survival.Status, patient.histological_type))


###
coxFit_tensorClini <- coxph(Surv(time = Overall.Survival..Months.,
                                 event = Overall.Survival.Status)~
                              factor(component) +
                              # age.f +
                              #factor(Normalized.Stage) +
                              #factor(Node.Coded) +
                              #factor(Metastasis.Coded)
                              #factor(Converted.Tumor)#+
                              # factor(patient.histological_type)+
                              # factor(HER2.Final.Status)+
                              # factor(patient.breast_carcinoma_estrogen_receptor_status)+
                               factor(patient.breast_carcinoma_progesterone_receptor_status)
                              ,x=T, y=T, model=T, method = "breslow",
                              data=patiAfi_NTF_and_cliniVars)
summary(coxFit_tensorClini)
###

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
#iAUC = 0.4972065
#######



# TEST 2 = JUST COMPONENTs + factor(patient.breast_carcinoma_progesterone_receptor_status)----
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
  while (min(table(train$component))==0 || min(table(test$component))==0 ||
         min(table(train$Overall.Survival.Status))==0 || min(table(test$Overall.Survival.Status))==0 ||
         any(levels(train$patient.breast_carcinoma_progesterone_receptor_status) != levels(test$patient.breast_carcinoma_progesterone_receptor_status))
         )
  {
    obsInTest <- perrySplits(totalNumOfSamples, splitControl(m = smp_size, R = numSplits))
    test_ind <- obsInTest$subsets[,i]
    train <- patiF[-test_ind, ]
    test <- patiF[test_ind, ]
    print("*** non zero components ****")
  }



  #model
  coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                            event = Overall.Survival.Status)~
                    factor(component)+
                    factor(patient.breast_carcinoma_progesterone_receptor_status),
                  model=T, data=train)
  scF <- summary(coxFit)
  if( any( scF$conf.int[,4] == Inf) ){ #resample
    obsInTest <- perrySplits(totalNumOfSamples, splitControl(m = smp_size, R = numSplits))
    i=lasti
    print("*******CoxFit = beta may be infinite. **********")
    next
  }

  lp <- predict(coxFit)
  lpnew <- predict(coxFit, newdata=test)
  if(anyNA(lpnew)) { #resample
    obsInTest <- perrySplits(totalNumOfSamples, splitControl(m = smp_size, R = numSplits))
    i=lasti
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
  lasti = i
}
mean(AUC_CD_K)
#baseline iAUC = 0.4972065
#new iAUC = 0.5586578
