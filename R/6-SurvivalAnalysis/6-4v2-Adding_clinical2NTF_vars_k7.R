#version 2 (Jun/8/2018)
# Adding AUC

rm(list = ls())


load("temp/6-4v1.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)

# top 5 components alone p-value=0.03779645
# AUC = 0.5389344


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

coxFit_tensorClini <- coxph(Surv(time = Overall.Survival..Months.,
                                 event = Overall.Survival.Status)~
                              factor(component)+
                              # age.f +
                              factor(Normalized.Stage) +
                              factor(Node.Coded) +
                              #factor(Metastasis.Coded)
                              factor(Converted.Tumor)#+
                              # factor(patient.histological_type)+
                              # factor(HER2.Final.Status)+
                              # factor(patient.breast_carcinoma_estrogen_receptor_status)+
                              # factor(patient.breast_carcinoma_progesterone_receptor_status)
                              ,x=T, y=T, model=T, method = "breslow",
                              data=patiAfi_NTF_and_cliniVars)
summary(coxFit_tensorClini)


# AUC -----
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

  #training
  coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                       event = Overall.Survival.Status)~component,
                  model=T, data=train)
  logp_val <- summary(coxFit)
  logp_val <- logp_val$sctest[3]


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
  LogrankP_K[i] <- logp_val
  i=i+1
}

mean(AUC_CD_K)

#AUC when adding clinical variables = 0.3952659
