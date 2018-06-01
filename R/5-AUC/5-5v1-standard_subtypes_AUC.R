#' version 1 (May/31/2018)
#' Subtypes:
#' Luminal A (ER+ and/or PR+, HER2-)
# Luminal B (ER+ and/or PR+, HER2+)
# HER2+ overexpressing (ER-, PR- and HER2+)
# Basal-like/Triple negative (ER-, PR-, HER2- CK 5/6+ and/or HER1+)


rm(list = ls())

#LoadMyData <- function(){ #-------
  load("data/survivalClinical-5-1_30.RData")
  load("data/myLib.RData" ) #load my libraries
  loadls("plyr survival Rcpp survAUC perry",F)
#}

#standard_subtypes_AUC <- function() { #-------
  numSplits = 10
  smp_size <- floor(0.2 * numberOfPatiens) #testing sample size

  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == 1,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]

  row.names(survClinical) <- survClinical$patient.bcr_patient_barcode
  survClinical <- survClinical[patients,]
  table(survClinical[patients,"PAM50.Subtype"] )

  patiF <- survClinical[,c("patient.bcr_patient_barcode","PAM50.Subtype")]
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiF <- patiF[,-1]
  patiF$component <- factor(patiF$PAM50.Subtype)
  #levels(patiF$component) <- c(levels(patiF$component) , c(1:k))

  AUC_CD_a <- rep(0,numSplits)
  Cstat_K <- rep(0,numSplits)
  LogrankP_K <- rep(0,numSplits)

  set.seed(1234)  # set seed for reproducibility
  obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))

  #cox proportional hazard model
  coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                       event = Overall.Survival.Status)~component,
                  model=T, x=T, y=T, data=patiF)
  logp_val <- summary(coxFit)
  logp_val <- logp_val$sctest[3]


  for (i in seq(numSplits) ) {
    test_ind <- obsInTest$subsets[,i]
    train <- patiF[-test_ind, ]
    test <- patiF[test_ind, ]

    #cox proportional hazard model
    coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                         event = Overall.Survival.Status)~component,
                    model=T, x=T, y=T, data=train)
    logp_val <- summary(coxFit)
    logp_val <- logp_val$sctest[3]

    #AUC
    lp <- predict(coxFit, type= "expected")
    lpnew <- predict(coxFit, newdata=test, type= "expected")
    Surv.rsp <- Surv(time = train$Overall.Survival..Months.,
                     event = train$Overall.Survival.Status)
    Surv.rsp.new <- Surv(time = test$Overall.Survival..Months.,
                         event = test$Overall.Survival.Status)
    Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
    AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
    #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
    AUC_CD_a[i] <- AUC_CD$iauc
    Cstat_K[i] <- Cstat
    LogrankP_K[i] <- logp_val
  }
  PAM50_AUC <- AUC_CD_a
  PAM50_Cstat <- Cstat_K
  PAM50_logRank <- LogrankP_K
  save(PAM50_AUC, PAM50_Cstat, PAM50_logRank, file = "output4paper/PAM50_AUC.RData")
  mean(PAM50_AUC)
  sd(PAM50_AUC)

  mean(PAM50_Cstat)
  sd(PAM50_Cstat)

  mean(PAM50_logRank)
  sd(PAM50_logRank)

#}

sess <- sessionInfo() #save session on variable
save.image("temp/5-2v6.RData")
