#' version 7 (May/31/2018)
#' principal verotype as covariate

rm(list = ls())
load(file="temp/5-2v6-data.RData")
loadls("plyr survival Rcpp survAUC perry",F)

LoadMyData <- function(){
  load("data/survivalClinical-5-1_30.RData")
  rm(list=ls()[-which(ls() %in% c("survivalClinical",
                                  "patients","kMax", "numberOfPatiens"))])
  load("data/myLib.RData" ) #load my libraries
  loadls("plyr survival Rcpp survAUC perry",F)
  save.image(file="temp/5-2v6-data.RData")
}


#k=15

NTF_AUC <- function()
{
  kMax=30
  numSplits = 10
  AUC_CD_all <- rep(list(),kMax)
  Cstat_allK <- rep(list(),kMax)
  LogrankP_allK_train <- rep(list(),kMax)
  LogrankP_allK <- rep(0,kMax)

  smp_size <- floor(0.2 * numberOfPatiens) #testing sample size

  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == 1,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]



  for (k in c(2:kMax) ){
    load( file = paste("data/patientAfiliation/patientAfiliation_",k,".RData", sep='') )
    row.names(patientAfiliation) <- patientAfiliation$patient
    patiF <- merge(as.data.frame(patientAfiliation), as.data.frame(survivalClinical), by='row.names', all=TRUE)
    patiF <- patiF[,-1]
    patiF$component <- factor(patiF$component)
    #levels(patiF$component) <- c(levels(patiF$component) , c(1:k))

    AUC_CD_K <- rep(0,numSplits)
    Cstat_K <- rep(0,numSplits)
    LogrankP_K <- rep(0,numSplits)

    set.seed(k*1234)  # set seed for reproducibility
    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))

    for (i in seq(1,numSplits)) {
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
      #Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K[i] <- AUC_CD$iauc
      #Cstat_K[i] <- Cstat
      LogrankP_K[i] <- logp_val
    }

    coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                         event = Overall.Survival.Status)~component,
                    model=T, x=T, y=T, data=patiF)
    logp_val <- summary(coxFit)
    logp_val <- logp_val$sctest[3]
    LogrankP_allK[k] <- logp_val

    save.image(paste("temp/5-2v5-AUC_NTF_k",k,".RData",sep=''))
    AUC_CD_all[k] = list(AUC_CD_K)
    #Cstat_allK[k] <- list(Cstat_K)
    LogrankP_allK_train[k] <- list(LogrankP_K)

  }
  NTF_AUC <- AUC_CD_all
  #NTF_Cstat <- Cstat_allK
  NTF_logRank <- LogrankP_allK
  save(NTF_AUC, file = "output4paper/NTF_AUC.RData")
  #save(NTF_Cstat, file = "output4paper/NTF_Cstat.RData")
  save(NTF_logRank, file = "output4paper/NTF_logRank.RData")
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-2v6.RData")
