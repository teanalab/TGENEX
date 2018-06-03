#' version 4 (May/30/2018)
#' AUC with random clustering

rm(list = ls())
load(file="data/survivalClinical-5-1_30.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)

#line by line
randomAUCs <- function()
{
  prefixPatAfilFile <- "data/random/"
  ouputFileEachK <- "temp/5-1v3-AUC_Random_k"
  folderOutput <- "output4paper/random_"
  kMax=12
  numSplits = 10
  AUC_CD_all <- rep(list(),kMax)
  Cstat_allK <- rep(list(),kMax)
  LogrankP_allK_train <- rep(list(),kMax)
  LogrankP_allK <- rep(0,kMax)

  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == T,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]

  for (k in c(2:kMax) ){
    AUC_CD_K <- rep(0,numSplits)
    Cstat_K <- rep(0,numSplits)
    LogrankP_K <- rep(0,numSplits)
    set.seed(k*1234)  # set seed for reproducibility
    patVarName <- paste("randomF_R_",k, sep='')
    load(file= paste(prefixPatAfilFile,patVarName,".RData",sep='')  )
    patientAfiliation <- get(patVarName)
    row.names(patientAfiliation) <- patientAfiliation$patient
    #remove patients that belong to components with only one patient
    disCompo <- table(patientAfiliation$component)
    compo2del <- names(disCompo[which(disCompo==1)])
    if( length(compo2del) != 0){
      patients2keep <- patientAfiliation [-which(patientAfiliation$component %in% compo2del), "patient"]
      patients2keep <- as.character(patients2keep)
      patientAfiliation <- patientAfiliation[patients2keep,]
      survivalClinical2 <- survivalClinical[patients2keep,]
      #TODO
      patiF <- cbind(patientAfiliation,survivalClinical)
      patiF$component <- factor(patiF$component)
      #levels(patiF$component) <- c(levels(patiF$component) , c(1:k))
      numOfRemainingP <- length(patients2keep)
      smp_size <- floor(0.2 * numOfRemainingP ) #testing sample size
      obsInTest <- perrySplits( numOfRemainingP, splitControl(m = smp_size, R = numSplits))
    } else {
      #TODO
      patiF <- cbind(patientAfiliation,survivalClinical)
      patiF$component <- factor(patiF$component)
      #levels(patiF$component) <- c(levels(patiF$component) , c(1:k))
      smp_size <- floor(0.2 * numberOfPatiens) #testing sample size
      obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))
      numOfRemainingP <- numberOfPatiens
    }

    i=1
    while (i <= numSplits) {
      test_ind <- obsInTest$subsets[,i]
      train <- patiF[-test_ind, ]
      test <- patiF[test_ind, ]

      #prediction on test data gives NA if there is atleast 1 component without samples
      #Error in Cox if there are no events on training or testing

      # table(train$component)
      # table(test$component)
      # table(train$Overall.Survival.Status)
      # table(test$Overall.Survival.Status)


      while (min(table(train$component))==0 || min(table(test$component))==0 ||
             min(table(train$Overall.Survival.Status))==0 || min(table(test$Overall.Survival.Status))==0)
      {
        obsInTest <- perrySplits(numOfRemainingP, splitControl(m = smp_size, R = numSplits))
        test_ind <- obsInTest$subsets[,i]
        train <- patiF[-test_ind, ]
        test <- patiF[test_ind, ]
        print("resample")
      }
      #cox proportional hazard model
      coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~component,
                      model=T, data=train)
      logp_val <- summary(coxFit)
      logp_val <- logp_val$sctest[3]

      #AUC
      lp <- predict(coxFit, type= "expected")
      lpnew <- predict(coxFit, newdata=test, type= "expected")
      if(anyNA(lpnew)) {
        obsInTest <- perrySplits(numOfRemainingP, splitControl(m = smp_size, R = numSplits))
        i=i-1
        print("*******lpnew**********")
        next
      }
      Surv.rsp <- Surv(time = train$Overall.Survival..Months.,
                       event = train$Overall.Survival.Status)
      Surv.rsp.new <- Surv(time = test$Overall.Survival..Months.,
                           event = test$Overall.Survival.Status)
      Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
      try(if(Cstat < 0.01) stop("****Cstat*****"))
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K[i] <- AUC_CD$iauc
      Cstat_K[i] <- Cstat
      LogrankP_K[i] <- logp_val
      i=i+1
    }

    ## overall patients
    coxFit <- coxph(Surv(time = Overall.Survival..Months., event = Overall.Survival.Status)~component,
                    model=T, x=T, y=T, data=patiF)
    logp_val <- summary(coxFit)
    logp_val <- logp_val$sctest[3]
    LogrankP_allK[k] <- logp_val

    AUC_CD_all[k] = list(AUC_CD_K)
    Cstat_allK[k] <- list(Cstat_K)
    LogrankP_allK_train[k] <- list(LogrankP_K)
    save.image(paste(ouputFileEachK,k,".RData",sep=''))
  }
  random_AUC <- AUC_CD_all
  random_Cstat <- Cstat_allK
  random_logRank <- LogrankP_allK
  random_logRank_train <- LogrankP_allK_train
  save(random_AUC, file = paste(folderOutput, "AUC.RData", sep='') )
  save(random_Cstat, file =  paste(folderOutput, "Cstat.RData", sep='') )
  save(random_logRank, file = paste(folderOutput, "logRank.RData", sep='') )
  save(random_logRank_train, file = paste(folderOutput, "logRank_train.RData", sep='') )
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-1v5.RData")
