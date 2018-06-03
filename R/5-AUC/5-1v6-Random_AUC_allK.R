#' version 4 (May/30/2018)
#' AUC with random clustering

rm(list = ls())
load(file="data/survivalClinical-5-1_30.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)

LoadMyData <- function(){
  #load my libraries
  loadlib<-c("Packages.R")
  loadlib<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", loadlib,sep='')
  sapply(loadlib, function(u) {source(u)})
  #required for the cox proportional hazard model
  loadls("plyr survival Rcpp missForest survAUC",FALSE)

  load("data/survClinical.RData")
  load("data/patients.RData")

  survivalClinical <- survClinical
  survivalClinical[,"Patient.ID"] <- survivalClinical$patient.bcr_patient_barcode
  survivalClinical[,"Overall.Survival.Status"] <- (survivalClinical$patient.vital_status=="dead")
  survivalClinical[,"Overall.Survival..Months."] <- survivalClinical$OS_MONTHS

  survivalClinical <-  survivalClinical[which(survivalClinical$Patient.ID %in% patients),]
  numberOfPatiens <- dim(survivalClinical)[[1]]
  kMax <- 30
  row.names(survivalClinical) <- survivalClinical$Patient.ID
  survivalClinical <- survivalClinical[,c("Overall.Survival.Status","Overall.Survival..Months.")]
  save.image(file="data/survivalClinical-5-1_30.RData")
}

#line by line
#randomAUCs <- function()
#{
  kMax = 9
  AUC_CD_allK <- rep(list(),kMax)
  Cstat_allK <- rep(list(),kMax)
  LogrankP_allK_train <- rep(list(),kMax)
  LogrankP_allK <- rep(0,kMax)

  numSplits = 10
  percenTesting = 0.2

  set.seed(100)

  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == 1,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]

  for (k in c(2:kMax) ){
    patVarName <- paste("randomF_R_",k, sep='')
    load(file= paste("data/random/",patVarName,".RData",sep='')  )
    patientAfiliation <- get(patVarName)
    patiF <- cbind(patientAfiliation,survivalClinical)
    patiF$component <- factor(patiF$component)

    AUC_CD_K <- rep(0,numSplits)
    Cstat_K <- rep(0,numSplits)
    LogrankP_K <- rep(0,numSplits)

    set.seed(k*1234)  # set seed for reproducibility
    smp_size <- floor(0.2 * numberOfPatiens)  #testing sample

    ## data folds for K-fold cross-validation
    # perrySplits(20, foldControl(K = 5))
    # perrySplits(20, foldControl(K = 5, R = 10))
    ## random data splits
    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))
    ## bootstrap samples
    # obsInTest <- perrySplits(numberOfPatiens, bootControl())

    i=1
    while (i <= numSplits) {
      set.seed(i)
      test_ind <- obsInTest$subsets[,i]
      train <- patiF[-test_ind, ]
      test <- patiF[test_ind, ]


      while (min(table(train$component))==0 || min(table(test$component))==0 ||
             min(table(train$Overall.Survival.Status))==0 || min(table(test$Overall.Survival.Status))==0)
      {
        obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))
        test_ind <- obsInTest$subsets[,i]
        train <- patiF[-test_ind, ]
        test <- patiF[test_ind, ]
        print("resample")
      }


      #cox proportional hazard model
      coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~component,
                      model=TRUE,x=T,y=T, data=train)
      logp_val <- summary(coxFit)
      logp_val <- logp_val$sctest[3]

      ## AUC
      lp <- predict(coxFit)
      lpnew <- predict(coxFit, newdata=test)
      Surv.rsp <- Surv(time = train$Overall.Survival..Months.,
                       event = train$Overall.Survival.Status)
      Surv.rsp.new <- Surv(time = test$Overall.Survival..Months.,
                           event = test$Overall.Survival.Status)
      Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K[i] <- AUC_CD$iauc
      Cstat_K[i] <- Cstat
      LogrankP_K[i] <- logp_val

    }

    ## overall patients
    coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                         event = Overall.Survival.Status)~component,
                    model=T, x=T, y=T, data=patiF)
    logp_val <- summary(coxFit)
    logp_val <- logp_val$sctest[3]
    LogrankP_allK[k] <- logp_val
    i=i+1

    save.image(paste("temp/5-1v3-AUC_Random_k",k,".RData",sep=''))
    AUC_CD_allK[k] = list(AUC_CD_K)
    Cstat_allK[k] <- list(Cstat_K)
    LogrankP_allK_train[k] <- list(LogrankP_K)

  }
  randomAUC <- AUC_CD_allK
  randomCstat <- Cstat_allK
  random_logRank_train <- LogrankP_allK_train
  random_logRank <- LogrankP_allK
  save(randomAUC, file = "output4paper/randomAUC.RData")
  save(randomCstat, file = "output4paper/randomCstat.RData")
  save(random_logRank, file = "output4paper/random_logRank.RData")
  save(random_logRank_train, file = "output4paper/random_logRank_train.RData")
#}

sess <- sessionInfo() #save session on variable
save.image("temp/5-1v5.RData")
