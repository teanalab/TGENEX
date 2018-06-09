
getAUCaddVar2component <- function(coxFunction, one_var_name){
  if( anyNA(patiAfi_NTF_and_cliniVars[,variabl]) )
    stop("******* NAs on variable ****")
  if(min(table(patiAfi_NTF_and_cliniVars[,variabl], useNA = "ifany")) == 0 )
    stop("******* NAs on variable ****")

  totalNumOfSamples <- dim(patiAfi_NTF_and_cliniVars)[[1]]
  numSplits = 10
  smp_size <- floor(0.2 * totalNumOfSamples ) #testing sample size
  set.seed(123)
  obsInTest <- perrySplits( totalNumOfSamples, splitControl(m = smp_size, R = numSplits) )
  AUC_CD_K <- rep(0, numSplits)
  Cstat_K <- rep(0, numSplits)
  LogrankP_K <- rep(0,numSplits)
  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( patiAfi_NTF_and_cliniVars[patiAfi_NTF_and_cliniVars$Overall.Survival.Status == T,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]
  i=1
  lasti=1
  while (i <= numSplits) {
    test_ind <- obsInTest$subsets[,i]
    train <- patiAfi_NTF_and_cliniVars[-test_ind, ]
    test <- patiAfi_NTF_and_cliniVars[test_ind, ]
    while (min(table(train$component))==0 || min(table(test$component))==0 ||
           min(table(train[,one_var_name]))==0 || min(table(test[,one_var_name]))==0 ||
           min(table(train$Overall.Survival.Status))==0 || min(table(test$Overall.Survival.Status))==0 ||
           any( levels(train[,one_var_name]) != levels(test[,one_var_name]) ) )   {
      obsInTest <- perrySplits(totalNumOfSamples, splitControl(m = smp_size, R = numSplits))
      test_ind <- obsInTest$subsets[,i]
      train <- patiAfi_NTF_and_cliniVars[-test_ind, ]
      test <- patiAfi_NTF_and_cliniVars[test_ind, ]
      print("*** non zero components && levels ****")
    }

    #model
    coxFit <- coxph(foo, model=T, data=train)
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
}


#' GetAUC
#' only one variable at the time
getAUC <- function(coxFunction, one_var_name){
  if( anyNA(patiAfi_NTF_and_cliniVars[,variabl]) )
    stop("******* NAs on variable ****")
  if(min(table(patiAfi_NTF_and_cliniVars[,variabl], useNA = "ifany")) == 0 )
    stop("******* NAs on variable ****")

  totalNumOfSamples <- dim(patiAfi_NTF_and_cliniVars)[[1]]
  numSplits = 10
  smp_size <- floor(0.2 * totalNumOfSamples ) #testing sample size
  set.seed(123)
  obsInTest <- perrySplits( totalNumOfSamples, splitControl(m = smp_size, R = numSplits) )
  AUC_CD_K <- rep(0, numSplits)
  Cstat_K <- rep(0, numSplits)
  LogrankP_K <- rep(0,numSplits)
  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( patiAfi_NTF_and_cliniVars[patiAfi_NTF_and_cliniVars$Overall.Survival.Status == T,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]
  i=1
  lasti=1
  while (i <= numSplits) {
    test_ind <- obsInTest$subsets[,i]
    train <- patiAfi_NTF_and_cliniVars[-test_ind, ]
    test <- patiAfi_NTF_and_cliniVars[test_ind, ]
    while (min(table(train[,one_var_name]))==0 || min(table(test[,one_var_name]))==0 ||
           min(table(train$Overall.Survival.Status))==0 || min(table(test$Overall.Survival.Status))==0 ||
           any( levels(train[,one_var_name]) != levels(test[,one_var_name]) ) )   {
      obsInTest <- perrySplits(totalNumOfSamples, splitControl(m = smp_size, R = numSplits))
      test_ind <- obsInTest$subsets[,i]
      train <- patiAfi_NTF_and_cliniVars[-test_ind, ]
      test <- patiAfi_NTF_and_cliniVars[test_ind, ]
      print("*** non zero components && levels ****")
    }

    #model
    coxFit <- coxph(foo, model=T, data=train)
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
}
