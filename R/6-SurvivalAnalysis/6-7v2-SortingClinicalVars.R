
# ## check for empty values
# sum(table(patiAfi_NTF_and_cliniVars[,variabl], useNA = "always"))
# # check NAs
# # View(patiAfi_NTF_and_cliniVars[which(is.na(patiAfi_NTF_and_cliniVars[,variabl])),])
# ## remove patients that belong to factors with only one patient
# disCompo <- table(patiAfi_NTF_and_cliniVars[,variabl])
# patiAfi_NTF_and_cliniVars[,"patient"] <- row.names(patiAfi_NTF_and_cliniVars)
# compo2del <- names(disCompo[which(disCompo==1)])
# if( length(compo2del) != 0){
#   patients2keep <- patiAfi_NTF_and_cliniVars [-which(patiAfi_NTF_and_cliniVars[,variabl] %in% compo2del), "patient"]
#   patients2keep <- as.character(patients2keep)
#   patiAfi_NTF_and_cliniVars <- patiAfi_NTF_and_cliniVars[patients2keep,]
# }

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

AUCsingleCovariates <- rep(0.0,10)

#ten variables
VarNames <- c("factor(component)",
              "age.f",
              "factor(Normalized.Stage)",
              "factor(Node.Coded)",
              "factor(Metastasis.Coded)",
              "factor(Converted.Tumor)",
              "factor(patient.histological_type)",
              "factor(HER2.Final.Status)",
              "factor(patient.breast_carcinoma_estrogen_receptor_status)",
              "factor(patient.breast_carcinoma_progesterone_receptor_status)")

foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component) )
AUCsingleCovariates[1] <- getAUC(foo,"component")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                    age.f)
AUCsingleCovariates[2] <- getAUC(foo,"age.f")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(Normalized.Stage))
AUCsingleCovariates[3] <- getAUC(foo,"Normalized.Stage")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(Node.Coded)  )
AUCsingleCovariates[4] <- getAUC(foo, "Node.Coded")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(Metastasis.Coded))
AUCsingleCovariates[5] <- getAUC(foo, "Metastasis.Coded")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(Converted.Tumor))
AUCsingleCovariates[6] <- getAUC(foo,"Converted.Tumor")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(patient.histological_type)
               )
AUCsingleCovariates[7] <- getAUC(foo,"patient.histological_type")



foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(HER2.Final.Status)
)
AUCsingleCovariates[8] <- getAUC(foo, "HER2.Final.Status")



foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(patient.breast_carcinoma_estrogen_receptor_status)
)
AUCsingleCovariates[9] <- getAUC(foo, "patient.breast_carcinoma_estrogen_receptor_status")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(patient.breast_carcinoma_progesterone_receptor_status)
)
AUCsingleCovariates[10] <- getAUC(foo, "patient.breast_carcinoma_progesterone_receptor_status")

names(AUCsingleCovariates) <- VarNames

sortedAUCs <- sort(AUCsingleCovariates, decreasing = T)

factor(component)
0.51819707
factor(patient.breast_carcinoma_progesterone_receptor_status)
0.37719690
factor(Normalized.Stage)
0.37284256
age.f
0.36107346
factor(patient.breast_carcinoma_estrogen_receptor_status)
0.32727184
factor(Node.Coded)
0.30950702
factor(Converted.Tumor)
0.20862196
factor(HER2.Final.Status)
0.17987748
factor(patient.histological_type)
0.12567533
factor(Metastasis.Coded)
0.08980807

