
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

