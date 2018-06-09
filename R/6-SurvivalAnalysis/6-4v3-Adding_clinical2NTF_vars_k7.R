#version 2 (Jun/9/2018)
# combining each clinical variable

rm(list = ls())


load("temp/6-4v2.RData")
load("data/myLib.RData")
source("R/6-SurvivalAnalysis/getAUCs.R")
loadls("plyr survival missForest survAUC prodlim survminer perry", F)


# baseline
# 0.5426689
# #0.5389344 slicely different to the AUC obtained with all patients


AUCadding1Covariate <- rep(0.0,9)

#ten variables
VarNames <- c("age.f",
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
                 factor(component)+
                 age.f)
AUCadding1Covariate[1] <- getAUCaddVar2component(foo,"age.f")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(Normalized.Stage))
AUCadding1Covariate[2] <- getAUCaddVar2component(foo,"Normalized.Stage")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(Node.Coded)  )
AUCadding1Covariate[3] <- getAUCaddVar2component(foo, "Node.Coded")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(Metastasis.Coded))
AUCadding1Covariate[4] <- getAUCaddVar2component(foo, "Metastasis.Coded")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(Converted.Tumor))
AUCadding1Covariate[5] <- getAUCaddVar2component(foo,"Converted.Tumor")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(patient.histological_type)
)
AUCadding1Covariate[6] <- getAUCaddVar2component(foo,"patient.histological_type")



foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(HER2.Final.Status)
)
AUCadding1Covariate[7] <- getAUCaddVar2component(foo, "HER2.Final.Status")



foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(patient.breast_carcinoma_estrogen_receptor_status)
)
AUCadding1Covariate[8] <- getAUCaddVar2component(foo, "patient.breast_carcinoma_estrogen_receptor_status")


foo <- formula(Surv(time = Overall.Survival..Months.,
                    event = Overall.Survival.Status)~
                 factor(component)+
                 factor(patient.breast_carcinoma_progesterone_receptor_status)
)
AUCadding1Covariate[9] <- getAUCaddVar2component(foo, "patient.breast_carcinoma_progesterone_receptor_status")

names(AUCadding1Covariate) <- VarNames

sortedAUCs <- sort(AUCadding1Covariate, decreasing = T)

factor(component)+age.f
0.6644742
factor(component)+factor(Normalized.Stage)
0.6509436
factor(component)+factor(patient.breast_carcinoma_estrogen_receptor_status)
0.6089686
factor(component)+factor(patient.breast_carcinoma_progesterone_receptor_status)
0.6039567
factor(component)+factor(Node.Coded)
0.6021760
factor(component)+factor(HER2.Final.Status)
0.5907964
factor(component)+factor(Metastasis.Coded)
0.5711561
factor(component)+factor(Converted.Tumor)
0.5674793
factor(component)+factor(patient.histological_type)
0.5542064
