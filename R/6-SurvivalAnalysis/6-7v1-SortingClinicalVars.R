
SignificanceOfEachVariable <- rep(list(),10)

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



coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                                 event = Overall.Survival.Status)~
                          factor(component),
                            data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[1] <- list(summary(coxFitOneByOne))



coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                                 event = Overall.Survival.Status)~age.f
                            ,x=T, y=T, model=T, method = "breslow",
                            data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[2] <- list(summary(coxFitOneByOne))




coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(Normalized.Stage)
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[3] <- list(summary(coxFitOneByOne))


coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(Node.Coded)
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[4] <- list(summary(coxFitOneByOne))


coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(Metastasis.Coded)
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[5] <- list(summary(coxFitOneByOne))


coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(Converted.Tumor)
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[6] <- list(summary(coxFitOneByOne))


coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(patient.histological_type)+
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[7] <- list(summary(coxFitOneByOne))


coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(HER2.Final.Status)
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[8] <- list(summary(coxFitOneByOne))


coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(patient.breast_carcinoma_estrogen_receptor_status)
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[9] <- list(summary(coxFitOneByOne))


coxFitOneByOne <- coxph(Surv(time = Overall.Survival..Months.,
                             event = Overall.Survival.Status)~
                          factor(patient.breast_carcinoma_progesterone_receptor_status)
                        ,x=T, y=T, model=T, method = "breslow",
                        data=patiAfi_NTF_and_cliniVars)
SignificanceOfEachVariable[10] <- list(summary(coxFitOneByOne))



x <- SignificanceOfEachVariable[1]
y <- x[[1]]
z <- y$coefficients
All_Pr_z <- z[,5]

for (i in (2:10))
{
  x <- SignificanceOfEachVariable[i]
  y <- x[[1]]
  z <- y$coefficients
  Pr_z <- z[,5]
  if( is.null(names(Pr_z)) )
  {
    names(Pr_z) <- VarNames[i]
  }
  All_Pr_z <- c(All_Pr_z,Pr_z)
}


sort(All_Pr_z, decreasing = T)



factor(patient.breast_carcinoma_progesterone_receptor_status)negative
9.779495e-01
factor(component)6
8.472098e-01
factor(component)5
7.173838e-01
age.f.Q
6.974507e-01
factor(component)7
5.917183e-01
factor(patient.breast_carcinoma_progesterone_receptor_status)positive
5.701255e-01
age.f.L
3.879412e-01
factor(component)4
3.240272e-01
factor(Converted.Tumor)
3.218340e-01
factor(patient.histological_type)
3.218340e-01
factor(Normalized.Stage)II
2.700847e-01
factor(Normalized.Stage)III
1.851674e-01
factor(HER2.Final.Status)Negative
1.473431e-01
factor(HER2.Final.Status)Positive
1.123228e-01
factor(Node.Coded)
2.281950e-02
factor(patient.breast_carcinoma_estrogen_receptor_status)negative
2.809841e-03
factor(Normalized.Stage)IV
5.132184e-04
factor(Metastasis.Coded)
1.333741e-04
factor(patient.breast_carcinoma_estrogen_receptor_status)positive
5.767354e-05
