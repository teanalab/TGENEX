options(warn=1) #normal
options(warn=2) #as error
coxFit <-  try(coxph(Surv(time = Overall.Survival..Months.,
                                                event = Overall.Survival.Status)~
                                             factor(component) +
                                             # age.f +
                                             #factor(Normalized.Stage) +
                                             #factor(Node.Coded) +
                                             #factor(Metastasis.Coded)
                                             #factor(Converted.Tumor)#+
                                             # factor(patient.histological_type)+
                                             # factor(HER2.Final.Status)+
                                             # factor(patient.breast_carcinoma_estrogen_receptor_status)+
                                             factor(patient.breast_carcinoma_progesterone_receptor_status)
                                           ,x=T, y=T, model=T, method = "breslow",
                                           data=patiAfi_NTF_and_cliniVars)
         )
if( inherits(coxFit,'try-error') )
{
  stop(" something went wrong  ")
}

summary(coxFit)
