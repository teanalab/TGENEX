#reading survival

load("dataExample/survival/survivalClinical.RData")

rm(list = ls()[-which(ls() %in% c("survClinical")  )])

save(survClinical, file="dataExample/survival/survClinical.RData")


row.names(clinicalLong) <- as.character(clinicalLong$patient.bcr_patient_barcode)
clinicalLong <- clinicalLong[as.character(patients),c( "patient.vital_status", "OS_MONTHS")]

clinicalLong$patient.vital_status
clinicalLong$patient.year_of_form_completion
clinicalLong$patient.year_of_initial_pathologic_diagnosis

clinicalLong$patient.month_of_form_completion
clinicalLong$patient.days_to_last_followup
clinicalLong$patient.days_to_death
View(clinicalLong[,c("patient.days_to_death","patient.days_to_last_followup")])

surv2 <- clinicalLong[,c("patient.days_to_death","patient.days_to_last_followup")]
surv2 <- as.matrix(surv2)
as.numeric(surv2[1,2])/30.431
View(surv2)

#NA to 0
surv2[which(is.na(surv2[,1])), 1] = 0
surv2[which(is.na(surv2[,2])), 2] = 0

sumSurv = as.numeric(surv2[,1]) + as.numeric(surv2[,2])

which(sumSurv ==0 )

View(sumSurv/30.431)

surv3 <- clinicalLong[,c("patient.days_to_death","patient.days_to_last_followup")]
surv3[76,]
row.names(surv3[which(sumSurv ==0 ), ])

pp1 <- as.character(survClinical$patient.bcr_patient_barcode)
pp2 <- as.character(clinicalLong$patient.bcr_patient_barcode)
pp2 <- as.character(clinicalLong$patient.bcr_patient_barcode)
anyNA(survClinical)

length(intersect(pp1,pp2))
missingP <- setdiff(pp2,pp1)
surv3[missingP,]

View(survClinical[row.names(surv3[which(sumSurv ==0 ), ]),])


survivalClinical[missingP[1], ] <- c("alive", round(as.numeric(surv2[missingP[1],"patient.days_to_last_followup"])/30.4,2)
)
survivalClinical[missingP[2], ] <- c("alive", round(as.numeric(surv2[missingP[2],"patient.days_to_last_followup"])/30.4,2)
)


row.names(survClinical) <- as.character(survClinical$patient.bcr_patient_barcode)

survivalClinical <- survClinical[,c( "patient.vital_status", "OS_MONTHS")]
survivalClinical$patient.vital_status = survivalClinical$patient.vital_status != "alive"
names(survivalClinical) <- c( "Overall.Survival.Status", "Overall.Survival..Months.")
View(survivalClinical)


survivalClinical <- survivalClinical[patients,]

survivalClinical$Overall.Survival..Months. <- as.numeric(survivalClinical$Overall.Survival..Months.)

save(survivalClinical, file="dataExample/Filtered/survivalClinical.RData")

######