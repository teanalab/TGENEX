rm(list = ls())
load(file="temp/boolMutation_OVA.RData")


patients



patientsOVA <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_output_sample_id4.csv", sep = "\t", quote = "'",
                          header = FALSE, stringsAsFactors = FALSE)

patientsOVA <- tolower(as.character(unlist(patientsOVA)))


ppp <- intersect(patientsOVA,patients)


clinical<-read.table(file="/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_clinical_patient.txt", sep = "\t",
                     header = FALSE, stringsAsFactors = FALSE,  fill = TRUE)
names(clinical) <- clinical[1,]
clinical <- clinical[-1,]

clinicalSample<-read.table(file="/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_clinical_sample.txt", sep = "\t",
                           header = FALSE, stringsAsFactors = FALSE,  fill = TRUE)
names(clinicalSample) <- clinicalSample[1,]
clinicalSample <- clinicalSample[-1,]
clinical <- merge(clinical,clinicalSample)

clinicalPatientsSample <- substring(tolower(clinicalSample$PATIENT_ID), 1, 12)

#phenotype
ppp <- intersect(patientsOVA,clinicalPatientsSample)

clinicalSample$shortID <- clinicalPatientsSample
clinicalSample2 <- clinicalSample[!duplicated(clinicalSample$shortID),]
row.names(clinicalSample2) <- clinicalSample2$shortID

clinicalSample3 <- clinicalSample2[ppp,-c(1,2,12)]


#survival
clinicalPatients <- substring(tolower(clinical$PATIENT_ID), 1, 12)
ppp <- intersect(patientsOVA,clinicalPatients)

clinical$shortID <- clinicalPatients
clinical2 <- clinical[!duplicated(clinical$shortID),]
row.names(clinical2) <- clinical2$shortID

clinical3 <- clinical2[ppp,-c(1,2,12)]


#firebrowse
clinicalFb<-read.table(file="/Users/diamac/10-Research/6Preprocessed_DataSets/firebrowser/ov_tcga_firebrowser/gdac.broadinstitute.org_OV.Merge_Clinical.Level_1.2016012800.0.0/OV.clin.merged.txt",
                       sep = "\t", header = FALSE, stringsAsFactors = FALSE,  fill = TRUE)
row.names(clinicalFb) <- clinicalFb[,1]
clinicalFb <- clinicalFb[,-1]

pppp <- intersect(as.character(clinicalFb[17,]), patientsOVA)

row.names(clinicalFb)[which(regexpr("status",row.names(clinicalFb)) != -1)]

clinicalFb["patient.vital_status",]

row.names(clinicalFb)[which(regexpr("dea",row.names(clinicalFb)) != -1)]

clinicalFb["patient.days_to_death", ]


row.names(clinicalFb)[which(regexpr("days",row.names(clinicalFb)) != -1)]

clinicalFb["patient.days_to_death", ]





####### medium --------



x = data.frame(clinicalFb[which(regexpr("days",row.names(clinicalFb)) != -1), ])
x <- data.frame(t(x))

mediumClinical =  data.frame(patient =as.character(unlist(clinicalFb[17,])), days2d = as.numeric(unlist(clinicalFb["patient.days_to_death", ])),
                             status = as.character(unlist(clinicalFb["patient.vital_status",])),
                             statusfollowup1 = as.character(unlist(clinicalFb["patient.follow_ups.follow_up.days_to_death",])),
                             statusfollowup2 = as.character(unlist(clinicalFb["patient.follow_ups.follow_up-2.days_to_death",])),
                             statusfollowup3 = as.character(unlist(clinicalFb["patient.follow_ups.follow_up-3.days_to_death",])),
                             statusfollowup4 = as.character(unlist(clinicalFb["patient.follow_ups.follow_up-4.days_to_death",])),
                             x,
                             stringsAsFactors = FALSE)
row.names(mediumClinical) <- mediumClinical$patient


mediumClinical <- mediumClinical[pppp,]
anyNA(mediumClinical)
which(is.na(mediumClinical[,3]))
View(mediumClinical[241,])
mediumClinical[241, which(!is.na(mediumClinical[241,]))]

mediumClinical <- mediumClinical[-241,]

View(mediumClinical[which((mediumClinical[,3] == "alive")), c("patient.days_to_last_followup","patient.follow_ups.follow_up.days_to_last_followup")])



####short -----


shortClinical =  data.frame(patient =as.character(unlist(clinicalFb[17,])),
                            days2d = as.numeric(unlist(clinicalFb["patient.days_to_death", ])),
                            status = as.character(unlist(clinicalFb["patient.vital_status",])),
                            lastfollowup = as.character(unlist(clinicalFb["patient.days_to_last_followup",])),
                            stringsAsFactors = FALSE)
row.names(shortClinical) <- shortClinical$patient
head(shortClinical)


shortClinical <- shortClinical[pppp,]
anyNA(shortClinical)
which(is.na(shortClinical[,3]))
shortClinical[241,]
# patient days2d status
# tcga-36-2533 tcga-36-2533     NA   <NA>
shortClinical <- shortClinical[-241,]

which(is.na(shortClinical[,2]))
which((shortClinical[,3] == "alive"))
shortClinical[ which(is.na(shortClinical[,2])), ]



### testing RTCGA ------
#https://cran.r-project.org/web/packages/survminer/vignettes/Informative_Survival_Plots.html
library(survminer)
library(RTCGA.clinical)
survivalTCGA(BRCA.clinical, OV.clinical,
             extract.cols = "admin.disease_code") -> BRCAOV.survInfo
library(survival)
fit <- survfit(Surv(times, patient.vital_status) ~ admin.disease_code,
               data = BRCAOV.survInfo)
# Visualize with survminer
ggsurvplot(fit, data = BRCAOV.survInfo, risk.table = TRUE)


View(OV.clinical)


survivalTCGA(OV.clinical) -> OV.survInfo

p5 <- intersect( tolower(as.character(OV.survInfo$bcr_patient_barcode)), as.character(unlist(patients)) )

OV.survInfo5 <- OV.survInfo[which(tolower(as.character(OV.survInfo$bcr_patient_barcode)) %in% p5),]
