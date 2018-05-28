# 5 components after selecting 15
# Component 5 has only 2 patients and ruins the analysis
# but makes not huge difference


#analize the data
rm(list = ls())

load("temp/FinalFrom71.RData")

#source("R/plotSurvivalUpdate.R")
#load("results/BRCA15.RData")

#load components by patient, generated with 41-AnalizeCP_BRCA

#save.image(file="temp/components15.RData")
#load("temp/components15.RData")

#A <- Boolcp_3comp$A #patients vs components

patientVsComponents <- as.data.frame(A)

patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)

lcmps <- c(1, 6, 9, 11) # c(1, 5, 6, 9, 11)

patientVsComponents <- patientVsComponents [which(patientVsComponents$whichmax %in% lcmps), ]


patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
                                component=patientVsComponents$whichmax)

library(plyr)
#frequency of components
count(patientAfiliation, 'component')


#attach survival data
#clinical<-read.table(file="../../12Datasets/BRCA/clinical.txt", sep = "\t",
#                     header = TRUE, stringsAsFactors = FALSE)

#remove unwanted columns
#[10] "Integrated.Clusters..no.exp."
#[11] "Integrated.Clusters..unsup.exp."  "Integrated.Clusters..with.PAM50."
#[13] "MIRNA.Cluster"
#[16] "Methylation.Cluster"
#[21] "PAM50.subtype"
#[24] "RPPA.Cluster"
#[25] "SigClust.Intrinsic.mRNA"          "SigClust.Unsupervised.mRNA"


#include survival columns
#[19] "Overall.Survival..Months."        "Overall.Survival.Status"
#[27] "Survival.Data.Form"


#select 5 components




survivalClinical = clinical[,c(1,19,20)]

survivalClinical = merge(patientAfiliation, survivalClinical, by.x="patient", by.y= "Patient.ID")

count(survivalClinical, 'Overall.Survival.Status')


#load my libraries
#libs<-c("Packages.R", "permutation.R", "clustering.R", "Plots.R", "survival_analysis.R")
#libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
#sapply(libs, function(u) {source(u)})


groups = factor(survivalClinical$component)



survivalData = data.frame(PatientID=survivalClinical$patient,
                          Survival=survivalClinical$Overall.Survival..Months.,
                          Death=(survivalClinical$Overall.Survival.Status=="DECEASED"))

labelsGroups <- c(paste("group 1 (",table(groups)[2], ")", sep="" ),
                  paste("group 2 (",table(groups)[1], ")", sep="" ),
                  paste("group 3 (",table(groups)[3], ")", sep="" ) )
colorsLines <- c("black","red","blue")
colorsLabels <- c("red","black","blue")
 pdf( file = "temp/74_Survival4from15components.pdf",  onefile = TRUE, width = 9, height = 7)
 plot.SurvivalKN (groups, "mainTitle", survivalData,"bottomright",
                  labelClu = labelsGroups,
                  colorsP = colorsLines,
                  colorsL = colorsLabels,
                  centerT = 0.4 )
 dev.off()



