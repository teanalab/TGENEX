# Construct Binary Matrices from RTCGA

#1- Example how to get Breast Cancer mutation data using RTCGA
source('../code/R/1-createBoolMutation.R', echo=F)
load("data/loadls.RData")
mutationFromRTCGA(disease="BRCA")

#2- Example how to get Breast Cancer clinical data using RTCGA
source('~/10-Research/Current/CLINGEN/9source/CLIGEN/publicPackage/CLIGEN/R/1-createBoolClinical.R', echo=F)
clinicalFromRTCGA(disease="BRCA")

