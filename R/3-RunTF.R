rm(list = ls())

#Load tensor and matrices
# load("dataExample/noFiltered/MatrixOrderPxGxC.RData")
# load("dataExample/noFiltered/boolMutation.RData")
load("dataExample/Filtered/MatrixOrderPxGxC.RData")
load("dataExample/Filtered/binaMatrices.RData")
load("dataExample/Filtered/mutationSmooth.RData")

load("data/loadls.RData")
labelsPA <- patients
labelsGenes <- names(binaMutation)
labelsClinical <- names(binaClinical)


#manuallychange the value of r in line 30
r=10

library(ThreeWay) #lib for tensor factorizarion

con <- file(paste("temp/log_r",r,".log", sep='') )
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")

bcCP <- CP(MatrixOrderPxGxC, labelsPA, labelsGenes, labelsClinical)
482
499
32
0
0
0
0
0
10
0
1e-6
0
10000
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0

# Restore output to console
sink()
sink(type="message")

# And look at the log...
#cat(readLines("test.log"), sep="\n")

save(bcCP, file=paste("dataExample/Filtered/bcCPk",r,".Rd", sep='') )

patientsF <- bcCP$A  #p x k
genesF <- bcCP$B #genes x k
clinicalF <- bcCP$C #clinical x k.


LoadFiles <- function(r){
  load(paste("dataExample/Filtered/bcCPk",r,".Rd", sep=''), envir = globalenv())
  patientsF <- bcCP$A  #p x k
  genesF <- bcCP$B #genes x k
  clinicalF <- bcCP$C #clinical x k.
  save(list=c("patientsF", "genesF", "genesF"), file=paste("temp/factors_r",r,".RData", sep='') )
}

LoadFiles(r)

