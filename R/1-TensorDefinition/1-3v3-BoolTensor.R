# version 3 (5/20/2018)
#1) using firebrowse and cBiolite data

rm(list = ls())
load(file="data/boolClinical.RData")
load("data/boolMutation.RData")

#Get all needed data
LoadMyData <- funtion()
{
  #Get Boolean Matrices
  load("temp/10-BoolMatrices_v2.RData")
  rm(list=ls()[-which(ls() %in% c("boolClinical", "boolMutation"))])

  # libs<-c("Packages.R")
  # libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
  # sapply(libs, function(u) {source(u)})

  save.image(file="temp/20-data.RData")
}


#from boolean matrices to boolean 3- tensor
createTensor <-function(boolMutation,boolClinical)
{
  #make sure the matrices have the same patients in the same order
  #sort data by patient
  boolClinical <- boolClinical[ order(row.names(boolClinical)), ]
  boolMutation <- boolMutation[ order(row.names(boolMutation)), ]

  bC <- row.names(boolClinical)
  bM <- row.names(boolMutation)

  #same patients
  length(bC)
  length(bM)
  length(intersect(bC,bM))

  boolMutationNoNames <- boolMutation
  boolClinicalNoNames <- boolClinical

  names(boolMutationNoNames)  <- seq_len( dim(boolMutation)[2] )
  names(boolClinicalNoNames) <- seq_len( dim(boolClinical)[2] )

  #same patients
  View(cbind(row.names(boolClinicalNoNames),row.names(boolMutationNoNames)))

  ###add first table
  AboolMutationNoNames <- boolMutationNoNames
  for ( j in seq_along(row.names(AboolMutationNoNames)))
  {
    AboolMutationNoNames[j,] <- boolClinicalNoNames[j,1] & as.logical(AboolMutationNoNames[j,])
  }
  write.csv(AboolMutationNoNames,file = "temp/20-tensorMatrices-01.csv")
  MatrixOrderPxGxC <- AboolMutationNoNames


  for ( i in (2:length(names(boolClinicalNoNames))))
  {
    AboolMutationNoNames <- boolMutationNoNames
    for ( j in seq_along(row.names(AboolMutationNoNames)))
    {
      AboolMutationNoNames[j,] <- boolClinicalNoNames[j,i] & as.logical(AboolMutationNoNames[j,])
    }
    ##CONCATANATE MATRICES
    write.csv(AboolMutationNoNames,file = paste("temp/20-tensorMatrices-",
                                    i,".csv", sep=''))
    MatrixOrderPxGxC <- cbind(MatrixOrderPxGxC,AboolMutationNoNames)
  }
  save.image(file = "temp/1-3v3.RData")
  load(file = "temp/1-3v3.RData")
}


fromLogic2Numeric <- function(MatrixOrderPxGxC)
{
  #to integer
  IntegerPxGxC <- as.data.frame(lapply(MatrixOrderPxGxC, as.numeric), stringsAsFactors = FALSE)
  save.image(file = "temp/20-1and2_86clinicalBRCA_v2.RData")
}

write.csv2(IntegerPxGxC, file = "temp/1-3v3-IntegerPxGxC.csv",row.names = FALSE, col.names = FALSE)
write.csv2(MatrixOrderPxGxC, file = "temp/20-MatrixOrderPxGxC.csv",row.names = FALSE)

tIPGC <- t(IntegerPxGxC)
write.csv2(tIPGC, file = "temp/20-tIntegerPxGxC.csv",row.names = FALSE)
write.csv2(t(MatrixOrderPxGxC), file = "temp/20-tMatrixOrderPxGxC.csv",row.names = FALSE)


load(file = "temp/20-1and2_86clinicalBRCA_v2.RData")

#test
boolMutation[1,1]
boolClinical[1,1]
MatrixOrderPxGxC[1,1,1]
sum(IntegerPxGxC)
