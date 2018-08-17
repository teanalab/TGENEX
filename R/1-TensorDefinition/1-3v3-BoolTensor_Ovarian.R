# version 3 (8/15/2018)

rm(list = ls())

#Get all needed data
LoadMyData <- function()
{
  load(file="data/boolClinical.RData", envir = globalenv())
  load("data/boolMutation.RData", envir = globalenv())
  load("loadls.RData", envir = globalenv())
}
LoadMyData()

#from boolean matrices to boolean 3- tensor
createTensor <-function(boolMutation,boolClinical)
{
  bC <- row.names(boolClinical)
  bM <- row.names(boolMutation)

  #same patients
  patients <- as.character(unlist( intersect(bC,bM) ))
  boolClinical <- boolClinical[ patients, ]
  boolMutation <- boolMutation[ patients, ]


  boolMutationNoNames <- boolMutation
  boolClinicalNoNames <- boolClinical

  names(boolMutationNoNames)  <- seq_len( dim(boolMutation)[2] )
  names(boolClinicalNoNames) <- seq_len( dim(boolClinical)[2] )


  ###add first table
  AboolMutationNoNames <- boolMutationNoNames
  for ( j in seq_along(row.names(AboolMutationNoNames)))
  {
    AboolMutationNoNames[j,] <- as.logical(boolClinicalNoNames[j,1]) & as.logical(AboolMutationNoNames[j,])
  }
  #write.csv(AboolMutationNoNames,file = "temp/1-3-tensorMatrices-01.csv")
  MatrixOrderPxGxC <- AboolMutationNoNames

  for ( i in (2:length(names(boolClinicalNoNames))))
  {
    AboolMutationNoNames <- boolMutationNoNames
    for ( j in seq_along(row.names(AboolMutationNoNames)))
    {
      AboolMutationNoNames[j,] <- as.logical(boolClinicalNoNames[j,i]) & as.logical(AboolMutationNoNames[j,])
    }
    ##CONCATANATE MATRICES
    # write.csv(AboolMutationNoNames,file = paste("temp/20-tensorMatrices-",
    #                                 i,".csv", sep=''))
    MatrixOrderPxGxC <- cbind(MatrixOrderPxGxC, AboolMutationNoNames)
  }
  save.image(file = "temp/1-3v3.RData")
  load(file = "temp/1-3v3.RData")
}

createTensor(boolMutation, boolClinical)

IntegerPxGxC <- fromLogic2Numeric(MatrixOrderPxGxC)

save(IntegerPxGxC, file ="data/IntegerPxGxC.RData")
save(MatrixOrderPxGxC, file="data/MatrixOrderPxGxC.RData")

#write csv files for rubik matlab ------------

load(file = "data/IntegerPxGxC.RData")

#write.csv2(IntegerPxGxC, file = "temp/1-3v3-IntegerPxGxC.csv",row.names = FALSE, col.names = FALSE)
#write.csv2(MatrixOrderPxGxC, file = "temp/20-MatrixOrderPxGxC.csv",row.names = FALSE)

tIPGC <- t(IntegerPxGxC)
write.csv2(tIPGC, file = "temp/1-3-OV-tIntegerPxGxC.csv",row.names = FALSE)
#write.csv2(t(MatrixOrderPxGxC), file = "temp/20-tMatrixOrderPxGxC.csv",row.names = FALSE)


load(file = "temp/20-1and2_86clinicalBRCA_v2.RData")

#test -------
boolMutation[1,1]
boolClinical[1,1]
MatrixOrderPxGxC[10:15,10:15]
sum(IntegerPxGxC)
IntegerPxGxC[1:5,1:5]




