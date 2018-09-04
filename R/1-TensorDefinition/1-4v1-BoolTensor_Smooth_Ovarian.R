# version 1 (9/4/2018)

rm(list = ls())

#Get all needed data
LoadMyData <- function()
{
  load("data/binaClinicalSmooth.Rd", envir = globalenv())
  load("data/matrixMutationSmooth.Rd", envir = globalenv())
  load("loadls.RData", envir = globalenv())
}
LoadMyData()

#from boolean matrices to boolean 3- tensor
createTensor <-function(binaMutation,binaClinical)
{
  bC <- row.names(binaClinical)
  bM <- row.names(binaMutation)

  #same patients
  patients <- as.character(unlist( intersect(bC,bM) ))
  binaClinical <- binaClinical[ patients, ]
  binaMutation <- binaMutation[ patients, ]


  binaMutationNoNames <- binaMutation
  binaClinicalNoNames <- binaClinical

  names(binaMutationNoNames)  <- seq_len( dim(binaMutation)[2] )
  names(binaClinicalNoNames) <- seq_len( dim(binaClinical)[2] )



  ###add first table
  AbinaMutationNoNames <- binaMutationNoNames
  for ( j in seq_along(row.names(AbinaMutationNoNames)))
  {
    AbinaMutationNoNames[j,] <- binaClinicalNoNames[j,1] * AbinaMutationNoNames[j,]
  }
  #write.csv(AbinaMutationNoNames,file = "temp/1-3-tensorMatrices-01.csv")
  MatrixOrderPxGxC <- AbinaMutationNoNames

  for ( i in (2:length(names(binaClinicalNoNames))))
  {
    AbinaMutationNoNames <- binaMutationNoNames
    for ( j in seq_along(row.names(AbinaMutationNoNames)))
    {
      AbinaMutationNoNames[j,] <- binaClinicalNoNames[j,i] * AbinaMutationNoNames[j,]
    }
    ##CONCATANATE MATRICES
    # write.csv(AbinaMutationNoNames,file = paste("temp/20-tensorMatrices-",
    #                                 i,".csv", sep=''))
    MatrixOrderPxGxC <- cbind(MatrixOrderPxGxC, AbinaMutationNoNames)
  }
  save.image(file = "temp/1-4v1.RData")

}

createTensor(binaMutation, binaClinical)

save(MatrixOrderPxGxC, file="data/MatrixSmoothOrderPxGxC.RData")

#test -------
binaMutation[1,1]
binaClinical[1,1]
MatrixOrderPxGxC[100:105,100:105]
sum(MatrixOrderPxGxC)





