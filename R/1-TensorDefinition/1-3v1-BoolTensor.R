# version 1 (Jul-20-2017)

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
  names(boolClinicalNoNames) <- seq_len( dim(boolClinical)[2])

  #same patients
  View(cbind(row.names(boolClinicalNoNames),row.names(boolMutationNoNames)))

  ###add first table
  AboolMutationNoNames <- boolMutationNoNames
  for ( j in seq_along(row.names(AboolMutationNoNames)))
  {
    AboolMutationNoNames[j,] <- boolClinicalNoNames[j,1] & as.logical(AboolMutationNoNames[j,])
  }
  write.csv(AboolMutationNoNames,file = "tensorMatrices/1.csv")
  MatrixOrderPxGxC <- AboolMutationNoNames


  for ( i in (2:length(names(boolClinicalNoNames))))
  {
    AboolMutationNoNames <- boolMutationNoNames
    for ( j in seq_along(row.names(AboolMutationNoNames)))
    {
      AboolMutationNoNames[j,] <- boolClinicalNoNames[j,i] & as.logical(AboolMutationNoNames[j,])
    }
    ##CONCATANATE MATRICES
    write.csv(AboolMutationNoNames,file = paste("tensorMatrices/",
                                    i,".csv", sep=''))
    MatrixOrderPxGxC <- cbind(MatrixOrderPxGxC,AboolMutationNoNames)
  }

  save(list=c("boolClinical","boolMutation","MatrixOrderPxGxC"),
       file = "1and2_86clinicalBRCA.RData")
}



fromLogic2Numeric <- function(MatrixOrderPxGxC)
{
  #integer
  IntegerPxGxC <- as.data.frame(as.numeric(MatrixOrderPxGxC),
                                stringsAsFactors = FALSE)

}
