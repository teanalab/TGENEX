#' createTensor
#' from boolean matrices to boolean 3- tensor
#'
#' @param mutationSmooth 
#' @param binaClinical 
#'
#' @return
#' @export
#'
#' @examples
#' load("data/mutationSmooth.RData")
#' load("data/patients.RData")
#' load("data/genes.RData")
#' load("data/binaClinical.RData")
#' MatrixOrderPxGxC <- createTensor(mutationSmooth,binaClinical)
#' save(MatrixOrderPxGxC, file="data/MatrixOrderPxGxC.RData")
createTensor <- function(mutationSmooth,binaClinical,patients)
{
  binaClinical <- binaClinical[ patients, ]
  mutationSmooth <- mutationSmooth[ patients, ]

  mutationSmoothNoNames <- mutationSmooth
  binaClinicalNoNames <- binaClinical

  names(mutationSmoothNoNames)  <- seq_len( dim(mutationSmooth)[2] )
  names(binaClinicalNoNames) <- seq_len( dim(binaClinical)[2] )

  ###add first table
  AmutationSmoothNoNames <- mutationSmoothNoNames
  for ( j in seq_along(row.names(AmutationSmoothNoNames)))
  {
    AmutationSmoothNoNames[j,] <- binaClinicalNoNames[j,1] * AmutationSmoothNoNames[j,]
  }
  #write.csv(AmutationSmoothNoNames,file = "temp/1-3-tensorMatrices-01.csv")
  MatrixOrderPxGxC <- AmutationSmoothNoNames

  for ( i in (2:length(names(binaClinicalNoNames))))
  {
    AmutationSmoothNoNames <- mutationSmoothNoNames
    for ( j in seq_along(row.names(AmutationSmoothNoNames)))
    {
      AmutationSmoothNoNames[j,] <- binaClinicalNoNames[j,i] * AmutationSmoothNoNames[j,]
    }
    ##CONCATANATE MATRICES
    # write.csv(AmutationSmoothNoNames,file = paste("temp/20-tensorMatrices-",
    #                                 i,".csv", sep=''))
    MatrixOrderPxGxC <- cbind(MatrixOrderPxGxC, AmutationSmoothNoNames)
  }
  MatrixOrderPxGxC
}


