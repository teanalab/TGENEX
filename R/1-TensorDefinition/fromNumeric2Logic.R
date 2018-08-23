
#' from Numeric to Logic
#'
#' @param aDataFrame
#'
#' @return
#' @export
#'
#' @examples
#' boolMutation <- fromNumeric2Logic(boolMutation)
fromNumeric2Logic <- function(aDataFrame)
{
  bDataFrame <- data.frame(matrix(F,nrow=dim(aDataFrame)[[1]],ncol=dim(aDataFrame)[[2]]),
                           row.names = row.names(aDataFrame), stringsAsFactors = FALSE)
  for ( j in seq_along(bDataFrame)){
    for ( i in seq_along(bDataFrame[[1]])){
      bDataFrame[i,j] <- as.logical(aDataFrame[i,j])
    }
  }
  names(bDataFrame) <- names(aDataFrame)
  bDataFrame
}


#' from Numeric to Logic
#'
#' @param aDataFrame
#'
#' @return
#' @export
#'
#' @examples
#' boolMutation <- fromNumeric2Logic(boolMutation)
fromNumeric2LogicFast <- function(aDataFrame)
{
  bDataFrame <- as.data.frame(lapply(aDataFrame, as.logical), stringsAsFactors = FALSE)
  row.names(bDataFrame) <- row.names(aDataFrame)
  names(bDataFrame) <- names(aDataFrame)
  bDataFrame
}


