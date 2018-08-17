

#' from Logic to Numeric
#'
#' @param aDataFrame
#'
#' @return
#' @export
#'
#' @examples
#' binaMutation <- fromLogic2Numeric(boolMutation)
fromLogic2Numeric <- function(aDataFrame)
{
  bDataFrame <- as.data.frame(lapply(aDataFrame, as.numeric), stringsAsFactors = FALSE)
  row.names(bDataFrame) <- row.names(aDataFrame)
  names(bDataFrame) <- names(aDataFrame)
  bDataFrame
}

