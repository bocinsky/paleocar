#' Get a matrix of predictors for a reconstruction
#'
#' This subsets a chronologies matrix to a set of reconstruction years.
#'
#' @param chronologies A matrix of predictors to be subset.
#' Columns are predictors, rows are years.
#' @param reconstruction.years A vector of years to be reconstructed.
#' If missing, will use all years in \code{chronologies}.
#' @return A matrix of predictors for a reconstruction.
getReconstructionMatrix <- function(chronologies,reconstruction.years=NULL){
  if(is.null(reconstruction.years)){
    reconstruction.years <- chronologies$YEAR
  }
  
  reconstruction.matrix <- chronologies[chronologies$YEAR %in% reconstruction.years,]
  rownames(reconstruction.matrix) <- reconstruction.matrix$YEAR
  reconstruction.matrix <- reconstruction.matrix[,-1]
  return(as.matrix(reconstruction.matrix))  
}
