#' Get a matrix of predictors for a reconstruction
#'
#' This subsets a chronologies matrix to a set of reconstruction years.
#'
#' @param chronologies A matrix of predictors to be subset.
#' Columns are predictors, rows are years.
#' @param reconstruction.years A vector of years to be reconstructed.
#' If missing, will use all years in \code{chronologies}.
#' @param min.width integer, indicating the minimum number of tree-ring samples allowed for that year of a chronology to be valid.
#' @return A matrix of predictors for a reconstruction.
getReconstructionMatrix <- function(chronologies,reconstruction.years=NULL,min.width=NULL){
  YEARS <- as.numeric(rownames(chronologies[['widths']]))
  
  if(is.null(reconstruction.years)){
    reconstruction.years <- YEARS
  }
  
  if(!is.null(min.width)){
    chronologies[['widths']][chronologies[['depths']]<min.width] <- NA
  }
  
  reconstruction.matrix <- chronologies[['widths']][YEARS %in% reconstruction.years,]
#   rownames(reconstruction.matrix) <- reconstruction.matrix$YEAR
#   reconstruction.matrix <- reconstruction.matrix[,-1]
  return(as.matrix(reconstruction.matrix))  
}
