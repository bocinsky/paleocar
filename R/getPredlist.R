#' Get list of stepwise predictors from a reconstruction matrix
#'
#' This creates a list of available predictors along the years in 
#' the reconstruction matrix.
#'
#' @param predictor.matrix A matrix of predictors for a reconstruction.
#' Columns are predictors, rows are years.
#' @return A named list containing the names of available predictors 
#' for the years beginning with the name of element "i" and
#' ending with element "i+1"-1.
getPredlist <- function(reconstruction.matrix){
  reconstruction.matrix <- reconstruction.matrix
  reconstruction.matrix[!is.na(reconstruction.matrix)] <- 1
  reconstruction.matrix[is.na(reconstruction.matrix)] <- NA
  reconstruction.matrix <- unique(reconstruction.matrix)
#   predlist <- lapply(1:nrow(reconstruction.matrix),function(i){colnames(reconstruction.matrix)[which(reconstruction.matrix[i,]==1)]})
#   names(predlist) <- rownames(reconstruction.matrix)
  
  return(reconstruction.matrix)
}