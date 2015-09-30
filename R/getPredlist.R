#' Get list of stepwise predictors from a reconstruction matrix
#'
#' This creates a list of available predictors along the years in 
#' the reconstruction matrix.
#'
#' @param predictor.matrix A matrix of predictors for a reconstruction.
#' Columns are predictors, rows are years.
#' @return A matrix of unique combinations of predictors through time.
getPredlist <- function(reconstruction.matrix){
  reconstruction.matrix.present <- !is.na(reconstruction.matrix)
  # reconstruction.matrix.present[!reconstruction.matrix.present] <- NA
  reconstruction.matrix.present <- unique(reconstruction.matrix.present)

  return(reconstruction.matrix.present)
}