#' Get a matrix of predictors for a calibration
#'
#' This subsets a chronologies matrix to a set of calibration years.
#'
#' @param chronologies A matrix of predictors to be subset.
#' Columns are predictors, rows are years.
#' @param reconstruction.years A vector of years to be used for calibration
#' @return A matrix of predictors for a calibration
getPredictorMatrix <- function(chronologies,calibration.years){
  predictor.matrix <- chronologies[chronologies$YEAR %in% calibration.years,]
  predictor.matrix <- predictor.matrix[,t(complete.cases(t(predictor.matrix)))]
  rownames(predictor.matrix) <- predictor.matrix$YEAR
  predictor.matrix <- predictor.matrix[,-1]
  return(as.matrix(predictor.matrix))
}