#' Get a matrix of predictors for a calibration
#'
#' This subsets a chronologies matrix to a set of calibration years.
#'
#' @param chronologies A matrix of predictors to be subset. Should be output from FedData::get_itrdb().
#' Columns are predictors, rows are years.
#' @param calibration.years A vector of years to be used for calibration
#' @param min.width integer, indicating the minimum number of tree-ring samples allowed for that year of a chronology to be valid.
#' @return A matrix of predictors for a calibration
getPredictorMatrix <- function(chronologies,calibration.years,min.width=NULL){
  YEARS <- as.numeric(rownames(chronologies[['widths']]))
  
  if(!is.null(min.width)){
    chronologies[['widths']][chronologies[['depths']]<min.width] <- NA
  }
  
  predictor.matrix <- chronologies[['widths']][YEARS %in% calibration.years,]
  predictor.matrix <- predictor.matrix[,t(complete.cases(t(predictor.matrix)))]
  predictor.matrix <- predictor.matrix[,!duplicated(t(predictor.matrix))]
  return(as.matrix(predictor.matrix))
}