#' Computation of CAR scores over a matrix of predictands
#'
#' Takes a matrix of predictands and a matrix of predictors (each predictand/predictor gets a column)
#' and returns a list of CAR scores, one for each columns in the predictand.matrix.
#'
#' @param predictand.matrix A matrix of predictand (response variables). Columns are predictands
#' @param predictor.matrix A matrix of predictors. Columns are predictors, number of rows must equal number of layers in \code{predictand.matrix}.
#' @return A matrix of carscores.
carscore.batch <- function(predictand.matrix, predictor.matrix){
  
  carscores <- apply(predictand.matrix,2,function(predictands){
    if(any(is.na(predictands)){
      out <- predictor.matrix[1,]
      out[] <- NA
      return(out)
    }
    carscore(Ytrain=predictands,Xtrain=predictor.matrix,verbose=F)
  })
  
  return(carscores)
}