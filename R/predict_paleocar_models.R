#' Get a reconstruction from a PaleoCAR simple model
#'
#' This generates a reconstruction vector from a PaleoCAR simple model.
#' 
#' Optionally, mean-variance matching (scaling/transformation) can be performed.
#' Let \bold{\eqn{X_c}} be the vector of calibration values over the calibration period,
#' \bold{\eqn{X_r}} be the vector of reconstructed values over the calibration period,
#' \bold{\eqn{A_r}} be the vector of all reconstructed values (of which \bold{\emph{X_r}} 
#' is a part), and \bold{\emph{A_r*}} vector of all reconstructed values after 
#' mean and variance matching. Then,
#' 
#' \bold{\emph{A_r*}} = \eqn{\alpha} x \bold{\eqn{A_r}}) + \bold{\eqn{\beta}},
#' 
#' where, \eqn{\alpha} = \eqn{\sigma}(\bold{\emph{X_c}})/\eqn{\sigma}(\bold{\emph{X_r}}),
#' 
#' or a scalar that is the ratio of the standard deviations of the calibration and 
#' reconstructed vectors over the calibration period,
#' 
#' and, \bold{\eqn{\beta}} = \eqn{\mu}(\bold{\emph{X_c}})-[\eqn{\alpha} x \eqn{\mu}(\bold{\emph{X_r}})],
#' 
#' or a transformation to the mean of the calibration vector corrected by the 
#' scaled mean of the reconstructed data over the calibration period.
#' Thus, every reconstruction through time for will have the same mean and 
#' variance over the calibration period.
#'
#' @param models A PaleoCAR simple model, as returned from \code{\link{paleocar_models}}.
#' @param meanVarMatch Whether or not to perform mean-variance matching.
#' @param prediction.years The set of years over which to generate a reconstructions. Optional.
#' @return A numeric vector containing the predictions for each year.
predict_paleocar_models <- function(models, meanVarMatch = TRUE, years=NULL){
  if(is.null(years)) years <- as.numeric(rownames(models[['reconstruction.matrix']]))
  
  if(!all(years %in% as.numeric(row.names(models[['reconstruction.matrix']])))){
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    years <- years[years %in% as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  coefficients <- data.matrix(models[['models']][,c("Intercept",colnames(models[['reconstruction.matrix']])),with=F])
  
  coefficients[is.na(coefficients)] <- 0
  rownames(coefficients) <- models[['models']]$year
  
  predictand.mean <- mean(models[['predictand']])
  predictand.sd <- sd(models[['predictand']])
  
  breaks <- as.numeric(rownames(coefficients))
  coefficients.expand <- matrix_expand(coefficients,1:2000)
  predictions <- rowSums(coefficients.expand * cbind(1,models[['reconstruction.matrix']]), na.rm=T)
  
  if(meanVarMatch){
    # Mean and variance matching over calibration period
    calibration.predictions <- apply(coefficients,1,function(coefs){rowSums(sweep(cbind(1,models[['predictor.matrix']]),MARGIN=2,coefs,`*`))})
    calibration.predictions.mean <- apply(calibration.predictions,2,mean)
    calibration.predictions.sd <- apply(calibration.predictions,2,sd)
    
    scalar <- predictand.sd/calibration.predictions.sd
    scalar.expand <- rep(scalar,diff(c(breaks,2001)))
    predictions.scaled <- predictions * scalar.expand
    
    predictions.scaled.meaned <- predictions.scaled + rep(predictand.mean-(scalar*calibration.predictions.mean),diff(c(breaks,2001)))
    
    return(predictions.scaled.meaned)
  }

  return(predictions)

}
