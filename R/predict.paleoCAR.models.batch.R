#' Get a RasterBrick reconstruction from a PaleoCAR batch model
#'
#' This generates a reconstruction RasterBrick from a PaleoCAR batch model.
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
#' Thus, every reconstruction for a particular cell will have the same mean and 
#' variance over the calibration period.
#'
#' @param models A PaleoCAR batch model, as returned from \code{\link{paleoCAR.models.batch}}.
#' @param meanVarMatch Whether or not to perform mean-variance matching.
#' @param chained.meanVar Logical, indicating whether to chain mean-variance matching if performed. See \code{\link{predict.paleocar.models.batch}}.
#' @param prediction.years The set of years over which to generate reconstruction rasters. Optional.
#' @return A RasterBrick containing the predictions for each year.
predict.paleocar.models.batch <- function(models, meanVarMatch = TRUE, chained.meanVar=FALSE, prediction.years=NULL){
  if(is.null(prediction.years)) prediction.years <- as.numeric(rownames(models[['reconstruction.matrix']]))
  
  if(!all(prediction.years %in% as.numeric(row.names(models[['reconstruction.matrix']])))){
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    prediction.years <- prediction.years[prediction.years %in% as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  models[['reconstruction.matrix']] <- models[['reconstruction.matrix']][as.numeric(rownames(models[['reconstruction.matrix']])) %in% prediction.years,]

  newx <- data.table(cbind(1,models[['reconstruction.matrix']]))
  setnames(newx,c("Intercept",colnames(models[['reconstruction.matrix']])))
  rownames(newx) <- prediction.years
  newx.present <- !is.na(as.matrix(newx)[,-1])
  
  predictands.matrix <- t(raster::as.matrix(models[["predictands"]]))
  
  newx.calib <- as.matrix(cbind(data.table(Intercept=1),models[['predictor.matrix']]))
  rownames(newx.calib) <- rownames(models[['predictor.matrix']])
  
  
  predictions <- lapply(unique(models[['models']][['cell']]),function(this.cell){
    cat(this.cell,"\n")
    this.models <- models[['models']][cell==this.cell & year %in% prediction.years]
    coefficients <- rbindlist(lapply(this.models$coefs,function(x){data.table(matrix(data=x,ncol=length(x),byrow=T,dimnames=list(NA,names(x))))}),fill=T)
    this.newx <- newx[,names(coefficients),with=F]
    
    model.rows <- rep(1:nrow(coefficients),times=diff(c(this.models$year,tail(prediction.years,1)+1)))
    
    this.predictions <- rowSums(coefficients[model.rows]*this.newx, na.rm=T)
    
    if(chained.meanVar){
      coefficients.available <- newx.present[,names(coefficients),with=F]
      match.periods <- 
      
      
    }else if(meanVarMatch){
      this.newx.calib <- newx.calib[,names(coefficients)]
      calibration.years <- rownames(this.newx.calib)
      
      calibration.coefficients <- coefficients[,colnames(this.newx.calib),with=F]
      
      calibration.predictors <- as.matrix(this.newx[as.numeric(calibration.years)])[rep(1:length(calibration.years),nrow(calibration.coefficients)),]
      
      calibration.coefficients <- as.matrix(calibration.coefficients)[rep(1:nrow(calibration.coefficients),each=length(calibration.years)),]
      
      calibration.predictions <- rowSums(calibration.coefficients*calibration.predictors, na.rm=T)
      calibration.predictions <- do.call(rbind,split(calibration.predictions,rep(1:nrow(coefficients),each=length(calibration.years))))
      
      calibration.predictands <- predictands.matrix[,this.cell]
      
      scalars <- sd(calibration.predictands)/rowSds(calibration.predictions)
      transforms <- mean(calibration.predictands)-(scalars*rowMeans(calibration.predictions))
      
      this.predictions <- (this.predictions*scalars[model.rows]) + transforms[model.rows]
      
    }
    
    # Errors
    errors <- (this.models$CV)[model.rows]
    sizes <- (this.models$numPreds)[model.rows]
    
    return(list(predictions=this.predictions,errors=errors,sizes=sizes))
  })
  names(predictions) <- unique(models[['models']][['cell']])
  
  errors <- sapply(predictions,'[[','errors')
  sizes <- sapply(predictions,'[[','sizes')
  predictions <- sapply(predictions,'[[','predictions')
  
  if(class(models[["predictands"]]) %in% c("RasterBrick","RasterStack")){
    predictions <- setValues(models[["predictands"]],t(predictions))
    errors <- setValues(models[["predictands"]],t(errors))
    sizes <- setValues(models[["predictands"]],t(sizes))
  }
  return(list(predictions=predictions,errors=errors,sizes=sizes))
  
}
