#' Get a RasterBrick of reconstruction LOOCV errors
#'
#' This extracts the leave-one-out cross-validation (LOOCV) error information 
#' from a set of PaleoCAR models and represents them as a RasterBrick with
#' layers for each reconstructed year. The LOOCV statistic is the sum of cross-validated 
#' squared errors for the model over the calibration period. This method returns the average errors
#' by dividing the LOOCV statistic by the number of predictand observations, and then taking the square root.
#' Resulting error is delivered in the units of the predictand.
#'
#' @param models A PaleoCAR batch model, as returned from \link{paleocar_models_batch}.
#' @param prediction.years The set of years over which to generate error rasters.
#' @import raster
#' @return A RasterBrick containing the LOOCV statistics per year.
errors_paleocar_models_batch <- function(models, prediction.years=NULL){
  if(is.null(prediction.years)) prediction.years <- as.numeric(rownames(models[['reconstruction.matrix']]))
  
  if(!all(prediction.years %in% as.numeric(row.names(models[['reconstruction.matrix']])))){
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    prediction.years <- prediction.years[prediction.years %in% as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  models[['models']][,endYear := c(year[-1]-1,2000),by=cell]
  
  errors <- sapply(prediction.years,function(this.year){ 
    return(models[['models']][year<=this.year & endYear>=this.year,CV])
  })
  
  if(class(models[["predictands"]]) %in% c("RasterBrick","RasterStack")){
    errors <- raster::setValues(models[["predictands"]],errors)
  }
  
  errors <- sqrt(errors/nrow(models[['predictor.matrix']]))
  
  return(errors)
}
