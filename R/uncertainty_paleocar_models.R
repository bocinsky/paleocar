globalVariables(c("year","endYear", "cell"))
#' Get reconstruction LOOCV uncertainty
#'
#' This extracts the leave-one-out cross-validation (LOOCV) uncertainty information 
#' from a set of PaleoCAR models and represents them as a RasterBrick with
#' layers for each reconstructed year. The LOOCV statistic is the sum of cross-validated 
#' squared errors for the model over the calibration period. This method returns the average uncertainty
#' by dividing the LOOCV statistic by the number of predictand observations, and then taking the square root.
#' Resulting uncertainty is delivered in the units of the predictand.
#'
#' @param models A PaleoCAR batch model, as returned from \link{paleocar_models}.
#' @param prediction.years The set of years over which to generate error estimates.
#' @param ... Further arguments to be passed to other functions.
#' @return The LOOCV statistics per year.
#' @importFrom data.table :=
#' @export
uncertainty_paleocar_models <- function(models,
                                   prediction.years = NULL,
                                   ...){
  if(is.null(prediction.years)) prediction.years <- as.numeric(rownames(models$reconstruction.matrix))
  
  if(!all(prediction.years %in% as.numeric(row.names(models$reconstruction.matrix)))){
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    prediction.years <- prediction.years[prediction.years %in% as.numeric(rownames(models$reconstruction.matrix))]
  }
  
  models$models[,endYear := c(year[-1]-1,tail(as.numeric(rownames(models$reconstruction.matrix)),1)), by=cell]
  
  errors <- sapply(prediction.years,
                   function(this.year){ 
                     return(models$models[year<=this.year & endYear>=this.year,CV])
                   })
  
  if(is.matrix(errors)){
    errors <- t(errors)
  }
  
  if(is.vector(errors)){
    errors <- matrix(errors)
  }
  
  if(class(models$predictands) %in% c("RasterBrick","RasterStack")){
    errors <- raster::setValues(models$predictands,t(errors))
    names(errors) <- prediction.years
  }
  
  errors <- sqrt(errors/nrow(models$predictor.matrix))
  
  if(is.vector(models$predictands)){
    errors <- as.vector(errors)
    names(errors) <- prediction.years
  } 
  
  return(errors)
}
