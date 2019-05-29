globalVariables(c("year","endYear", "cell"))
#' Get reconstruction standard errors of the models, through time
#'
#' This extracts the standard errors 
#' from a set of PaleoCAR models and represents them as a RasterBrick with
#' layers for each reconstructed year.
#'
#' @param models A PaleoCAR batch model, as returned from \link{paleocar_models}.
#' @param prediction.years The set of years over which to generate error estimates.
#' @param ... Further arguments to be passed to other functions.
#' @return The standard error per year.
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
    na.errors <- matrix(data=NA,
                        nrow=length(prediction.years),
                        ncol=raster::ncell(models$predictands))
    
    na.errors <- raster::setValues(models$predictands,
                                values = t(na.errors))
    
    na.errors[models$models$cell %>% unique()] <- t(errors)
    
    errors <- na.errors
    
    names(errors) <- prediction.years
  }
  
  errors <- sqrt(errors/nrow(models$predictor.matrix))
  
  if(is.vector(models$predictands)){
    errors <- as.vector(errors)
    names(errors) <- prediction.years
  }
  
  return(errors)
}
