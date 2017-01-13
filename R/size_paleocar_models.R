globalVariables(c("size"))
#' Get a RasterBrick of PaleoCAR model sizes
#'
#' This extracts the model size information 
#' from a set of PaleoCAR models and represents them as a RasterBrick with
#' layers for each reconstructed year.
#'
#' @param models A PaleoCAR batch model, as returned from \link{paleocar_models}.
#' @param prediction.years The set of years over which to generate model size rasters.
#' @param ... Further arguments to be passed to other functions.
#' @return A RasterBrick containing the model sizes per year.
size_paleocar_models <- function(models,
                                 prediction.years=NULL,
                                 ...){
  if(is.null(prediction.years)) prediction.years <- as.numeric(rownames(models[['reconstruction.matrix']]))
  
  if(!all(prediction.years %in% as.numeric(row.names(models[['reconstruction.matrix']])))){
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    prediction.years <- prediction.years[prediction.years %in% as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  models[['models']][,endYear := c(year[-1]-1,2000),by=cell]
  sizes <- data.matrix(models[['models']][,colnames(models[['reconstruction.matrix']]),with=F])
  sizes <- !is.na(sizes)
  sizes <- rowSums(sizes)
  models[['models']][,size:=sizes]
  
  sizes <- sapply(prediction.years,function(this.year){
    return(models[['models']][year<=this.year & endYear>=this.year,size])
  })
  
  if(class(models[["predictands"]]) %in% c("RasterBrick","RasterStack")){
    sizes <- raster::setValues(models[["predictands"]],sizes)
  }
  
  return(sizes)
}
