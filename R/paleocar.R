#' Fit PaleoCAR models to a RasterBrick of predictands, and generate prediction and uncertainty RasterBricks.
#'
#' This is a wrapper function to the four primary methods in the PaleoCAR package, specifically for
#' processing predictands held in a RasterBrick.
#'
#' @param chronologies An ITRDB object, as in from FedData::get_itrdb.
#' @param predictands A RasterBrick or RasterStack of the numeric predictand (response) variable.
#' @param calibration.years An integer vector of years corresponding to the layers in the \code{predictands} brick.
#' @param prediction.years An optional integer vector of years for the reconstruction.
#' If missing, defaults to the total years present in \code{chronologies}.
#' @param label A character label for the reconstruction, for saving.
#' @param out.dir The directory to which output is to be saved.
#' @param verbose Logical, display status messages during run.
#' @param force.redo Logical, should all computations be re-computed?
#' @param ... Further arguments to be passed to other functions.
#' @return A named list containing
#' \itemize{
#'   \item{\code{models}  The PaleoCAR models, as computed by \code{\link{paleocar_models}}.}
#'   \item{\code{predictions}  The PaleoCAR reconstruction, as computed by \code{\link{predict_paleocar_models}}.}
#' }
#' @export
paleocar <- function(chronologies,
                     predictands,
                     calibration.years,
                     prediction.years = NULL,
                     label,
                     out.dir = "./OUTPUT/",
                     force.redo = F,
                     verbose = F,
                     ...){
  t <- Sys.time()
  
  if(verbose) cat("\nCalculating all models\n")
  if(!force.redo & file.exists(paste0(out.dir,label,'.models.Rds'))){
    models <- readr::read_rds(paste0(out.dir,label,'.models.Rds'))
    
  }else{
    unlink(paste0(out.dir,label,".models.Rds"), 
           recursive = TRUE, 
           force = TRUE)
    
    models <- paleocar_models(chronologies = chronologies,
                              predictands = predictands,
                              calibration.years = calibration.years,
                              prediction.years = prediction.years,
                              verbose = verbose,
                              ...)
    
    readr::write_rds(models,
                     file = paste0(out.dir,label,".models.Rds"),
                     compress = "gz")
  }
  
  if(verbose) cat("\nGenerating prediction\n")
  if(!force.redo & file.exists(paste0(out.dir,label,".prediction.Rds"))){
    recon <- readr::read_rds(paste0(out.dir,label,".prediction.Rds"))
    
  }else{
    unlink(paste0(out.dir,label,".prediction.Rds"), 
           recursive = TRUE, 
           force = TRUE)
    
    recon <- predict_paleocar_models(models = models)
    readr::write_rds(recon,
                     file = paste0(out.dir,label,".prediction.Rds"),
                     compress = "gz")
    
  }
  
  if(verbose) message("\nThe entire reconstruction took ", round(difftime(Sys.time(),t,units='mins'), digits = 2)," minutes")
  
  return(
    list(models = models,
         predictions = recon)
  )
  
}
