#' Fit PaleoCAR models to a RasterBrick of predictands, and generate prediction and uncertainty RasterBricks.
#'
#' This is a wrapper function to the four primary methods in the PaleoCAR package, specifically for
#' processing predictands held in a RasterBrick. Additionally, it allows for the specification of
#' a \code{floor} and/or \code{ceiling} corrections for a reconstruction (e.g., of precipitation).
#' It returns a list containing the models, predictions, and cross-validated uncertainty estimates.
#'
#' @param chronologies An ITRDB object, as in from FedData::get_itrdb.
#' @param predictands A RasterBrick or RasterStack of the numeric predictand (response) variable.
#' @param calibration.years An integer vector of years corresponding to the layers in the \code{predictands} brick.
#' @param prediction.years An optional integer vector of years for the reconstruction.
#' If missing, defaults to the total years present in \code{chronologies}.
#' @param label A character label for the reconstruction, for saving.
#' @param out.dir The directory to which output is to be saved.
#' @param return.predictions Logical, should the predictions be generated?
#' @param return.uncertainty Logical, should the uncertainty estimates be generated?
#' @param return.objects Logical, should objects be returned (as opposed to only writing them to disk)?
#' @param verbose Logical, display status messages during run.
#' @param force.redo Logical, should all computations be re-computed?
#' @param ... Further arguments to be passed to other functions.
#' @return A named list containing
#' \itemize{
#'   \item{\code{models}  The PaleoCAR models, as computed by \code{\link{paleocar_models}}.}
#'   \item{\code{recon}  The PaleoCAR reconstruction, as computed by \code{\link{predict_paleocar_models}}.}
#'   \item{\code{errors}  The PaleoCAR reconstruction average LOOCV error, as computed by \code{\link{uncertainty_paleocar_models}}.}
#'   \item{\code{sizes}  The PaleoCAR model sizes, as computed by \code{\link{size_paleocar_models}}.}
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
                     return.predictions = T,
                     return.uncertainty = T,
                     return.objects = T,
                     ...){
  t <- Sys.time()

  if(verbose) cat("\nCalculating all models\n")
  if(!force.redo & file.exists(paste0(out.dir,label,'.models.Rds'))){
    models <- readr::read_rds(paste0(out.dir,label,'.models.Rds'))
  
  }else{
    models <- paleocar_models(chronologies = chronologies,
                              predictands = predictands,
                              calibration.years = calibration.years,
                              verbose = verbose,
                              ...)
    readr::write_rds(models,
                     path = paste0(out.dir,label,".models.Rds"),
                     compress = "gz")
  }
  
  if(return.predictions){
    if(verbose) cat("\nGenerating prediction\n")
    if(!force.redo & file.exists(paste0(out.dir,label,".prediction.Rds"))){
      recon <- readr::read_rds(paste0(out.dir,label,".prediction.Rds"))
      
    }else{
      recon <- predict_paleocar_models(models = models,
                                       ...)
      readr::write_rds(recon,
                       path = paste0(out.dir,label,".prediction.Rds"),
                       compress = "gz")
      
    }
  }else{
    recon = NULL
  }
  
  if(return.uncertainty){
    if(verbose) cat("\nGenerating uncertainty predictions")
    if(!force.redo & file.exists(paste0(out.dir,label,".uncertainty.Rds"))){
      uncertainty <- readr::read_rds(paste0(out.dir,label,".uncertainty.Rds"))
      
    }else{
      uncertainty <- uncertainty_paleocar_models(models = models,
                                       ...)
      readr::write_rds(uncertainty,
                       path = paste0(out.dir,label,".uncertainty.Rds"),
                       compress = "gz")
      
    }
  }else{
    uncertainty = NULL
  }
  
  if(verbose) message("\nThe entire reconstruction took ", round(difftime(Sys.time(),t,units='mins'), digits = 2)," minutes")
  
  if(return.objects){
    return(
      list(models = models,
           predictions = recon,
           uncertainty = uncertainty)
    ) 
  }else{
    return(
      message("Reconstruction finished")
    )
  }
  
}
