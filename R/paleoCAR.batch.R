#' Fit PaleoCAR models to a RasterBrick of predictands, and generate reconstruction, error, and model size RasterBricks
#'
#' This is a wrapper function to the four primary methods in the PaleoCAR package, specifically for
#' processing predictands held in a RasterBrick. Additionally, it allows for the specification of
#' a \code{floor} and/or \code{ceiling} corrections for a reconstruction (e.g., of precipitation).
#' It returns a list containing the models, predictions, cross-validated errors, and model sizes.
#'
#' @param chronologies A matrix of tree ring chronologies, indexed annually.
#' Each chronology is a column. The first column must be labeled "YEAR" and is the calendar year.
#' @param predictands A RasterBrick or RasterStack of the numeric predictand (response) variable.
#' @param calibration.years An integer vector of years corresponding to the layers in the \code{predictands} brick.
#' @param prediction.years An optional integer vector of years for the reconstruction.
#' If missing, defaults to the total years present in \code{chronologies}.
#' @param label A character label for the reconstruction, for saving.
#' @param out.dir The directory to which output is to be saved.
#' @param min.width integer, indicating the minimum number of tree-ring samples allowed for that year of a chronology to be valid.
#' @param meanVar A character string indicating the type of mean-variance matching to perform: either "none" (default), "calibration", or "chained".
#' @param chained.meanVar Logical, indicating whether to chain mean-variance matching if performed. See \code{\link{predict.paleocar.models.batch}}.
#' @param floor Numeric, an optional lower bound for reconstructed values, such as \code{0} for precipitation reconstructions.
#' @param ceiling Numeric, an optional upper bound for reconstructed values.
#' @param asInt Logical, should reconstructed values be rounded to integers for saving?
#' @param force.redo Logical, should all computations be re-computed?
#' @param generate.reconstruction Logical, should the reconstruction be generated (as opposed to only the models)?
#' @param return.objects Logical, should objects be returned (as opposed to only writing them to disk)?
#' @param verbose Logical, display status messages during run.
#' @return A named list containing
#' \itemize{
#'   \item{\code{models}  The PaleoCAR models, as computed by \code{\link{paleoCAR.models.batch}}.}
#'   \item{\code{recon}  The PaleoCAR reconstruction, as computed by \code{\link{predict.paleocar.models.batch}}.}
#'   \item{\code{errors}  The PaleoCAR reconstruction average LOOCV error, as computed by \code{\link{errors.paleocar.models.batch}}.}
#'   \item{\code{sizes}  The PaleoCAR model sizes, as computed by \code{\link{size.paleocar.models.batch}}.}
#' }
paleoCAR.batch <- function(chronologies, predictands, calibration.years, prediction.years=NULL, label, out.dir="./OUTPUT/", min.width=NULL, meanVar = "none", floor=NULL, ceiling=NULL, asInt=F, force.redo=F, verbose=F, generate.reconstruction=T, return.objects=T){
  t <- Sys.time()
  if(verbose) cat("\nCalculating all models")
  models <- paleoCAR.models.batch(chronologies=chronologies, predictands=predictands, calibration.years=calibration.years, prediction.years=prediction.years, label=label, out.dir=out.dir, min.width=min.width, force.redo=force.redo, verbose=verbose)
  
  if(generate.reconstruction){
    if(verbose) cat("\nGenerating reconstruction")
    if(!force.redo & file.exists(paste(out.dir,label,".recon.tif",sep=''))){
      
      recon <- raster::brick(paste(out.dir,label,".recon.tif",sep=''))
      
    }else{
      recon <- predict.paleocar.models.batch(models=models, meanVar=meanVar, prediction.years=prediction.years)
      
      if(!is.null(floor)){
        recon <- raster::calc(recon,function(x){x[x<floor] <- floor; return(x)})
      }
      if(!is.null(ceiling)){
        recon <- raster::calc(recon,function(x){x[x>ceiling] <- ceiling; return(x)})
      }
      if(asInt){
        recon <- raster::calc(recon,function(x){round(x,digits=0)})
        type <- "INT2S"
      }else{type="FLT4S"}
      raster::writeRaster(recon,paste(out.dir,label,".recon.tif",sep=''), datatype=type, options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
      
      #     if(asInt){
      #       recon$errors <- calc(recon$errors,function(x){round(x,digits=0)})
      #       type="INT2U"
      #     }else{type="FLT4S"}
      #     raster::writeRaster(recon$errors,paste(out.dir,label,".errors.tif",sep=''), datatype=type, options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
      #     
      #     raster::writeRaster(recon$sizes,paste(out.dir,label,".size.tif",sep=''), datatype="INT2U", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
      #     
      # saveRDS(recon,paste(out.dir,label,".recon.Rds",sep=''),compress='xz')
    }
  }
  if(verbose) cat("\nThe entire reconstruction took", difftime(Sys.time(),t,units='mins'),"minutes")
  
  if(return.objects){
    if(generate.reconstruction){
      return(list(models=models,recon=recon)) 
    }else{
      return(models)
    }
  }else{
    rm(models,recon); gc(); gc()
    return(NULL)
  }
}

