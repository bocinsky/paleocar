globalVariables(c("end.year","coefs"))
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
#' @param models A PaleoCAR batch model, as returned from \code{\link{paleocar_models_batch}}.
#' @param meanVar A character string indicating the type of mean-variance matching to perform: either "none" (default), "calibration", or "chained".
#' @param prediction.years The set of years over which to generate reconstruction rasters. Optional.
#' @return A RasterBrick containing the predictions for each year.
#' @export
predict_paleocar_models_batch <- function(models, meanVar = "none", prediction.years=NULL){
  if(is.null(prediction.years)) prediction.years <- as.numeric(rownames(models[['reconstruction.matrix']]))
  
  if(!all(prediction.years %in% as.numeric(row.names(models[['reconstruction.matrix']])))){
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    prediction.years <- prediction.years[prediction.years %in% as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  models[['reconstruction.matrix']] <- models[['reconstruction.matrix']][as.numeric(rownames(models[['reconstruction.matrix']])) %in% prediction.years,]
  
  newx <- data.table::data.table(cbind(1,models[['reconstruction.matrix']]))
  setnames(newx,c("Intercept",colnames(models[['reconstruction.matrix']])))
  rownames(newx) <- prediction.years
  newx.present <- !is.na(as.matrix(newx)[,-1])
  
  predictands.matrix <- t(raster::as.matrix(models[["predictands"]]))
  
  newx.calib <- as.matrix(cbind(data.table(Intercept=1),models[['predictor.matrix']]))
  rownames(newx.calib) <- rownames(models[['predictor.matrix']])
  calib.years <- as.numeric(row.names(newx.calib))
  
  predictions <- sapply(unique(models[['models']][['cell']]),function(this.cell){
    # cat(this.cell,"\n")
    this.models <- models[['models']][cell==this.cell & year %in% prediction.years]
    this.models[,end.year:=c(year[-1]-1,tail(prediction.years,1))]
    
    coefficients <- data.table::rbindlist(lapply(this.models$coefs,function(x){data.table::data.table(matrix(data=x,ncol=length(x),byrow=T,dimnames=list(NA,names(x))))}),fill=T)
    this.newx <- newx[,names(coefficients),with=F]
    
    model.rows <- rep(1:nrow(coefficients),times=diff(c(this.models$year,tail(prediction.years,1)+1)))
    
    this.predictions <- rowSums(coefficients[model.rows]*this.newx, na.rm=T)
    
    if(meanVar == "chained"){
      
      spans <- lapply(1:length(this.models$year),function(i){this.models$year[i]:(c(this.models$year,tail(prediction.years,1)+1)[i+1]-1)})
      pivot <- which(sapply(spans,function(span){all(calib.years %in% span)}))
      
      if(pivot>1){
        available.below <- apply(this.models[1:(pivot-1),list(end.year,coefs)],1,function(d){
          out <- prediction.years[which(rowSums(!is.na(newx[,names(d$coefs),with=F]))==length(names(d$coefs)))]
          out <- out[out>d[1]]
          out.rle <- rle(diff(out))
          out.rle.start <- which(cumsum(out.rle$lengths) == cumsum(out.rle$lengths)[head(which(out.rle$lengths>=length(calib.years) & out.rle$values==1),1)])
          if(length(out.rle.start) > 0){
            out <- out[out.rle.start:(out.rle.start+length(calib.years)-1)]
          }
          return(out)
        })
      }else available.below <- NULL
      
      if(length(spans)>pivot){
        available.above <- apply(this.models[(pivot+1):nrow(this.models),list(year,coefs)],1,function(d){
          out <- which(rowSums(!is.na(newx[,names(d$coefs),with=F]))==length(names(d$coefs)))
          out <- out[out<d[1]]
          out.rle <- rle(diff(out))
          out.rle.end <- which(cumsum(out.rle$lengths)==cumsum(out.rle$lengths)[tail(which(out.rle$lengths>=length(calib.years) & out.rle$values==1),1)])
          out <- rev(rev(out)[out.rle.end:(out.rle.end+length(calib.years)-1)])
          return(out)
        })
      }else available.above <- NULL
      
      
      availables <- list(available.below,as.matrix(calib.years),available.above)
      availables <- availables[!sapply(availables,is.null)]
      match.periods <- do.call(cbind,availables)
      
      match.newx <- apply(match.periods,2,function(x){
        this.newx[x]
      })
      
      match.predictions <- matrix(rowSums(coefficients[rep(1:nrow(coefficients),each=length(calib.years))]*this.newx[match(as.vector(match.periods),prediction.years)], na.rm=T),nrow=length(calib.years))
      
      scale.center <- mean_var_match(calib.vector=predictands.matrix[,this.cell],out.calib.vector=match.predictions[,pivot], out.vector=this.predictions[match(spans[[pivot]],prediction.years)])
      names(scale.center) <- spans[[pivot]]
      
      if(pivot>1){
        for(i in (pivot-1):1){
          new.names <- c(as.character(spans[[i]]),names(scale.center))
          scale.center <- c(mean_var_match(calib.vector=scale.center[as.character(match.periods[,i])],out.calib.vector=match.predictions[,i], out.vector=this.predictions[match(spans[[i]],prediction.years)]),scale.center)
          names(scale.center) <- new.names
        }
      }
      
      if(length(spans)>pivot){
        for(i in (pivot+1):length(spans)){
          new.names <- c(names(scale.center),as.character(spans[[i]]))
          scale.center <- c(scale.center,mean_var_match(calib.vector=scale.center[as.character(match.periods[,i])],out.calib.vector=match.predictions[,i], out.vector=this.predictions[match(spans[[i]],prediction.years)]))
          names(scale.center) <- new.names
        }
      }
      
      this.predictions <- scale.center
      
    }else if(meanVar == "calibration"){
      this.newx.calib <- newx.calib[,names(coefficients)]
      calibration.years <- rownames(this.newx.calib)
      
      calibration.coefficients <- coefficients[,colnames(this.newx.calib),with=F]
      
      calibration.predictors <- as.matrix(this.newx[which(prediction.years %in% as.numeric(calibration.years))])[rep(1:length(calibration.years),nrow(calibration.coefficients)),]
      
      calibration.coefficients <- as.matrix(calibration.coefficients)[rep(1:nrow(calibration.coefficients),each=length(calibration.years)),]
      
      calibration.predictions <- rowSums(calibration.coefficients*calibration.predictors, na.rm=T)
      calibration.predictions <- do.call(rbind,split(calibration.predictions,rep(1:nrow(coefficients),each=length(calibration.years))))
      
      calibration.predictands <- predictands.matrix[,this.cell]
      
      scalars <- sd(calibration.predictands)/rowSds(calibration.predictions)
      transforms <- mean(calibration.predictands)-(scalars*rowMeans(calibration.predictions))
      
      this.predictions <- (this.predictions*scalars[model.rows]) + transforms[model.rows]
      
    }
    return(this.predictions)
  })


  non_na_cells <- which(!is.na(models[["predictands"]][[1]][]))
  new.matrix <- matrix(NA,nrow=ncell(models[["predictands"]]),ncol=length(prediction.years))
  new.matrix[non_na_cells,] <- t(predictions)
  rm(predictions);gc();gc()
  predictions <- raster::setValues(models[["predictands"]],new.matrix)

  return(predictions)
  
}


