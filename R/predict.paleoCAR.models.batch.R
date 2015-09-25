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
#' @param prediction.years The set of years over which to generate reconstruction rasters. Optional.
#' @return A RasterBrick containing the predictions for each year.
predict.paleocar.models.batch <- function(models, meanVarMatch = TRUE, prediction.years=NULL){
  if(is.null(prediction.years)) prediction.years <- as.numeric(rownames(models[['reconstruction.matrix']]))
  
  if(!all(prediction.years %in% as.numeric(row.names(models[['reconstruction.matrix']])))){
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    prediction.years <- prediction.years[prediction.years %in% as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  models[['reconstruction.matrix']] <- models[['reconstruction.matrix']][as.numeric(rownames(models[['reconstruction.matrix']])) %in% prediction.years,]
  
  # models[['models']][,endYear := c(year[-1]-1,2000),by=cell]
  
  newx <- data.table(cbind(1,models[['reconstruction.matrix']]))
  setnames(newx,c("Intercept",colnames(models[['reconstruction.matrix']])))
  rownames(newx) <- prediction.years
  
  predictands <- t(raster::as.matrix(models[["predictands"]]))
#   predictands.mean <- colMeans(predictands)
#   predictands.sd <- colSds(predictands)
#   predictands <- data.table(cell=1:ncol(predictands),mean=predictands.mean,sd=predictands.sd)
#   setkey(predictands,cell)
  
  newx.calib <- data.table(cbind(1,models[['predictor.matrix']]))
  setnames(newx.calib,c("Intercept",colnames(models[['predictor.matrix']])))
  rownames(newx.calib) <- rownames(models[['predictor.matrix']])
  newx.calib <- as.matrix(newx.calib)
  
  predictions <- sapply(unique(models[['models']][['cell']]),function(this.cell){
    coefficients <- models[['models']][cell==this.cell]
    coefficients <- coefficients[base::rep(1:nrow(coefficients),times=diff(c(coefficients[['year']],tail(prediction.years,1)+1))),]
    coefficients <- coefficients[,c("Intercept",colnames(models[['reconstruction.matrix']])),with=F]
    this.predictions <- rowSums(coefficients*newx, na.rm=T)
    
    if(meanVarMatch){
      calibration.years <- as.numeric(gsub("X","",names(models[["predictands"]])))
      
      calibration.coefficients <- models[['models']][cell==this.cell]
      calibration.coefficients <- calibration.coefficients[,c("Intercept",colnames(models[['reconstruction.matrix']])),with=F]
      
      calibration.predictors <- as.matrix(newx[calibration.years])[rep(1:length(calibration.years),nrow(calibration.coefficients)),]
      
      calibration.coefficients <- as.matrix(calibration.coefficients)[rep(1:nrow(calibration.coefficients),each=nlayers(models[["predictands"]])),]
      
      calibration.predictions <- rowSums(calibration.coefficients*calibration.predictors, na.rm=T)
      calibration.predictions <- do.call(cbind,split(calibration.predictions,rep(1:nrow(models[['models']][cell==this.cell]),each=length(calibration.years))))
      
      calibration.predictands <- predictands[,this.cell]
      
      scalar <- sd(calibration.predictands)/sd(calibration.predictions)
      
      predictions.meanvar <- calibration.predictions*scalar + (mean(calibration.predictands)-(scalar*mean(calibration.predictions)))
      return(this.predictions.scaled.meaned)
    }
    
    return(this.predictions)
  })
  colnames(predictions) <- unique(models[['models']][['cell']])
  
  

  out <- rowSums(test*newx, na.rm=T)
  plot(as.numeric(recon[1]))
  points(out, col='red', pch=19)
  points(y=predictands[,1],x=1924:1983,col='blue',pch=19)
  
  #   predictions <- sapply(prediction.years,function(this.year){
  #     this.newx <- newx[this.year,]
  #     coefficients <- models[['models']][year<=this.year & endYear>=this.year,c("Intercept",colnames(models[['reconstruction.matrix']])),with=F]    
  #     return(rowSums(coefficients*this.newx[rep(1,nrow(coefficients))], na.rm=T))
  #   })
  #   colnames(predictions) <- prediction.years
  
#   if(meanVarMatch){
#     predictands <- t(raster::as.matrix(models[["predictands"]]))
#     predictands.mean <- colMeans(predictands)
#     predictands.sd <- colSds(predictands)
#     predictands <- data.table(cell=1:ncol(predictands),mean=predictands.mean,sd=predictands.sd)
#     setkey(predictands,cell)
#     
#     newx.calib <- data.table(cbind(1,models[['predictor.matrix']]))
#     setnames(newx.calib,c("Intercept",colnames(models[['predictor.matrix']])))
#     rownames(newx.calib) <- rownames(models[['predictor.matrix']])
#     newx.calib <- as.matrix(newx.calib)
#     
#     calibration.predictions <- cbind(models[['models']][,.(cell,year)],do.call(rbind,lapply(unique(models[['models']][['cell']]),function(this.cell){
#       coefficients <- models[['models']][cell==this.cell]
#       newx <- newx.calib[rep(1:nrow(newx.calib),nrow(coefficients)),]
#       coefficients <- coefficients[base::rep(1:nrow(coefficients),each=nrow(newx.calib)),]
#       coefficients <- as.matrix(coefficients[,c("Intercept",colnames(models[['reconstruction.matrix']])),with=F])
#       out <- do.call(cbind,split(rowSums(coefficients*newx, na.rm=T),rep(1:nrow(newx.calib),each=nrow(models[['models']][cell==this.cell]))))
#       out <- data.table(mean=rowMeans(out), sd=rowSds(out))
#       
#       return(out)
#     })))
#     
#     #     # Mean and variance matching over calibration period    
#     #     calibration.predictions <- lapply(1:nrow(newx.calib),function(this.row){
#     #       this.newx <- newx.calib[this.row,]
#     #       return(rowSums(coefficients*this.newx[rep(1,nrow(coefficients))], na.rm=T))
#     #     })
#     
#     #     calibration.predictions <- do.call(cbind,calibration.predictions)
#     #     calibration.predictions.mean <- rowMeans(calibration.predictions)
#     #     calibration.predictions.sd <- rowSds(calibration.predictions)
#     #     calibration.predictions <- models[['models']][,.(cell,year,endYear)]
#     #     calibration.predictions[,mean:=calibration.predictions.mean]
#     #     calibration.predictions[,sd:=calibration.predictions.sd]
#     
#     setkey(calibration.predictions,cell)
#     all <- predictands[calibration.predictions]
#     all[,scalar:=sd/i.sd]
#     
#     scalar.matrix <- sapply(unique(all[['cell']]),function(this.cell){
#       
#       coefficients <- models[['models']][cell==this.cell]
#       coefficients <- coefficients[base::rep(1:nrow(coefficients),times=diff(c(coefficients[['year']],tail(prediction.years,1)+1))),]
#       coefficients <- coefficients[,c("Intercept",colnames(models[['reconstruction.matrix']])),with=F]    
#       return(rowSums(coefficients*newx, na.rm=T))
#     })
#     
#     
#     predictions <- sapply(colnames(predictions), function(this.cell){
#       coefficients <- models[['models']][cell==this.cell]
#       coefficients <- coefficients[base::rep(1:nrow(coefficients),times=diff(c(coefficients[['year']],tail(prediction.years,1)+1))),]
#       coefficients <- coefficients[,c("Intercept",colnames(models[['reconstruction.matrix']])),with=F]    
#       return(rowSums(coefficients*newx, na.rm=T))
#     })
#     
#     predictions <- sapply(colnames(predictions), function(this.year){
#       this.predictions <- predictions[,this.year]
#       this.all <- all[year<=as.numeric(this.year) & endYear>=as.numeric(this.year)]
#       this.scalar <- this.all[,scalar]
#       this.predictions.scaled <- this.predictions*this.scalar
#       this.predictand.mean <- this.all[,mean]
#       this.calibration.predictions.mean <- this.all[,i.mean]
#       this.predictions.scaled.meaned <- this.predictions.scaled + (this.predictand.mean-(this.scalar*this.calibration.predictions.mean))
#       return(this.predictions.scaled.meaned)
#     })
#   }
  
  if(class(models[["predictands"]]) %in% c("RasterBrick","RasterStack")){
    predictions <- setValues(models[["predictands"]],t(predictions))
  }
  
  return(predictions)
}
