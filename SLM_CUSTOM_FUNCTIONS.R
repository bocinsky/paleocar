source('./UTILITY_FUNCTIONS.R')

slm.models.custom = function(Xtrain, Ytrain, predlist, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
{
  p = dim(Xtrain)[2]
  m = length(predlist)
  numpred = sapply(predlist, length)
  R2 = numeric(m)
  coeff = matrix(0, nrow = m, ncol = p + 1)
  modelnames = names(predlist)
  if (is.null(modelnames)) 
    modelnames = paste("SIZE.", numpred, sep = "")
  xnames = colnames(Xtrain)
  if (is.null(xnames)) 
    xnames = paste("X", 1:p, sep = "")
  rownames(coeff) = modelnames
  colnames(coeff) = c("(Intercept)", xnames)
  names(R2) = modelnames
  names(numpred) = modelnames
  for (i in 1:m) {
    if (verbose) 
      cat("Determine regression coefficients for", modelnames[i], 
          "model\n")
    idx = predlist[[i]]
    fit = slm(Xtrain[, idx, drop = FALSE], Ytrain, diagonal = diagonal, 
              lambda = lambda, lambda.var = lambda.var, verbose = verbose)
    coeff[i, 1 + idx] = fit$coefficients[-1]
    coeff[i, 1] = fit$coefficients[1]
    R2[i] = fit$R2
  }
  return(list(coefficients = coeff, numpred = numpred, R2 = R2))
}


slm.models.custom.brick.unique <- function(X, Y.brick, predlist, name, lambdas=NULL, lambda.vars=NULL, regularize=T, max.series=nrow(X)-5, out.dir="../OUTPUT/", cellNums=1:ncell(Y.brick), select.stat="RMSE", force.redo=F){
  if(!force.redo & file.exists(file=paste(out.dir,name,'_selected_lms_unique.rds',sep=''))){
    lms.out.unique <- readRDS(file=paste(out.dir,name,'_selected_lms_unique.rds',sep=''))
    return(lms.out.unique)
  }
  
  # Limit number of chronologies entering into liear model at the number of observations (n-1)
  if(max.series >= ncol(X)){
    max.series <- nrow(X)-5
  }
  
  hinge.year <- max(as.numeric(gsub("X","",names(Y.brick))))
  
  ## estimate correlation shrinkage intensity
  if(is.null(lambdas)){
    lambdas <-  mclapply(cellNums,function(i){estimate.lambda(cbind(as.numeric(Y.brick[i]),X), verbose=F)}, mc.cores=detectCores())
    names(lambdas) <- cellNums
  }else if(length(lambdas) != length(cellNums)){
    stop("Number of specified lambda values must equal number of cells")
  }
  
  ## estimate correlation shrinkage intensity
  if(is.null(lambda.vars)){
    lambda.vars <-  mclapply(cellNums,function(i){estimate.lambda.var(cbind(as.numeric(Y.brick[i]),X), verbose=F)}, mc.cores=detectCores())
    names(lambda.vars) <- cellNums
  }else if(length(lambda.vars) != length(cellNums)){
    stop("Number of specified lambda.var values must equal number of cells")
  }
  
  ## Calculate CAR scores across all cells using shrinkage estimates of the correlation matrices
  # 0.385 seconds/16
  carscores <- mclapply(cellNums,function(i){carscore(Xtrain=X, Ytrain=as.numeric(Y.brick[i]), lambda=lambdas[[as.character(i)]], verbose=FALSE)}, mc.cores=detectCores())
  saveRDS(carscores,file=paste(out.dir,name,'_carscores.rds',sep=''), compress="xz")
  
  # Output the shrinkage parameter if calculating on the entire brick
  if(length(cellNums)==ncell(Y.brick)){  
    lambdas.rast <- setValues(raster(Y.brick),unlist(lambdas))
    writeRaster(lambdas.rast,paste(out.dir,name,'_lambdas.tif',sep=''),datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
    
    lambda.vars.rast <- setValues(raster(Y.brick),unlist(lambda.vars))
    writeRaster(lambda.vars.rast,paste(out.dir,name,'_lambda_vars.tif',sep=''),datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
  }
  
  # Calculate the squared CAR scores, and sort in decreasing order
  # 0.026 seconds/16
  car2 <- mclapply(carscores,function(x){sort(x^2, decreasing=T)}, mc.cores=detectCores())
  
  # Generate cell-wise lists of predictor varibles sorted by squared CAR scores
  predlist.orders <- mclapply(car2,function(x){unlist(lapply(names(x), function(i,...){which(i == names(X))}))}, mc.cores=detectCores()) # 0.033 seconds/16
  predlists <- mclapply(predlist.orders, function(x,...){
    lapply(predlist,function(i){
      out <- x[x %in% i][1:max.series]
      return(out <- out[!is.na(out)])
      #       return(make.predlist(out,numpred=1:length(out)))
    })
  }, mc.cores=detectCores()) # 0.282 seconds/16
  names(predlists) <- cellNums
  
  lms.out.unique <- mclapply(cellNums, function(cell){
    cell.predlists <- predlists[[as.character(cell)]]
    cell.predlists <- cell.predlists[order(as.numeric(names(cell.predlists)))]
    #     cell.predlists <- cell.predlists[order(as.numeric(names(cell.predlists)), decreasing=T)]
    #     cell.predlists <- cell.predlists[order(unlist(lapply(cell.predlists,length)), decreasing=T)]
    #     
    numpred.max <- max.series
    cell.predlists.down <- cell.predlists[as.numeric(names(cell.predlists))<hinge.year]
    cell.predlists.down <- cell.predlists.down[order(as.numeric(names(cell.predlists.down)), decreasing=T)]
    cell.predlists.up <- cell.predlists[as.numeric(names(cell.predlists))>hinge.year]
    
    calcYear <- function(year,...){
      #       year.count <<- year.count+1
      if(!is.null(selected.preds) & all(selected.preds %in% year)){
        return(NA)
      }
      predlist.partial <- make.predlist(year,numpred=1:length(year))
      slms <- cv.select.slm.models(X=X, Y=as.numeric(Y.brick[cell]), predlist=predlist.partial, lambda=lambdas[[as.character(cell)]], lambda.var=lambda.vars[[as.character(cell)]], regularize=regularize)
      errors <- Reduce('+',slms$CV.ERRORS)/length(slms$CV.ERRORS)
      min.stat <- which(errors[,select.stat] == min(errors[,select.stat]))
      #       min.stat <- getFirstMin(errors[,select.stat])
      #       cat(min.stat,'\n')
      slms$MODELS <- slms$MODELS[min.stat,]
      slms$CV.ERRORS <- do.call(rbind,lapply(slms$CV.ERRORS,function(x){x[min.stat,]}))
      selected.preds <<- year[1:min.stat]
      rownames(slms$CV.ERRORS) <- 1:nrow(slms$CV.ERRORS)
      slms$SELECTED <- selected.preds
      #       slms$CARS <- carscore(Xtrain=X[,selected.preds,drop=F],Ytrain=as.numeric(Y.brick[cell]), lambda=lambdas[[as.character(cell)]], verbose=F)
      return(slms)
    }
    
    selected.preds <- NULL
    #     year.count <- 0
    down.out <- lapply(cell.predlists.down, calcYear)
    down.out <- down.out[!is.na(down.out)]
    true.names <- unlist(lapply(1:length(names(down.out)),function(i){names(cell.predlists)[which(names(cell.predlists) == names(down.out)[i+1])+1]}))
    names(down.out) <- c(true.names,1)
    selected.preds <- down.out[[1]]$SELECTED
    up.out <- lapply(cell.predlists.up, calcYear)
    up.out <- up.out[!is.na(up.out)]
    
    out <- c(down.out,up.out)
    out <- out[!is.na(out)]
    out <- out[order(as.numeric(names(out)))]
    
    coefs <- do.call(rbind,lapply(out,function(x){x$MODELS}))
    coefs <- coefs[order(as.numeric(rownames(coefs))),]
    duplicates <- sequential.duplicated(as.data.frame(coefs), rows=T)
    coefs <- coefs[!duplicates,]
    beta.weights <- apply(X,2,sd)/sd(as.numeric(Y.brick[cell]))
    betas <- t(apply(coefs[,-1],1,function(x){x*beta.weights}))
    
    choices <- apply(coefs,1,function(x){colnames(coefs)[which(x != 0)][-1]})
    if(regularize){
      cars <- mclapply(choices,function(choice,...){carscore(Xtrain=X[,choice,drop=F], Ytrain=as.numeric(Y.brick[cell]), lambda=lambdas[[as.character(cell)]], verbose=F)})
    }else{
      cars <- mclapply(choices,function(choice,...){carscore(Xtrain=X[,choice,drop=F], Ytrain=as.numeric(Y.brick[cell]), verbose=F)})
    }
    errors <- lapply(out,function(x){x$CV.ERRORS})
    errors <- errors[names(choices)]
    errors <- do.call(rbind,lapply(errors, colMeans))
    
    #     errors.means <- do.call(rbind,lapply(errors,colMeans))
    #     aics <- unlist(lapply(out,function(x){x$AIC}))
    #     sizes <- as.numeric(gsub(".*.SIZE.","",names(aics)))
    #     names(sizes) <- as.numeric(gsub(".SIZE.*","",names(aics)))
    #     sizes <- sizes[order(as.numeric(names(sizes)))]
    return(list(coefficients=coefs, betas=betas, cars=cars, errors=errors))
  }, mc.cores=detectCores())
  
  names(lms.out.unique) <- cellNums
  saveRDS(lms.out.unique,file=paste(out.dir,name,'_selected_lms_unique.rds',sep=''), compress="xz")
  
  return(lms.out.unique)
}

# A function to calculate a group of shrinkage linear models, and 
# calculate the AICs and cross-validated statistics of all models
# @returns MODELS A matrix of linear model coefficients
# @returns LAMBDA The correlation shrinkage intensity of the models
# @returns ALL.CAR A numeric vector of raw CAR scores for all Xtrain variables
# @returns PART.CAR A matrix of per-model CAR scores
# @returns AIC A numeric vector of AIC scores
# @returns CORRECTED A boolean term describing whether or not the AIC is corrected for small sample sizes
cv.select.slm.models <- function(X, Y, lambda=NULL, lambda.var=NULL, regularize=T, diagonal=F, predlist=NULL, meanVar=T, K=3, R=1, type="consecutive", fold.list=NULL){
  if(is.null(lambda)){
    lambda <- estimate.lambda(cbind(Y,X), verbose=F)
  }
  
  if(is.null(lambda.var)){
    lambda.var <- estimate.lambda.var(cbind(Y,X), verbose=F)
  }
  
  if(is.null(predlist)){
    car <- carscore(Xtrain=X, Ytrain=Y, lambda=lambda, diagonal=diagonal, verbose=F)
    predlist <- make.predlist(order((car)^2, decreasing=T), numpred = 1:(nrow(X)-5), name="CAR")
  }
  
  # Estimate coefficients for all linear models
  if(regularize){
    out.slms <- slm.models(Xtrain=X, Ytrain=Y, predlist=predlist, lambda=lambda, lambda.var=lambda.var, verbose=F)$coefficients
  }else{
    out.slms <- slm.models(Xtrain=X, Ytrain=Y, predlist=predlist, verbose=F)$coefficients
  }
  
  if(is.null(fold.list)){
    if(K==1){
      folds.out <- vector('list',1)
      folds.out[[1]]$TRAIN <- 1:nrow(X)
      folds.out[[1]]$VALIDATE <- 1:nrow(X)
    }else{
      folds <- cvFolds(nrow(X), K=K, R=R, type=type)
      folds.out <- unlist(lapply(1:folds$K,FUN=function(x,...){ getFolds(folds=folds,i=x) }), recursive=F)
      if(type!="consecutive" & length(folds.out)<=R & K!=nrow(X)){
        folds.out <- unlist(folds.out, recursive=F)
      }
    }
  }else{
    folds.out <- fold.list
  }
  
  errors <- lapply(folds.out,function(fold,...){
    Xtrain <- X[fold$TRAIN,]
    Ytrain <- Y[fold$TRAIN]
    Xvalidation <- X[fold$VALIDATE,]
    Yvalidation <- Y[fold$VALIDATE]
    
    errors <- cv.slm.models(Xtrain=Xtrain, Ytrain=Ytrain, Xvalidation=Xvalidation, Yvalidation=Yvalidation, car.predlist=predlist, lambda=lambda, lambda.var=lambda.var, regularize=regularize)
  })
  
  errors <- lapply(errors,"[[","ERROR")
  
  return(list(MODELS=out.slms, CV.ERRORS=errors))
}

## This function takes output from slm.models.custom.brick.redux()
## and creates a set of temporal models that utilize the temporal
## union of all selected t-r chronologies. The outcome is that
## the resulting cell-wise sets of linear models will each use the same
## chronologies, though the coefficients will differ from one another. This
## serves to remove spatial artifacts from the reconstructions that are created
## by inconsistant correlation structures between the chronologies.
slm.models.custom.brick.union <- function(models, X, Y.brick, name, lambdas=NULL, lambda.vars=NULL, regularize=T, out.dir="../OUTPUT/", force.redo=F){
  if(!force.redo & file.exists(file=paste(out.dir,name,'_selected_lms_union.rds',sep=''))){
    lms.out.union <- readRDS(file=paste(out.dir,name,'_selected_lms_union.rds',sep=''))
    return(lms.out.union)
  }
  
  cellNums <- as.numeric(names(models))
  
  if(regularize){
    ## estimate correlation shrinkage intensity
    if(is.null(lambdas)){
      lambdas <-  mclapply(cellNums,function(i){estimate.lambda(cbind(as.numeric(Y.brick[i]),X), verbose=F)}, mc.cores=detectCores())
      names(lambdas) <- cellNums
    }else if(length(lambdas) != length(cellNums)){
      stop("Number of specified lambda values must equal number of cells")
    }
    
    ## estimate correlation shrinkage intensity
    if(is.null(lambda.vars)){
      lambda.vars <-  mclapply(cellNums,function(i){estimate.lambda.var(cbind(as.numeric(Y.brick[i]),X), verbose=F)}, mc.cores=detectCores())
      names(lambda.vars) <- cellNums
    }else if(length(lambda.vars) != length(cellNums)){
      stop("Number of specified lambda.var values must equal number of cells")
    }
  }
  
  models.coefficients <- lapply(models,'[[','coefficients')
  
  union.breaks <- unique(unlist(lapply(models.coefficients,function(cell){as.numeric(rownames(cell))})))
  union.breaks <- union.breaks[order(union.breaks)]
  
  selected <- mclapply(models.coefficients, function(cell){
    coefficients.expand <- matrix.expand(cell,union.breaks)
    coefficients.expand[coefficients.expand != 0] <- 1
    return(coefficients.expand)
  }, mc.cores=detectCores())
  
  # convert to 3D matrix
  selected <- abind(selected, along=3)
  
  # collapse along z-axis
  selected <- apply(selected, c(1,2), max)[,-1]
  
  # convert to predlist
  selected.predlist <- apply(selected,1,function(year){
    return(which(colnames(X) %in% colnames(selected)[year==1]))
  })
  
  # Calculate new series of linear models for each cell
  lms.out.union <- mclapply(as.numeric(names(models)), function(cell){
    
    out <- cv.select.slm.models(X=X, Y=as.numeric(Y.brick[cell]), predlist=selected.predlist, lambda=lambdas[[as.character(cell)]], lambda.var=lambda.vars[[as.character(cell)]], regularize=regularize)
    
    coefs <- out$MODELS
    coefs <- coefs[order(as.numeric(rownames(coefs))),]
    
    beta.weights <- apply(X,2,sd)/sd(as.numeric(Y.brick[cell]))
    betas <- t(apply(coefs[,-1],1,function(x){x*beta.weights}))
    
    choices <- apply(coefs,1,function(x){colnames(coefs)[which(x != 0)][-1]})
    
    if(regularize){
      cars <- mclapply(choices,function(choice,...){carscore(Xtrain=X[,choice,drop=F], Ytrain=as.numeric(Y.brick[cell]), lambda=lambdas[[as.character(cell)]], verbose=F)})
    }else{
      cars <- mclapply(choices,function(choice,...){carscore(Xtrain=X[,choice,drop=F], Ytrain=as.numeric(Y.brick[cell]), verbose=F)})
    }
    
    errors <- out$CV.ERRORS
    errors <- Reduce(x=errors, f='+')/3
    
    return(list(coefficients=coefs, betas=betas, cars=cars, errors=errors))
  }, mc.cores=detectCores())
  
  names(lms.out.union) <- cellNums
  saveRDS(lms.out.union,file=paste(out.dir,name,'_selected_lms_union.rds',sep=''), compress="xz")
  
  return(lms.out.union)
  
}

recon <- function(models, Ytrain.brick, recon.series, training.series, name, out.dir="../OUTPUT/", force.redo=F){
  if(length(models)==ncell(Ytrain.brick)){
    raster.out=T
  }else{
    raster.out=F
  }
  
  
  if(!force.redo){
    if(raster.out & file.exists(paste(out.dir,name,'_recons.tif',sep=''))){
      recons.brick <- brick(paste(out.dir,name,'_recons.tif',sep=''))
      return(recons.brick)
    }else if(!raster.out & file.exists(paste(out.dir,name,'_recons.rds',sep=''))){
      recons <- readRDS(file=paste(out.dir,name,'_recons.rds',sep=''))
      return(recons)
    }
  } 
  
  coefficients <- lapply(models,'[[','coefficients')
  
  Ytrain.mean <- calc(Ytrain.brick,mean)
  Ytrain.sd <- calc(Ytrain.brick, sd)
  # Pure reconstruction
  recons <- mclapply(1:length(coefficients), function(i){
    breaks <- as.numeric(rownames(coefficients[[i]]))
    coefficients.expand <- matrix.expand(coefficients[[i]],1:2000)
    predictions <- rowSums(coefficients.expand * cbind(1,recon.series[,-1]), na.rm=T)
    
    # Mean and variance matching over calibration period
    coefficients.calibration <- coefficients[[i]][rep(1:length(breaks),rep(nlayers(Ytrain.brick),length(breaks))),]
    calibration.predictions <- matrix(rowSums(coefficients.calibration * cbind(1,training.series[rep(1:nrow(training.series),length(breaks)),])), ncol=nlayers(Ytrain.brick), byrow=T)
    calibration.predictions.mean <- rowMeans(calibration.predictions)
    calibration.predictions.sd <- apply(calibration.predictions,1,sd)
    
    scalar <- Ytrain.sd[i]/calibration.predictions.sd
    scalar.expand <- rep(scalar,diff(c(breaks,2001)))
    predictions.scaled <- predictions * scalar.expand
    
    predictions.scaled.meaned <- predictions.scaled + rep(Ytrain.mean[i]-(scalar*calibration.predictions.mean),diff(c(breaks,2001)))
    return(predictions.scaled.meaned)
  }, mc.cores=detectCores())
  
  names(recons) <- names(coefficients)
  
  if(raster.out){
    recons.brick <- setValues(brick(Ytrain.brick),do.call("rbind",recons))
    writeRaster(recons.brick,paste(out.dir,name,'_recons.tif',sep=''), datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
    return(recons.brick)
  }else{
    saveRDS(recons,file=paste(out.dir,name,'_recons.rds',sep=''), compress="xz")
    return(recons)
  }
}

recon.errors <- function(models, Ytrain.brick, recon.series, name, out.dir="../OUTPUT/", force.redo=F){
  errors <- lapply(models,'[[','errors')
  error.names <- colnames(errors[[1]])
  
  if(length(models)==ncell(Ytrain.brick)){
    raster.out=T
    if(!force.redo & file.exists(paste(out.dir,name,'_errors_',error.names[1],'.tif',sep=''))){
      return()
    }
  }else{
    raster.out=F
  }
  
  errors <- lapply(errors,function(x){matrix.expand(x,recon.series$YEAR)})
  
  for(error.name in error.names){
    error <- lapply(errors,function(x){x[,error.name]})
    if(raster.out){
      recons.brick <- setValues(brick(Ytrain.brick),do.call("rbind",error))
      writeRaster(recons.brick,paste(out.dir,name,'_errors_',error.name,'.tif',sep=''), datatype="FLT4S", options=c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"), overwrite=T, setStatistics=FALSE)
    }else{
      saveRDS(error, file=paste(out.dir,name,'_errors_',error.name,'.rds',sep=''), compress='xz')
    }
  }
  
  return()
}



getFirstMin <- function(vector){
  min <- which(diff(vector)>0)[1]
  
  if(!is.na(min)){
    return(min)
  }else{
    return(1)
  } 
}


matrix.expand <- function(matrix, new.rows){
  
  missing.rows <- new.rows[!(new.rows %in% intersect(as.numeric(rownames(matrix)),new.rows))]
  new.rows.matrix <- matrix(nrow=length(missing.rows),ncol=ncol(matrix),data=NA)
  colnames(new.rows.matrix) <- colnames(matrix)
  rownames(new.rows.matrix) <- missing.rows
  matrix.expand <- rbind(matrix,new.rows.matrix)
  matrix.expand <- matrix.expand[order(as.numeric(rownames(matrix.expand))),]
  matrix.expand <- na.locf(matrix.expand)
  
  return(matrix.expand)
}

