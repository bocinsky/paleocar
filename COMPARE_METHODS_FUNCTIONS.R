## Other methods tested include
## single-chronology regression (nearest chronology),
## and Principal Components Regression (after Cook et al. 1999).
cv.compare.methods <- function(X, Y, regularize=T, K=3, R=1, type="consecutive"){
  if(nrow(X)!=length(Y)){
    stop("Number of X observations does not match the length of the response variable!")
  }
  
  if(K==1){
    folds.out <- vector('list',1)
    folds.out[[1]]$TRAIN <- 1:length(Y)
    folds.out[[1]]$VALIDATE <- 1:length(Y)
  }else{
    folds <- cvFolds(length(Y), K=K, R=R, type=type)
    folds.out <- unlist(lapply(1:folds$K,FUN=function(x,...){ getFolds(folds=folds,i=x) }), recursive=F)
    if(type!="consecutive" & length(folds.out)<=R & K!=length(Y)){
      folds.out <- unlist(folds.out, recursive=F)
    }
  }
  
  
  
  robust.car <- cv.select.slm.models(X=X,Y=Y,fold.list=folds.out)$CV.ERRORS
  robust.car.mean <- Reduce('+',robust.car)/length(robust.car)
  robust.car.min.rmse <- robust.car.mean[robust.car.mean[,'RMSE'] == min(robust.car.mean[,'RMSE']),]
  robust.car.min.aic <- robust.car.mean[robust.car.mean[,'AIC'] == min(robust.car.mean[,'AIC']),]
  
  
  robust.mcor <- cv.select.slm.models(X=X, Y=Y, lambda=0, lambda.var=0, diagonal=T,fold.list=folds.out)$CV.ERRORS
  robust.mcor.mean <- Reduce('+',robust.mcor)/length(robust.mcor)
  robust.mcor.min.rmse <- robust.mcor.mean[robust.mcor.mean[,'RMSE'] == min(robust.mcor.mean[,'RMSE']),]
  robust.mcor.min.aic <- robust.mcor.mean[robust.mcor.mean[,'AIC'] == min(robust.mcor.mean[,'AIC']),]
  
  
  cook.pcr <- lapply(folds.out,function(fold){getError.pcr.robust(X=X, Y=Y, fold=fold, correct.aic=T, meanVar=F)})
  cook.pcr <- lapply(cook.pcr,'[[','ERROR')
  cook.pcr.mean <- Reduce('+',cook.pcr)/length(cook.pcr)
  cook.pcr.min.rmse <- cook.pcr.mean[cook.pcr.mean[,'RMSE'] == min(cook.pcr.mean[,'RMSE']),]
  cook.pcr.min.aic <- cook.pcr.mean[cook.pcr.mean[,'AIC'] == min(cook.pcr.mean[,'AIC']),]
  
  
  out.list <- list(robust.car.min.rmse,robust.car.min.aic,robust.mcor.min.rmse,robust.mcor.min.aic,cook.pcr.min.rmse,cook.pcr.min.aic)
  out.table <- do.call(rbind, out.list)
  rownames(out.table) <- c("robust.car.min.rmse","robust.car.min.aic","robust.mcor.min.rmse","robust.mcor.min.aic","cook.pcr.min.rmse","cook.pcr.min.aic")
  
  return(out.table)
}