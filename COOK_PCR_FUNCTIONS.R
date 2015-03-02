## A script for performing point-by-point principal components regression, 
## following Cook et al. 1999.
## Cook, Edward R, David M Meko, David W Stahle, and Malcolm K Cleaveland
## 1999 Drought reconstructions for the continental United States. Journal of Climate, 12(4):1145–1162.

getError.pcr.robust <- function(X, Y, fold, correct.aic=T, meanVar=T){
  ## This version performs all pooling using the total validation and calibration data, 
  ## then only performs cross-validation over the final set of principal components.
  ## Thus, the only thing that changes in each cross validation is how many components
  ## enter into the model. This is more in line with the CAR approach, which merely uses
  ## cross-validation to achieve a more accurate measure of model fit than can be achieved
  ## using calibration data alone.
  
  Xtotal <- X
  Ytotal <- Y
  ## Cook et al. 1999 "Level-1" pool: Chronologies within 450km
  ## We already just take the chronologies within the four-corner's states, so nothing to do here.
  
  ## Cook et al. 1999 "Level-2" pool: Significance of raw correlation (alpha = 0.10 cutoff)
  xy.corrs <- apply(Xtotal,2,function(x){cor.test(x,Ytotal,alternative="two.sided",conf.level=0.9)$p.value<=0.1})
  Xtotal.L2 <- Xtotal[,xy.corrs,drop=F]
  
  ## Principal components regression
  ## Calculate the principal components of the training t-r series
  Xtotal.prcomp <- prcomp(Xtotal.L2, center=T, scale=T, retx=T)

  ## Cook et al. 1999 "Level-3" pool: Kaiser–Guttman eigenvalue-1 criterion
  # Estimate the eigenvalues as the squared standard deviations of the eigenvectors (components)
  ev <- Xtotal.prcomp$sdev^2
  Xtotal.prcomp$x <- Xtotal.prcomp$x[,ev>=1]
  
  ## Cook et al. 1999 "Level-4" pool: Stepwise multiple regression, and minimization of the AIC
  # Get correlation of each rotated variable and the scaled outcome variable
  Ytotal.scale <- as.numeric(scale(Ytotal))
  xy.prcomp.corrs <- apply(Xtotal.prcomp$x,2,function(x){cor(x,Ytotal)})
  # Order the correlations
  xy.prcomp.corrs <- order(abs(xy.prcomp.corrs), decreasing=T)
  # Create a stepwise list
  predlist = make.predlist(xy.prcomp.corrs, numpred=1:length(xy.prcomp.corrs), name="CORR")
  
  Xtrain.prcomp <-  Xtotal.prcomp$x[fold$TRAIN,]
  Xvalidation.prcomp <- Xtotal.prcomp$x[fold$VALIDATE,]
  Ytrain=Y[fold$TRAIN]
  Yvalidation=Y[fold$VALIDATE]
  
  # Calculate the OLS estimate of model coefficients
  out.slms <- slm.models(Xtrain=Xtrain.prcomp, Ytrain=Ytrain, predlist, lambda=0, lambda.var=0, verbose=F)
  
  # Calculate the model estimates for all linear models
  Ytrain.estimate <- apply(out.slms$coefficients,1, function(x){as.numeric(apply(Xtrain.prcomp,1,function(y,...){as.numeric(sum(y*x[-1]) + x[1])}))})
  Yvalidation.estimate <- apply(out.slms$coefficients,1, function(x){as.numeric(apply(Xvalidation.prcomp,1,function(y,...){as.numeric(sum(y*x[-1]) + x[1])}))})
  
  if(meanVar){
    Ytrain.estimate <- meanVarianceTuning(prediction=as.matrix(Ytrain.estimate), training=Ytrain, x.validation=Ytrain.estimate)
    Yvalidation.estimate <- meanVarianceTuning(prediction=as.matrix(Ytrain.estimate), training=Ytrain, x.validation=Yvalidation.estimate)
  }
  
  # Calculate the residual sum of squared errors for all linear models
  out.RSS <- apply(Ytrain.estimate, 2, function(x){sum((Ytrain-x)^2)})
  
  # The sample size
  n <- length(Ytrain)
  
  # The number of free parameters
  # This is the number of variables,
  # + 1 for the intercept,
  # + 1 for the estimate of population variance
  K <- out.slms$numpred + 1 + 1
  
  # Calculate the AIC
  out.slms.AIC <- n*log(out.RSS/n) + 2*K
  
  # Calculate the corrected AIC
  if(correct.aic){
    out.slms.AIC <- out.slms.AIC + (2*K*(K+1))/(n-K-1)
  }
  
  # Calculate estimated cross-validated model errors
  R2c <- apply(Ytrain.estimate,2,function(x){cor(x,Ytrain)^2})
  R2v <- apply(Yvalidation.estimate,2,function(x){cor(x,Yvalidation)^2})
  RMSE <- apply(Yvalidation.estimate,2,function(x){sqrt(mean((Yvalidation-x)^2))})
  CV_RMSE <- RMSE/mean(Yvalidation)
  N_RMSE <- RMSE/abs(diff(range(Yvalidation)))
  RE <- apply(Yvalidation.estimate,2,function(x){1-(sum((Yvalidation-x)^2)/sum((Yvalidation-mean(Ytrain))^2))})
  CE <- apply(Yvalidation.estimate,2,function(x){1-(sum((Yvalidation-x)^2)/sum((Yvalidation-mean(Yvalidation))^2))})
  
  error.stats <- data.frame(R2c=R2c,R2v=R2v,AIC=out.slms.AIC,RMSE=RMSE,CV_RMSE=CV_RMSE,N_RMSE=N_RMSE,RE=RE,CE=CE)
  
  return(list(MODELS=out.slms$coefficients, ERROR=error.stats, CORRECTED=correct.aic))
}