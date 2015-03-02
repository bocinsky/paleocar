pkgTest("miscTools")

getFolds <- function(folds,i){
  train <- folds$subsets[folds$which != i,]
  validation <- folds$subsets[folds$which == i,]
  if(!is.null(ncol(train))){
    out.list <- lapply(1:ncol(train),FUN=function(x,...){list(TRAIN=train[,x],VALIDATE=validation[,x])})
  }else{
    out.list <- list(TRAIN=train,VALIDATE=validation)
  }
  
  return(list(out.list))
}

# A function to calculate a group of shrinkage linear models, and 
# calculate the AICs and cross-validated statistics of all models
# @returns MODELS A matrix of linear model coefficients
# @returns LAMBDA The correlation shrinkage intensity of the models
# @returns ALL.CAR A numeric vector of raw CAR scores for all Xtrain variables
# @returns PART.CAR A matrix of per-model CAR scores
# @returns AIC A numeric vector of AIC scores
# @returns CORRECTED A boolean term describing whether or not the AIC is corrected for small sample sizes
cv.slm.models <- function(Xtrain, Ytrain, Xvalidation, Yvalidation, lambda=NULL, lambda.var=NULL, regularize=T, max.series=NULL, car.predlist){
  if(is.null(lambda)){
    lambda <- estimate.lambda(cbind(Ytrain,Xtrain))
  }
  
  if(is.null(lambda.var)){
    lambda.var <- estimate.lambda.var(cbind(Ytrain,Xtrain))
  }
  
  # Set the maximum number of series if not provided
  max.series <- ncol(Xtrain)
  
  if(regularize){
    out.slms <- slm.models(Xtrain=Xtrain, Ytrain=Ytrain, predlist=car.predlist[1:length(car.predlist)], lambda=lambda, lambda.var=lambda.var, verbose=F)
  }else{
    out.slms <- slm.models(Xtrain=Xtrain, Ytrain=Ytrain, predlist=car.predlist[1:length(car.predlist)], verbose=F)
  }
  
  # Calculate model estimates for all linear models
  Ytrain.estimate <- apply(out.slms$coefficients,1, function(x){as.numeric(apply(Xtrain,1,function(y,...){as.numeric(sum(y*x[-1]) + x[1])}))})
  Yvalidation.estimate <- apply(out.slms$coefficients,1, function(x){as.numeric(apply(Xvalidation,1,function(y,...){as.numeric(sum(y*x[-1]) + x[1])}))})
  
  Yvalidation.estimate <- meanVarianceTuning(prediction=as.matrix(Ytrain.estimate), training=Ytrain, x.validation=as.matrix(Yvalidation.estimate))
  Ytrain.estimate <- meanVarianceTuning(prediction=as.matrix(Ytrain.estimate), training=Ytrain, x.validation=as.matrix(Ytrain.estimate))
  
  # Calculate the residual sum of squared errors for all linear models
  out.RSS <- apply(Ytrain.estimate,2, function(x){sum((Ytrain-x)^2)})
  
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
  out.slms.AIC <- out.slms.AIC + (2*K*(K+1))/(n-K-1)
  
  
  # Calculate estimated cross-validated model errors
  R2c <- apply(Ytrain.estimate,2,function(x){cor(x,Ytrain)^2})
  R2v <- apply(Yvalidation.estimate,2,function(x){cor(x,Yvalidation)^2})
  RMSE <- apply(Yvalidation.estimate,2,function(x){sqrt(mean((Yvalidation-x)^2))})
  CV_RMSE <- RMSE/mean(Yvalidation)
  N_RMSE <- RMSE/abs(diff(range(Yvalidation)))
  RE <- apply(Yvalidation.estimate,2,function(x){1-(sum((Yvalidation-x)^2)/sum((Yvalidation-mean(Ytrain))^2))})
  CE <- apply(Yvalidation.estimate,2,function(x){1-(sum((Yvalidation-x)^2)/sum((Yvalidation-mean(Yvalidation))^2))})
  
  error.stats <- data.frame(R2c=R2c,R2v=R2v,AIC=out.slms.AIC,RMSE=RMSE,CV_RMSE=CV_RMSE,N_RMSE=N_RMSE,RE=RE,CE=CE)
  
  return(list(MODELS=out.slms$coefficients, ERROR=error.stats))
}

meanVarianceTuning <- function(prediction, training, x.validation){
  prediction <- as.matrix(prediction)
  x.validation <- as.matrix(x.validation)
  
  prediction.mean <- as.numeric(apply(prediction,2,mean))
  training.mean <- mean(training)
  
  prediction.sd <- as.numeric(apply(prediction,2,sd))
  training.sd <- sd(training)
  
  scalar <- training.sd/prediction.sd
  
  y.validation <- t(apply(t(apply(x.validation,1,function(x){x*scalar})),1,function(x){x+(training.mean-(scalar*prediction.mean))}))
  
  return(y.validation)
}