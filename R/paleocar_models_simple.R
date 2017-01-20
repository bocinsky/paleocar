globalVariables(c("AICc"))
#' Fit PaleoCAR models to a single predictand
#'
#' This is the primary function for fitting PaleoCAR models to a single predictand using
#' a uniform set of predictors (tree-ring chronologies). \code{\link{carscore}s} are calculated for 
#' each predictor, and models are calculated by adding predictors stepwise to an ordinary least-squares
#' linear model using the \code{\link{lm}} function. Model selection is performed by minizing corrected AIC.
#' This occurs for every unique set of available predictors through time.
#' 
#' See \code{\link{paleocar_models}} for a batch algorithm that is more efficient (though more 
#' computationally complicated) than the algorithm provided here.
#'
#' @param chronologies A matrix of tree ring chronologies, indexed annually.
#' Each chronology is a column. The first column must be labeled "YEAR" and is the calendar year.
#' @param predictands A numeric vector of the predictand (response) variable.
#' @param calibration.years An integer vector of years corresponding to the layers in the \code{predictands} brick.
#' @param prediction.years An optional integer vector of years for the reconstruction.
#' If missing, defaults to the total years present in \code{chronologies}.
#' @param verbose Logical, display status messages during run.
#' @return A named list containing
#' \itemize{
#'   \item{\code{models}  A \code{\link{data.table}} giving the fitted coefficients, LOOCV errors, and AICc values through time for each predictand.}
#'   \item{\code{predictand}  A numeric vector of the predictand, as provided.}
#'   \item{\code{predictor.matrix}  A matrix of predictors for calibration; \code{chronologies} cropped to \code{calibration.years}.}
#'   \item{\code{reconstruction.matrix}  A matrix of predictors for reconstruction; \code{chronologies} cropped to \code{prediction.years}, or all of \code{chronologies} if \code{prediction.years==NULL}.}
#' }
#' @importFrom stats lm
#' @export
paleocar_models_simple <- function(predictands,
                            chronologies,
                            calibration.years,
                            prediction.years=NULL,
                            verbose=F){
  predictor.matrix <- get_predictor_matrix(chronologies, calibration.years)
  
  maxPreds <- nrow(predictor.matrix)-5
  
  reconstruction.matrix <- get_reconstruction_matrix(chronologies, prediction.years)
  reconstruction.matrix <- reconstruction.matrix[,colnames(predictor.matrix)]
  
  predlist <- get_predlist(reconstruction.matrix)
  prednums <- rowSums(predlist, na.rm=T)
  prednums[prednums>maxPreds] <- maxPreds
  
  predyears <- as.numeric(rownames(predlist))
  hinge.year <- min(predyears[predyears > max(calibration.years)], max(calibration.years))
  
  carscores <- care::carscore(Xtrain=predictor.matrix, Ytrain=predictands, verbose=F)
  carscores.ranks <- rank(1-(carscores^2))
  names(carscores.ranks) <- colnames(predictor.matrix)
  carscores <- carscores.ranks
  rm(carscores.ranks); gc(); gc()
  
  predlist <- predlist==1
  predlist[is.na(predlist)] <- F
  
  models <- lapply(1:nrow(predlist),function(i){
    preds <- predlist[i,]
    car.preds <- carscores
    car.preds[!preds] <- NA
    car.preds <- rank(car.preds, na.last='keep')
    
    numModels <- min(prednums[i],maxPreds)
    
    models <- matrix(FALSE,nrow=numModels,ncol=length(car.preds))
    for(i in 1:numModels){
      models[i,car.preds<=i] <- TRUE
    }
    colnames(models) <- names(car.preds)
    return(models)
  })
  names(models) <- rownames(predlist)
  
  new.t <- Sys.time()
  all.lms <- lapply(models,function(year.models){
    year.lms <- data.table::rbindlist(apply(year.models,1,function(this.model){
      model.lm <- stats::lm(predictands~predictor.matrix[,this.model,drop=F])
      # Get model errors
      model.errors <- data.table::data.table(numPreds=ncol(predictor.matrix[,this.model,drop=F]),t(forecast::CV(model.lm)[c("CV","AICc")]))
      
      # Get coefficients
      coefs <- data.table::data.table(t(as.matrix(model.lm$coefficients)))
      data.table::setnames(coefs,c("Intercept",colnames(predictor.matrix[,this.model,drop=F])))
      
      coefs <- lapply(1:nrow(coefs),function(j){
        out <- as.numeric(coefs[j])
        names(out) <- colnames(coefs)
        return(out)
      })
      # names(coefs) <- cells
      
      # coefs <- coefs[order(as.integer(names(coefs)))]
      
      if(length(coefs)==1){
        model.errors[,coefs:=list(coefs)]
      }else{
        model.errors[,coefs:=coefs]
      }
      model.errors[,model:=0]
      model.errors[,cell:=1]
      data.table::setcolorder(model.errors,c("cell","model","numPreds","CV","AICc","coefs"))
      data.table::setkey(model.errors,cell,model)
      
      
      
      
      
      
      # coefs[,setdiff(colnames(predictor.matrix),names(coefs)):=NA]
      # data.table::setcolorder(coefs,c("Intercept",colnames(predictor.matrix)))
      # 
      # if(length(coefs)==1){
      #   model.errors[,coefs:=list(coefs)]
      # }else{
      #   model.errors[,coefs:=coefs]
      # }
      # 
      # out <- cbind(model.errors,coefs)
      return(model.errors)
    }))
    
    return(year.lms)
  })
  
  all.lms <- lapply(1:length(all.lms),function(i){
    lms <- all.lms[[i]]
    lms[,year:=as.numeric(names(all.lms)[i])]
    data.table::setcolorder(lms,c("cell","year","model","numPreds","CV","AICc","coefs"))
    
    if(nrow(lms) == 1){
      lms <- lms
    }else if(nrow(lms) == 2){
      lms <- lms[which(lms$AICc == min(lms$AICc)),]
    }else if(nrow(lms) > 2){
      lms <- lms[which(sign(diff(lms$AICc)) == 1)[1],]
    }

    return(lms)
  })
  
  all.lms <- data.table::rbindlist(all.lms)
  data.table::setkey(all.lms,year,AICc)

#   if(verbose) cat("\nCalc lms:", Sys.time()-new.t)
  
  all.lms.below <- all.lms[year<hinge.year,]
  all.lms.above <- all.lms[year>=hinge.year,]
  data.table::setorder(all.lms.below, year,AICc)
  data.table::setorder(all.lms.above,-year,AICc)
  
  data.table::setkey(all.lms.below,NULL)
  data.table::setkey(all.lms.above,NULL)
  
  below.remove <- all.lms.below[,list(make_monotonic(AICc),year)]
  above.remove <- all.lms.above[,list(make_monotonic(AICc),year)]
  
  all.lms.below <- all.lms.below[below.remove$V1]
  all.lms.above <- all.lms.above[above.remove$V1]
  
  all.lms <- rbind(all.lms.below,all.lms.above)
  rm(all.lms.below,all.lms.above,below.remove,above.remove);gc();gc()
  
  data.table::setkey(all.lms,cell,year)
  
  allModels <- list(models = all.lms,
                    predictands = predictands,
                    predictor.matrix = predictor.matrix,
                    reconstruction.matrix = reconstruction.matrix,
                    carscores = data.table::data.table(t(carscores)))
    
  return(allModels)
}