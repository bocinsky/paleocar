#' Fit PaleoCAR models to a (potentially large) set of predictands
#'
#' This is the primary function for fitting PaleoCAR models to a large set of predictands using
#' a uniform set of predictors (tree-ring chronologies). For each predictand, \code{\link{carscore}s} are calculated for 
#' each predictor, and models are calculated by adding predictors stepwise to an ordinary least-squares
#' linear model using the \code{\link{lm}} function. Model selection is performed by minizing corrected AIC.
#' This occurs for every unique set of available predictors through time.
#' 
#' See \code{\link{paleoCAR.models}} for a simplified algorithm that fits a PaleoCAR model for 
#' a single predictand. The batch algorithm gains efficiencies by minimizing the number of unique model calculations
#' (using multiple-response models) and by halting calculations for individual predictands once they
#' stabilize.
#'
#' @param chronologies A matrix of tree ring chronologies, indexed annually.
#' Each chronology is a column. The first column must be labeled "YEAR" and is the calendar year.
#' @param predictands A RasterBrick, RasterStack, or matrix of the numeric predictand (response) variable.
#' @param calibration.years An integer vector of years corresponding to the layers in the \code{predictands} brick.
#' @param prediction.years An optional integer vector of years for the reconstruction.
#' If missing, defaults to the total years present in \code{chronologies}.
#' @param label A character label for the reconstruction, for saving.
#' @param out.dir The directory to which output is to be saved.
#' @param force.redo Logical, should all computations be re-computed?
#' @param verbose Logical, display status messages during run.
#' @return A named list containing
#' \itemize{
#'   \item{\code{models}  A \code{\link{data.table}} giving the fitted coefficients, LOOCV errors, and AICc values through time for each predictand.}
#'   \item{\code{predictands}  A matrix of RasterBrick of predictands, as provided.}
#'   \item{\code{predictor.matrix}  A matrix of predictors for calibration; \code{chronologies} cropped to \code{calibration.years}.}
#'   \item{\code{reconstruction.matrix}  A matrix of predictors for reconstruction; \code{chronologies} cropped to \code{prediction.years}, or all of \code{chronologies} if \code{prediction.years==NULL}.}
#' }
paleoCAR.models.batch <- function(chronologies, predictands, calibration.years, prediction.years=NULL, label, out.dir="./OUTPUT/", force.redo=F, verbose=F){
  if(!force.redo & file.exists(paste(out.dir,label,'.models.rds',sep=''))){
    allModels <- readRDS(paste(out.dir,label,'.models.rds',sep=''))
    return(allModels)
  }
  
  predictor.matrix <- getPredictorMatrix(chronologies, calibration.years)
  
  maxPreds <- nrow(predictor.matrix)-5
  
  if(class(predictands) %in% c("RasterBrick","RasterStack")){
    predictand.matrix <- t(raster::values(predictands))
    colnames(predictand.matrix) <- as.character(1:ncell(predictands)) 
  }

  reconstruction.matrix <- getReconstructionMatrix(chronologies, prediction.years)
  reconstruction.matrix <- reconstruction.matrix[,colnames(predictor.matrix)]
  
  predlist <- getPredlist(reconstruction.matrix)
  prednums <- rowSums(predlist, na.rm=T)
  prednums[prednums>maxPreds] <- maxPreds
  
  predyears <- as.numeric(rownames(predlist))
  hinge.year <- min(predyears[predyears>max(calibration.years)])
  
  if(!force.redo & file.exists(paste(out.dir,label,".carscores.Rds",sep='')) && ncol(readRDS(paste(out.dir,label,".carscores.Rds",sep='')))==ncol(predictand.matrix)){
    carscores <- readRDS(paste(out.dir,label,".carscores.Rds",sep=''))
  }else{
    carscores <- carscore.batch(predictand.matrix=predictand.matrix, predictor.matrix=predictor.matrix)
    carscores.ranks <- matrixStats::colRanks(1-(carscores^2), preserveShape=T)
    rownames(carscores.ranks) <- rownames(carscores)
    colnames(carscores.ranks) <- colnames(carscores)
    carscores <- carscores.ranks
    rm(carscores.ranks); gc(); gc()
    saveRDS(carscores, paste(out.dir,label,".carscores.Rds",sep=''), compress='xz')
  }
  
  carscores <- data.table(t(carscores))
  predlist <- predlist==1
  predlist[is.na(predlist)] <- F
  
  allModels <- data.table(cell=numeric(),year=numeric(),numPreds=numeric(),CV=numeric(),AICc=numeric(),Intercept=numeric(),predictor.matrix[0,])
  completed.cells <- vector("logical",ncol(predictand.matrix))
  times <- vector('numeric',max(prednums))
  
  for(i in 1:maxPreds){
    if(verbose) cat("\nCalculating models of size",i)
    t <- Sys.time()
    
    new.t <- Sys.time()
    preds <- predlist[prednums>=i,]
    
    ## MATCHES 1
    matches <- lapply(rownames(preds),function(year){
      pred.names <- colnames(preds)[preds[year,]]
      models <- data.table::data.table(matrixStats::rowRanks(as.matrix(carscores[,pred.names,with = F]))<=i)
      data.table::setnames(models,pred.names)
      models[,cell:=1:nrow(models)]
      models <- models[cell %in% which(!completed.cells)]
      data.table::setorderv(models,pred.names,rep(-1,length(pred.names)))
      models.duplicates <- !duplicated(models,by=pred.names)
      models.matches <- cumsum(models.duplicates)
      names(models.matches) <- models$cell
      models.matches <- models.matches[order(as.numeric(names(models.matches)))]
      models <- models[models.duplicates,]
      models[,cell:=NULL]
      
      return(list(models=models,matches=models.matches))
    })
    gc();gc()
    
    models <- lapply(matches,'[[','models')
    matches <- lapply(matches,'[[','matches')
    names(matches) <- rownames(preds)
    names(models) <- rownames(preds)
    
    nummodels <- c(0,cumsum(sapply(models,nrow)))[0-(length(models)+1)]
    names(nummodels) <- names(models)
    
    models <- rbindlist(models,fill=T)
    pred.names <- names(models)[order(names(models))]
    setcolorder(models,pred.names)
    models[,model:=1:nrow(models)]
    
    if(nrow(models)==0) break
    ## Replace all NAs with FALSE (we aren't going to use those predictors, right?)
    for (this.model in seq_along(models)){
      set(models, i=which(is.na(models[[this.model]])), j=this.model, value=F)
    }
    
    matches <- lapply(1:length(matches),function(count){return(list(year=rep(names(matches)[count],length(matches[[count]])),matches=matches[[count]]+nummodels[count]))})
    
    matches.years <- unlist(lapply(matches,'[[','year'))
    matches.matches <- unlist(lapply(matches,'[[','matches'))
    
    matches <- do.call(rbind,list(cell=as.integer(names(matches.matches)),year=as.integer(matches.years),model=as.integer(matches.matches)))
    colnames(matches) <- NULL
    
    rm(matches.years, matches.matches)
    gc();gc()
    
    models <- setorderv(models,pred.names,rep(-1,length(pred.names)))
    models.duplicates <- !duplicated(models,by=pred.names)
    models.matches <- cumsum(models.duplicates)
    names(models.matches) <- models$model
    models.matches <- models.matches[order(as.numeric(names(models.matches)))]
    models <- models[models.duplicates,]
    models[,model:=NULL]
    
    matches.match <- match(matches['model',],as.numeric(names(models.matches)))
    matches['model',] <- models.matches[as.character(matches.match)]
    
    matches <- data.table(t(matches))
    matches <- setorderv(matches,c('cell','year','model'),c(1,1,1))
    
    matches.below <- matches[year<hinge.year,]
    matches.below <- unique(matches.below,by=c('cell','model'))
    
    matches.above <- matches[year>=hinge.year,]
    matches.above <- setorderv(matches.above,c('cell','year'),c(1,-1))
    matches.above <- unique(matches.above,by=c('cell','model'))
    matches.above <- setorderv(matches.above,c('cell','year','model'),c(1,1,1))
    
    matches <- rbind(matches.below,matches.above)
    matches <- setorderv(matches,c('cell','year','model'),c(1,1,1))
    setkey(matches,model)
    
    rm(models.duplicates,models.matches,matches.match,matches.below,matches.above); gc(); gc()
    if(verbose) cat("\nDefine models:", Sys.time()-new.t, "nrow(models):",nrow(models),"ncells:",sum(!completed.cells))
    
    new.t <- Sys.time()
    all.lms <- lapply(1:nrow(models),function(this.model){
      #       cat("\nModel",this.model,"of",nrow(models))
      cells <- unique(matches[.(this.model),cell])
      model.mlm <- lm(predictand.matrix[,cells,drop=F]~predictor.matrix[,as.logical(models[this.model]),drop=F])
      
      # Get model errors
      model.errors <- data.table::data.table(CV.mlm(model.mlm))
      model.errors[,cell:=cells]
      #       model.errors[,model:=this.model]
      data.table::setkey(model.errors,cell)
      
      # Get coefficients
      coefs <- data.table::data.table(t(as.matrix(model.mlm$coefficients)))
      data.table::setnames(coefs,c("Intercept",colnames(predictor.matrix[,as.logical(models[this.model]),drop=F])))
      coefs[,cell:=cells]
      #       coefs[,model:=this.model]
      data.table::setkey(coefs,cell)
      
      out <- model.errors[coefs]
      out[,model:=this.model]
      data.table::setcolorder(out,c("cell","model","CV","AICc","Intercept",colnames(predictor.matrix[,as.logical(models[this.model]),drop=F])))
      data.table::setkey(out,cell,model)
      
      return(out)
    })
    if(verbose) cat("\nCalc lms:", Sys.time()-new.t)
    rm(models); gc(); gc()
    
    ## LMS SIMPLIFY PREP
    new.t <- Sys.time()
    all.lms <- rbindlist(all.lms,fill=T)
    setkey(all.lms,cell,model)
    
    setkey(matches,cell,model)
    
    all.lms <- all.lms[matches]
    rm(matches)
    
    setkey(all.lms,cell,year,model)
    setkey(all.lms,cell,year)
    
    all.lms.below <- all.lms[year<hinge.year,]
    all.lms.above <- all.lms[year>=hinge.year,]
    setorder(all.lms.above,cell,-year)
    
    setkey(all.lms.below,NULL)
    setkey(all.lms.above,NULL)
    
    below.remove <- all.lms.below[,.(makeMonotonic(AICc),year),by=cell]
    above.remove <- all.lms.above[,.(makeMonotonic(AICc),year),by=cell]
    
    all.lms <- rbind(all.lms.below[below.remove$V1],all.lms.above[above.remove$V1])
    rm(all.lms.below,all.lms.above,below.remove,above.remove);gc();gc()
    
    setkey(all.lms,cell,year)
    
    all.lms[,setdiff(colnames(predictor.matrix),names(all.lms)):=NA]
    all.lms[,model:=NULL]
    all.lms[,numPreds:=i]
    setcolorder(all.lms,c("cell","year","numPreds","CV","AICc","Intercept",colnames(predictor.matrix)))
    
    allModels <- rbind(allModels,all.lms)
    
    setkey(allModels,cell,year,numPreds)
    setorder(allModels,cell,year,AICc)
    allModels <- allModels[allModels[,!duplicated(year),by=cell]$V1,]
    
    allModels.below <- allModels[year<hinge.year,]
    allModels.above <- allModels[year>=hinge.year,]
    setorder(allModels.above,cell,-year)
    
    setkey(allModels.below,NULL)
    setkey(allModels.above,NULL)
    
    below.remove <- allModels.below[,.(makeMonotonic(AICc),year),by=cell]
    above.remove <- allModels.above[,.(makeMonotonic(AICc),year),by=cell]
    
    #     remove <- rbind(below.remove,above.remove)
    #     remove <- remove[V1==F]
    #     cellYears[cbind(remove$cell,match(as.character(remove$year),colnames(cellYears)))] <- F
    
    allModels <- rbind(allModels.below[below.remove$V1],allModels.above[above.remove$V1])
    rm(allModels.below,allModels.above,below.remove,above.remove);gc();gc()
    
    allModels[,cell:=as.numeric(cell)]
    setkey(allModels,cell,year,numPreds)
    if(verbose) cat("\nClean lms:", Sys.time()-new.t)
    
    if(i>1){
      new.t <- Sys.time()
      completed.cells <- merge(allModels[,1:5,with=F],models.last.iter, all=T)
      completed.cells[,CV.change:=CV.x-CV.y]
      completed.cells[,AICc.change:=AICc.x-AICc.y]
      completed.cells <- completed.cells[,sum(CV.change)+sum(AICc.change),by=cell]$V1==0
      completed.cells[is.na(completed.cells)] <- F
      if(verbose) cat("\nCalc completed:", Sys.time()-new.t)
    }  
    
    models.last.iter <- allModels[,1:5,with=F]
    
    time <- Sys.time()-t
    if(verbose) cat("\nTime:",time,"\n")
    times[i] <- time
    ## 
    
  }
  
  allModels <- list(models=allModels, predictands=predictands, predictor.matrix=predictor.matrix, reconstruction.matrix=reconstruction.matrix)
  
  saveRDS(allModels,file=paste(out.dir,label,'.models.rds',sep=''), compress='xz')
  
  return(allModels)
}
