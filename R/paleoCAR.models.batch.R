#' Fit PaleoCAR models to a (potentially large) set of predictands
#'
#' This is the primary function for fitting PaleoCAR models to a large set of predictands using
#' a uniform set of predictors (tree-ring chronologies). For each predictand, \code{\link{carscore}s} are calculated for 
#' each predictor, and models are calculated by adding predictors stepwise to an ordinary least-squares
#' linear model using the \code{\link{lm}} function. Model selection is performed by minimizing corrected AIC.
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
#' @param min.width integer, indicating the minimum number of tree-ring samples allowed for that year of a chronology to be valid.
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
paleoCAR.models.batch <- function(chronologies, predictands, calibration.years, prediction.years=NULL, min.width=NULL, label, out.dir="./OUTPUT/", force.redo=F, verbose=F){
  # if(verbose) cat("Calculating PaleoCAR models\n")
  if(!force.redo & file.exists(paste(out.dir,label,'.models.rds',sep=''))){
    allModels <- readRDS(paste(out.dir,label,'.models.rds',sep=''))
    return(allModels)
  }
  
  t <- Sys.time()
  predictor.matrix <- getPredictorMatrix(chronologies=chronologies, calibration.years=calibration.years, min.width=min.width)
  
  null.cells <- which(is.na(raster::getValues(predictands[[1]])))
  
  maxPreds <- nrow(predictor.matrix)-5
  
  if(class(predictands) %in% c("RasterBrick","RasterStack")){
    predictand.matrix <- t(raster::values(predictands))
    colnames(predictand.matrix) <- as.character(1:raster::ncell(predictands)) 
  }
  
  reconstruction.matrix <- getReconstructionMatrix(chronologies=chronologies, reconstruction.years=prediction.years, min.width=min.width)
  reconstruction.matrix <- reconstruction.matrix[,colnames(predictor.matrix)]
  
  predlist <- getPredlist(reconstruction.matrix)
  prednums <- rowSums(predlist, na.rm=T)
  prednums[prednums>maxPreds] <- maxPreds
  
  predyears <- as.numeric(rownames(predlist))
  hinge.year <- min(predyears[predyears>max(calibration.years)])
  
#   if(!force.redo & file.exists(paste(out.dir,label,".carscores.Rds",sep=''))){
#     carscores <- readRDS(paste(out.dir,label,".carscores.Rds",sep=''))
#   }else{
    carscores <- carscore.batch(predictand.matrix=predictand.matrix, predictor.matrix=predictor.matrix)
    carscores.ranks <- matrixStats::colRanks(1-(carscores^2), preserveShape=T)
    rownames(carscores.ranks) <- rownames(carscores)
    colnames(carscores.ranks) <- colnames(carscores)
    carscores <- carscores.ranks
    rm(carscores.ranks); gc(); gc()
    carscores <- data.table::data.table(t(carscores))
    saveRDS(carscores, paste(out.dir,label,".carscores.Rds",sep=''), compress='xz')
  # }
  
  if(verbose) cat("\nPrepare data and calculate CAR scores:", round(difftime(Sys.time(),t,units='mins'),digits=2),"minutes\n")
  
  allModels <- data.table(cell=numeric(),year=numeric(),model=numeric(),numPreds=numeric(),CV=numeric(),AICc=numeric(),coefs=numeric())
  complete.cell.years <- data.table(cell=numeric(),year=numeric())
  times <- vector('numeric',max(prednums))
  
  for(i in 1:maxPreds){
    if(verbose) cat("\nCalculating models of size",i)
    t <- Sys.time()
    
    new.t <- Sys.time()
    preds <- predlist[prednums>=i,]
    
    ## MATCHES 1
    matches <- lapply(rownames(preds),function(this.year){
      # cat(this.year,'\n')
      complete.cells <- c(complete.cell.years[year==as.numeric(this.year), cell],null.cells)
      
      pred.names <- colnames(preds)[preds[this.year,]]
      if(all((rownames(carscores) %in% complete.cells))) return(NULL)
      models <- data.table::data.table(matrixStats::rowRanks(as.matrix(carscores[!(rownames(carscores) %in% complete.cells),pred.names,with = F]))<=i)
      data.table::setnames(models,pred.names)
      models[,cell:=as.numeric(rownames(carscores)[!(rownames(carscores) %in% complete.cells)])]
      # models <- models[cell %in% which(!completed.cells)]
      data.table::setorderv(models,pred.names,rep(-1,length(pred.names)))
      models.duplicates <- !duplicated(models,by=pred.names)
      models.matches <- cumsum(models.duplicates)
      names(models.matches) <- models$cell
      models.matches <- models.matches[order(as.numeric(names(models.matches)))]
      models <- models[models.duplicates,]
      models[,cell:=NULL]
      blank.models <- !apply(models,2,any)
      blank.models <- names(blank.models)[blank.models]
      if(length(blank.models)>0){
        models[,blank.models:=NULL, with=F]
      }
      
      return(list(models=models,matches=models.matches))
    })
    gc();gc()
    
    match.names <- rownames(preds)[!sapply(matches,is.null)]
    matches <- matches[!sapply(matches,is.null)]
    
    models <- lapply(matches,'[[','models')
    matches <- lapply(matches,'[[','matches')
    names(matches) <- match.names
    names(models) <- match.names
    
    nummodels <- c(0,cumsum(sapply(models,nrow)))[0-(length(models)+1)]
    names(nummodels) <- names(models)
    
    models <- rbindlist(models,fill=T)
    pred.names <- sort(names(models))
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
    matches <- data.table(t(matches))
    
    rm(matches.years, matches.matches)
    gc();gc()
    
    models <- setorderv(models,pred.names,rep(-1,length(pred.names)))
    models.duplicates <- !duplicated(models,by=pred.names)
    models.matches <- cumsum(models.duplicates)
    names(models.matches) <- models$model
    models.matches <- models.matches[order(as.numeric(names(models.matches)))]
    models <- models[models.duplicates,]
    models[,model:=NULL]
    models.matches <- data.table(model=as.integer(names(models.matches)),new.model=models.matches)
    
    setkey(matches,model)
    setkey(models.matches,model)
    
    matches <- matches[models.matches]
    matches[,model:=NULL]
    setnames(matches,c('cell','year','model'))
    matches <- setorderv(matches,c('cell','year','model'),c(1,1,1))
    
    rm(models.duplicates,models.matches); gc(); gc()
    
    if(nrow(complete.cell.years)>0){
      which.matches <- merge(matches,complete.cell.years,by=c("cell","year"), all.x=T)[is.na(model.y)]
      matches <- which.matches[,.(cell,year,model.x)]
      setnames(matches,c("cell",'year','model'))
    }
    
    setkey(matches,model)
    
    if(verbose) cat("\nDefine models:", round(difftime(Sys.time(),new.t,units='mins'),digits=2),"minutes")
    
    new.t <- Sys.time()
    all.lms <- lapply(1:nrow(models),function(this.model){
      if(!(this.model %in% matches[['model']])) return(NULL)
      cells <- unique(matches[.(this.model),cell])
      
      predictand.names <- names(models[this.model])[which(as.logical(models[this.model]))]
      model.mlm <- lm(predictand.matrix[,cells,drop=F]~predictor.matrix[,predictand.names,drop=F])
      
      # Get model errors
      model.errors <- data.table::data.table(CV.mlm(model.mlm))
      model.errors[,cell:=cells]

      data.table::setkey(model.errors,cell)
      
      # Get coefficients
      
      coefs <- data.table::data.table(t(as.matrix(model.mlm$coefficients)))
      data.table::setnames(coefs,c("Intercept",colnames(predictor.matrix[,predictand.names,drop=F])))
      coefs <- lapply(1:nrow(coefs),function(j){
        out <- as.numeric(coefs[j])
        names(out) <- colnames(coefs)
        return(out)
      })
      names(coefs) <- cells

      coefs <- coefs[order(as.integer(names(coefs)))]

      if(length(coefs)==1){
        model.errors[,coefs:=list(coefs)]
      }else{
        model.errors[,coefs:=coefs]
      }
      model.errors[,model:=this.model]
      data.table::setcolorder(model.errors,c("cell","model","CV","AICc","coefs"))
      data.table::setkey(model.errors,cell,model)
      
      return(model.errors)
    })
    all.lms <- all.lms[!sapply(all.lms,is.null)]
    if(verbose) cat("\nCalculate",length(all.lms),"linear models:", round(difftime(Sys.time(),new.t,units='mins'),digits=2),"minutes")
    rm(models); gc(); gc()
    
    ## LMS SIMPLIFY PREP
    new.t <- Sys.time()

    all.lms <- rbindlist(all.lms,fill=T)
    setkey(all.lms,cell,model)
    all.lms[,numPreds:=i]
    
    if(nrow(allModels)==0){
      setkey(all.lms,cell,model)
      setkey(matches,cell,model)
      allModels <- matches[all.lms]

      setcolorder(allModels,c("cell","year","model","numPreds","CV","AICc","coefs"))
    }
    
    all.lms.stats <- all.lms[,.(cell,model,AICc)]
    setkey(all.lms.stats,cell,model)
    setkey(matches,cell,model)
    
    all.lms.stats <- matches[all.lms.stats]
    rm(matches)
    
    setkey(all.lms.stats,cell,year,model)
    
    allModels.stats <- allModels[,.(cell,year,model,AICc)]
    
    allModels.stats <- rbind(allModels.stats,all.lms.stats,fill=T)
    
    setkey(allModels.stats,cell,year)
    setorder(allModels.stats,cell,year,AICc)
    allModels.stats <- allModels.stats[allModels.stats[,!duplicated(year),by=cell]$V1,]
    
    setkey(allModels,cell,year,model)
    if(nrow(allModels.stats[model==0,.(cell,year,model)])>0){
      allModels.old <- merge(allModels.stats[model==0,.(cell,year,model)],allModels,by=c("cell","year","model"), all=T)
    }else{
      allModels.old <- data.table(cell=numeric(),year=numeric(),model=numeric(),numPreds=numeric(),CV=numeric(),AICc=numeric(),coefs=numeric())
    }
    
    if(nrow(allModels.stats[model!=0,.(cell,year,model)])>0){
      allModels.new <- merge(allModels.stats[model!=0,.(cell,year,model)],all.lms.stats,by=c("cell","year","model"), all=T)
      allModels.new <- merge(allModels.new[,.(cell,year,model)],all.lms,by=c("cell","model"),all=T)
      setcolorder(allModels.new,c("cell","year","model","numPreds","CV","AICc","coefs"))
    }else{
      allModels.new <- data.table(cell=numeric(),year=numeric(),model=numeric(),numPreds=numeric(),CV=numeric(),AICc=numeric(),coefs=numeric())
    }
    
    allModels <- rbind(allModels.old,allModels.new)
    
    setkey(allModels,cell,year)
    setorder(allModels,cell,year,AICc)
    allModels <- allModels[allModels[,!duplicated(year),by=cell]$V1,]
    
    complete.cell.years <- allModels[model==0,.(cell,year,model)]
    
    allModels[,model:=0]

    setkey(allModels,cell,year,numPreds)
    
    if(verbose) cat("\nClean linear models:", round(difftime(Sys.time(),new.t,units='mins'),digits=2),"minutes")
    
    time <- difftime(Sys.time(),t,units='mins')
    if(verbose) cat("\nTotal modeling time:",round(time,digits=2),"minutes\n")
    times[i] <- time
    if((nrow(allModels)-nrow(complete.cell.years))==0) break
    if(verbose) cat(nrow(allModels)-nrow(complete.cell.years),"cell-years remaining\n")
    ## 
    
  }
  
  if(verbose) cat("\nTotal Modeling Time:",sum(times),"minutes\n")
  
  t <- Sys.time()
  get.coef.names <- function(year,model,coefs,numPreds,CV,AICc){
    the.coefs <- lapply(coefs,function(x){names(x)[-1]})
    test.out <- lapply(the.coefs,function(x){which(rowSums(predlist[,x,drop=F])==length(x))})
    test.out <- data.table(model=rep(1:length(test.out),times=sapply(test.out,length)), year=unlist(test.out))
    # test.out <- lapply(1:length(test.out),function(i){data.table(model=i,year=test.out[[i]])})
    # test.out <- rbindlist(test.out)
    test.out <- test.out[which(!duplicated(test.out[,year]))]
    setkey(test.out,year)
    return(list(year=sort(year),model=model[test.out$model],numPreds=numPreds[test.out$model],CV=CV[test.out$model],AICc=AICc[test.out$model],coefs=coefs[test.out$model]))
  }
  
  allModels <- allModels[order(AICc)][,get.coef.names(year,model,coefs,numPreds,CV,AICc),by=cell]
  setkey(allModels,cell,year)
  
  allModels <- allModels[allModels[,!FedData::sequential_duplicated(CV),by=cell]$V1,]
  
  time <- difftime(Sys.time(),t,units='mins')
  if(verbose) cat("\nOptimizing models:",round(time,digits=2),"minutes\n")
  
  allModels <- list(models=allModels, predictands=predictands, predictor.matrix=predictor.matrix, reconstruction.matrix=reconstruction.matrix)
  
  saveRDS(allModels,file=paste(out.dir,label,'.models.rds',sep=''), compress='xz')

  return(allModels)
}
