globalVariables(c("model", "model.y", "model.x", "AdjR2", "AICc", "coefs", "numPreds", "."))
#' Fit PaleoCAR models to a (potentially large) set of predictands
#'
#' This is the primary function for fitting PaleoCAR models to a large set of predictands using
#' a uniform set of predictors (tree-ring chronologies). For each predictand, \code{\link{carscore}s} are calculated for
#' each predictor, and models are calculated by adding predictors stepwise to an ordinary least-squares
#' linear model using the \code{\link{lm}} function. Model selection is performed by minimizing corrected AIC.
#' This occurs for every unique set of available predictors through time.
#'
#' See \code{\link{paleocar_models}} for a simplified algorithm that fits a PaleoCAR model for
#' a single predictand. The batch algorithm gains efficiencies by minimizing the number of unique model calculations
#' (using multiple-response models) and by halting calculations for individual predictands once they
#' stabilize.
#'
#' @param chronologies A matrix of tree ring chronologies, indexed annually.
#' Each chronology is a column. The first column must be labeled "YEAR" and is the calendar year.
#' @param predictands A RasterBrick, RasterStack, matrix, or vector of the numeric predictand (response) variable.
#' @param calibration.years An integer vector of years corresponding to the layers in the \code{predictands} brick.
#' @param prediction.years An optional integer vector of years for the reconstruction.
#' If missing, defaults to the total years present in \code{chronologies}.
#' @param min.width integer, indicating the minimum number of tree-ring samples allowed for that year of a chronology to be valid.
#' @param verbose Logical, display status messages during run.
#' @param ... Further arguments to be passed to other functions.
#' @return A named list containing
#' \itemize{
#'   \item{\code{models}  A \code{\link{data.table}} giving the fitted coefficients, LOOCV errors, and AICc values through time for each predictand.}
#'   \item{\code{predictands}  A vector, matrix, or RasterBrick of predictands, as provided.}
#'   \item{\code{predictor.matrix}  A matrix of predictors for calibration; \code{chronologies} cropped to \code{calibration.years}.}
#'   \item{\code{reconstruction.matrix}  A matrix of predictors for reconstruction; \code{chronologies} cropped to \code{prediction.years}, or all of \code{chronologies} if \code{prediction.years==NULL}.}
#'   \item{\code{carscores}  A list of carscore matrices, one for each predictand.}
#' }
#' @importFrom data.table setkey setnames setcolorder data.table
#' @importFrom raster ncell nlayers
#' @importFrom matrixStats rowSds
#' @importFrom broom tidy
#' @export
## YesWorkflow markup!
# @BEGIN main
# @IN ITRDB_chronologies
# @IN predictands
# @IN calibration.years
# @IN prediction.years
# @PARAM min.width
# @PARAM verbose
# @OUT final.models
paleocar_models <- function(chronologies,
                            predictands,
                            calibration.years,
                            prediction.years = NULL,
                            min.width = NULL,
                            verbose = F,
                            ...) {
  if (verbose)
    cat("Calculating PaleoCAR models\n")
  # cat("...:\n")
  # print(list(...))
  # cat("params:\n")
  # print(list(min.width = min.width,
  #            verbose = verbose))
  
  # @BEGIN get_predictor_matrix
  # @IN chronologies
  # @IN calibration.years
  # @PARAM min.width
  # @OUT predictor.matrix
  # @OUT max.preds
  t <- Sys.time()
  predictor.matrix <-
    get_predictor_matrix(
      chronologies = chronologies,
      calibration.years = calibration.years,
      min.width = min.width
    )
  # @END get_predictor_matrix
  
  maxPreds <- nrow(predictor.matrix) - 5
  
  if (class(predictands) %in% c("RasterBrick", "RasterStack")) {
    if (raster::nlayers(predictands) != length(calibration.years)) {
      stop(
        "Predictand raster must have the same number of layers as the length of the calibration.years vector!"
      )
    }
    null.cells <- which(is.na(raster::getValues(predictands[[1]])))
    predictand.matrix <- t(raster::values(predictands))
    colnames(predictand.matrix) <-
      as.character(1:raster::ncell(predictands))
  } else if (is.vector(predictands)) {
    if (length(predictands) != length(calibration.years)) {
      stop("Predictand vector must be same length as calibration.years vector!")
    }
    null.cells <- vector(mode = "integer")
    predictand.matrix <- matrix(predictands)
    colnames(predictand.matrix) <- "1"
    rownames(predictand.matrix) <- paste0("X", calibration.years)
  } else if (is.matrix(predictands)) {
    if (nrow(predictands) != length(calibration.years)) {
      stop(
        "Predictand matrix must have the same number of rows as the length of the calibration.years vector!"
      )
    }
    null.cells <- vector(mode = "integer")
    predictand.matrix <- predictands
    colnames(predictand.matrix) <- 1:ncol(predictand.matrix)
    rownames(predictand.matrix) <- paste0("X", calibration.years)
  }
  
  # @BEGIN get_reconstruction_matrix
  # @IN chronologies
  # @IN reconstruction.years
  # @PARAM min.width
  # @OUT reconstruction.matrix
  reconstruction.matrix <-
    get_reconstruction_matrix(
      chronologies = chronologies,
      reconstruction.years = prediction.years,
      min.width = min.width
    )
  reconstruction.matrix <-
    reconstruction.matrix[, colnames(predictor.matrix)]
  # @END get_reconstruction_matrix
  
  # @BEGIN get_predlist
  # @IN reconstruction.matrix
  # @OUT predlist
  predlist <- get_predlist(reconstruction.matrix)
  # @END get_predlist
  
  prednums <- rowSums(predlist, na.rm = T)
  prednums[prednums > maxPreds] <- maxPreds
  
  predyears <- as.numeric(rownames(predlist))
  hinge.year <-
    min(predyears[predyears > max(calibration.years)], max(calibration.years))
  
  # @BEGIN getCarscores
  # @IN predictand.matrix
  # @IN predictor.matrix
  # @OUT carscores
  carscores <-
    carscore_batch(predictand.matrix = predictand.matrix, predictor.matrix =
                     predictor.matrix)
  carscores.ranks <-
    matrixStats::colRanks(1 - (carscores ^ 2), preserveShape = T)
  rownames(carscores.ranks) <- rownames(carscores)
  colnames(carscores.ranks) <- colnames(carscores)
  carscores <- carscores.ranks
  rm(carscores.ranks)
  gc()
  gc()
  carscores <- data.table::data.table(t(carscores))
  
  # @END getCarscores
  
  if (verbose)
    cat("\nPrepare data and calculate CAR scores:",
        round(difftime(Sys.time(), t, units = 'mins'), digits = 2),
        "minutes\n")
  
  # @BEGIN calculateModels
  # @IN predlist
  # @IN carscores
  # @IN max.preds
  # @OUT linear.models
  allModels <-
    data.table::data.table(
      cell = numeric(),
      year = numeric(),
      model = numeric(),
      numPreds = numeric(),
      CV = numeric(),
      AICc = numeric(),
      coefs = numeric()
    )
  complete.cell.years <-
    data.table::data.table(cell = numeric(), year = numeric())
  times <- vector('numeric', max(prednums))
  # @BEGIN defineLinearModels
  # @IN predlist
  # @IN carscores
  # @OUT models
  # @OUT matches
  for (i in 1:maxPreds) {
    if (verbose)
      cat("\nCalculating models of with", i, "input vectors.")
    t <- Sys.time()
    
    new.t <- Sys.time()
    preds <- predlist[prednums >= i, ]
    
    ## MATCHES 1
    matches <- lapply(rownames(preds), function(this.year) {
      # cat(this.year,'\n')
      complete.cells <-
        c(complete.cell.years[year == as.numeric(this.year), cell], null.cells)
      
      pred.names <- colnames(preds)[preds[this.year, ]]
      if (all((rownames(carscores) %in% complete.cells)))
        return(NULL)
      models <-
        data.table::data.table(matrixStats::rowRanks(as.matrix(carscores[!(rownames(carscores) %in% complete.cells), pred.names, with = F]), ties.method = 'min') <=
                                 i)
      data.table::setnames(models, pred.names)
      # which(apply(as.matrix(models),1,function(x){!any(x)}))
      models[, cell := as.numeric(rownames(carscores)[!(rownames(carscores) %in% complete.cells)])]
      # models <- models[cell %in% which(!completed.cells)]
      data.table::setorderv(models, pred.names, rep(-1, length(pred.names)))
      models.duplicates <- !duplicated(models, by = pred.names)
      models.matches <- cumsum(models.duplicates)
      names(models.matches) <- models$cell
      models.matches <-
        models.matches[order(as.numeric(names(models.matches)))]
      models <- models[models.duplicates, ]
      models[, cell := NULL]
      blank.models <- !apply(models, 2, any)
      blank.models <- names(blank.models)[blank.models]
      if (length(blank.models) > 0) {
        models[, (blank.models) := NULL]
      }
      # if(any(apply(as.matrix(models),1,function(x){!any(x)}))) error("PROBLEM!")
      return(list(models = models, matches = models.matches))
    })
    gc()
    gc()
    
    match.names <- rownames(preds)[!sapply(matches, is.null)]
    matches <- matches[!sapply(matches, is.null)]
    
    models <- lapply(matches, '[[', 'models')
    matches <- lapply(matches, '[[', 'matches')
    names(matches) <- match.names
    names(models) <- match.names
    
    nummodels <-
      c(0, cumsum(sapply(models, nrow)))[0 - (length(models) + 1)]
    names(nummodels) <- names(models)
    
    models <- data.table::rbindlist(models, fill = T)
    pred.names <- sort(names(models))
    data.table::setcolorder(models, pred.names)
    models[, model := 1:nrow(models)]
    
    if (nrow(models) == 0)
      break
    ## Replace all NAs with FALSE (we aren't going to use those predictors, right?)
    for (this.model in seq_along(models)) {
      data.table::set(models,
                      i = which(is.na(models[[this.model]])),
                      j = this.model,
                      value = F)
    }
    
    matches <-
      lapply(1:length(matches), function(count) {
        return(list(
          year = rep(names(matches)[count], length(matches[[count]])),
          matches = matches[[count]] + nummodels[count]
        ))
      })
    
    matches.years <- unlist(lapply(matches, '[[', 'year'))
    matches.matches <- unlist(lapply(matches, '[[', 'matches'))
    
    matches <-
      do.call(rbind,
              list(
                cell = as.integer(names(matches.matches)),
                year = as.integer(matches.years),
                model = as.integer(matches.matches)
              ))
    colnames(matches) <- NULL
    matches <- data.table::data.table(t(matches))
    
    rm(matches.years, matches.matches)
    gc()
    gc()
    
    models <-
      data.table::setorderv(models, pred.names, rep(-1, length(pred.names)))
    models.duplicates <- !duplicated(models, by = pred.names)
    models.matches <- cumsum(models.duplicates)
    names(models.matches) <- models$model
    models.matches <-
      models.matches[order(as.numeric(names(models.matches)))]
    models <- models[models.duplicates, ]
    models[, model := NULL]
    models.matches <-
      data.table::data.table(model = as.integer(names(models.matches)), new.model =
                               models.matches)
    
    data.table::setkey(matches, model)
    data.table::setkey(models.matches, model)
    
    matches <- matches[models.matches]
    matches[, model := NULL]
    data.table::setnames(matches, c('cell', 'year', 'model'))
    matches <-
      data.table::setorderv(matches, c('cell', 'year', 'model'), c(1, 1, 1))
    
    rm(models.duplicates, models.matches)
    gc()
    gc()
    
    if (nrow(complete.cell.years) > 0) {
      which.matches <-
        merge(
          matches,
          complete.cell.years,
          by = c("cell", "year"),
          all.x = T
        )[is.na(model.y)]
      matches <- which.matches[, list(cell, year, model.x)]
      setnames(matches, c("cell", 'year', 'model'))
    }
    
    data.table::setkey(matches, model)
    # @END defineLinearModels
    if (verbose)
      cat("\nDefine models:", round(difftime(Sys.time(), new.t, units = 'mins'), digits =
                                      2), "minutes")
    # @BEGIN calculateLinearModels
    # @IN models
    # @IN matches
    # @OUT coefficients
    # @OUT model.metrics
    new.t <- Sys.time()
    all.lms <- lapply(1:nrow(models), function(this.model) {
      # cat(this.model,'\n')
      if (!(this.model %in% matches[['model']]))
        return(NULL)
      cells <- unique(matches[list(this.model), cell])
      
      terms.names <-
        names(models[this.model])[which(as.logical(models[this.model]))]
      terms <- predictor.matrix[, terms.names, drop = F] %>%
        as.data.frame()
      response <- predictand.matrix[, cells, drop = F]
      
      model.mlm <- lm(response ~ .,
                      data = terms)
      
      # Get model cross-validation error, AICc, AdjR2
      model.metrics <- CV_mlm(model.mlm) %>%
        tibble::as_tibble(rownames = "cell") %>%
        dplyr::mutate(model = this.model)
      
      if ("mlm" %in% class(model.mlm)) {
        model.mlm %<>%
          unlist_mlm() %>%
          purrr::map(strip_lm) %>%
          tibble::tibble(cell = names(.),
                         coefs = .)
      } else {
        model.mlm %<>%
          strip_lm() %>%
          list() %>%
          magrittr::set_names(colnames(response)) %>%
          tibble::tibble(cell = names(.),
                         coefs = .)
      }
      
      model.metrics %<>%
        dplyr::full_join(model.mlm,
                         by = "cell") %>%
        dplyr::select(cell, model, dplyr::everything()) %>%
        dplyr::mutate(cell = as.integer(cell)) %>%
        dplyr::arrange(cell, model) %>%
        data.table::data.table()
      
      data.table::setkey(model.metrics, cell, model)
      
      return(model.metrics)
    })
    all.lms <- all.lms[!sapply(all.lms, is.null)]
    if (verbose)
      cat(
        "\nCalculate",
        length(all.lms),
        "linear models:",
        round(difftime(Sys.time(),
                       new.t,
                       units = 'mins'),
              digits = 2),
        "minutes"
      )
    rm(models)
    gc()
    gc()
    # @END calculateLinearModels
    
    # @BEGIN simplifyLinearModels
    # @IN coefficients
    # @IN model.metrics
    # @OUT final.models
    ## LMS SIMPLIFY PREP
    new.t <- Sys.time()
    
    all.lms <- data.table::rbindlist(all.lms, fill = T)
    setkey(all.lms, cell, model)
    all.lms[, numPreds := i]
    
    if (nrow(allModels) == 0) {
      data.table::setkey(all.lms, cell, model)
      data.table::setkey(matches, cell, model)
      allModels <- matches[all.lms]
      
      setcolorder(allModels,
                  c("cell", "year", "model", "numPreds", "CV", "AICc", "AdjR2", "coefs"))
    }
    
    all.lms.stats <- all.lms[, list(cell, model, AICc)]
    setkey(all.lms.stats, cell, model)
    setkey(matches, cell, model)
    
    all.lms.stats <- matches[all.lms.stats]
    rm(matches)
    
    data.table::setkey(all.lms.stats, cell, year, model)
    
    allModels.stats <- allModels[, list(cell, year, model, AICc)]
    
    allModels.stats <- rbind(allModels.stats, all.lms.stats, fill = T)
    
    data.table::setkey(allModels.stats, cell, year)
    data.table::setorder(allModels.stats, cell, year, AICc)
    allModels.stats <-
      allModels.stats[allModels.stats[, !duplicated(year), by = cell]$V1, ]
    
    data.table::setkey(allModels, cell, year, model)
    if (nrow(allModels.stats[model == 0, list(cell, year, model)]) > 0) {
      allModels.old <-
        merge(
          allModels.stats[model == 0, list(cell, year, model)],
          allModels,
          by = c("cell", "year", "model"),
          all = T
        )
    } else{
      allModels.old <- data.table(
        cell = numeric(),
        year = numeric(),
        model = numeric(),
        numPreds = numeric(),
        CV = numeric(),
        AICc = numeric(),
        AdjR2 = numeric(),
        coefs = list()
      )
    }
    
    if (nrow(allModels.stats[model != 0, list(cell, year, model)]) > 0) {
      allModels.new <-
        merge(
          allModels.stats[model != 0, list(cell, year, model)],
          all.lms.stats,
          by = c("cell", "year", "model"),
          all = T
        )
      allModels.new <-
        merge(allModels.new[, list(cell, year, model)],
              all.lms,
              by = c("cell", "model"),
              all = T)
      data.table::setcolorder(allModels.new,
                              c("cell", "year", "model", "numPreds", "CV", "AICc", "AdjR2", "coefs"))
    } else{
      allModels.new <-
        data.table::data.table(
          cell = numeric(),
          year = numeric(),
          model = numeric(),
          numPreds = numeric(),
          CV = numeric(),
          AICc = numeric(),
          AdjR2 = numeric(),
          coefs = list()
        )
    }
    
    allModels <- rbind(allModels.old, allModels.new)
    
    data.table::setkey(allModels, cell, year)
    data.table::setorder(allModels, cell, year, AICc)
    allModels <-
      allModels[allModels[, !duplicated(year), by = cell]$V1, ]
    
    complete.cell.years <-
      allModels[model == 0, list(cell, year, model)]
    
    allModels[, model := 0]
    
    setkey(allModels, cell, year, numPreds)
    
    if (verbose)
      cat("\nClean linear models:",
          round(difftime(Sys.time(), new.t, units = 'mins'), digits = 2),
          "minutes")
    # @END simplifyLinearModels
    
    time <- difftime(Sys.time(), t, units = 'mins')
    if (verbose)
      cat("\nTotal modeling time:",
          round(time, digits = 2),
          "minutes\n")
    times[i] <- time
    if ((nrow(allModels) - nrow(complete.cell.years)) == 0)
      break
    if (verbose)
      cat(nrow(allModels) - nrow(complete.cell.years),
          "cell-years remaining\n")
    ##
    
  }
  
  if (verbose)
    cat("\nTotal Modeling Time:", sum(times), "minutes\n")
  # @END calculateModels
  
  # @BEGIN optimizeModels
  # @IN linear.models
  # @OUT final.models
  t <- Sys.time()
  get.coef.names <- function(year, model, coefs, numPreds, CV, AICc, AdjR2) {
    the.coefs <- lapply(coefs, 
                        function(x) {
      names(stats::coefficients(x))[-1]
    })
    test.out <-
      lapply(the.coefs, function(x) {
        which(rowSums(predlist[, x, drop = F]) == length(x))
      })
    test.out <-
      data.table::data.table(model = rep(1:length(test.out), 
                                         times = sapply(test.out, length)),
                             year = unlist(test.out))
    # test.out <- lapply(1:length(test.out),function(i){data.table(model=i,year=test.out[[i]])})
    # test.out <- rbindlist(test.out)
    test.out <- test.out[which(!duplicated(test.out[, year]))]
    setkey(test.out, year)
    return(
      list(
        year = sort(year),
        model = model[test.out$model],
        numPreds = numPreds[test.out$model],
        CV = CV[test.out$model],
        AICc = AICc[test.out$model],
        AdjR2 = AdjR2[test.out$model],
        coefs = coefs[test.out$model]
      )
    )
  }
  
  allModels <-
    allModels[order(AICc)][, get.coef.names(year, model, coefs, numPreds, CV, AICc, AdjR2), 
                           by = cell]
  setkey(allModels, cell, year)
  
  allModels <-
    allModels[allModels[, !FedData::sequential_duplicated(CV), by = cell]$V1, ] %>%
    dplyr::select(-model) %>%
    dplyr::rename(model = coefs)
  
  time <- difftime(Sys.time(), t, units = 'mins')
  if (verbose)
    cat("\nOptimizing models:", round(time, digits = 2), "minutes\n")
  
  allModels <- list(
    models = allModels %>% 
      tibble::as_tibble(),
    predictands = predictands,
    predictor.matrix = predictor.matrix,
    reconstruction.matrix = reconstruction.matrix,
    carscores = carscores %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(., cell = 1:nrow(.)) %>%
      tidyr::gather("series", "car_score", -cell)
  )
  
  # @END optimizeModels
  return(allModels)
  # @END main
}
