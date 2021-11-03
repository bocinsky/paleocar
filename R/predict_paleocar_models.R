#' Get a RasterBrick reconstruction from a PaleoCAR batch model
#'
#' This generates a reconstruction RasterBrick from a PaleoCAR batch model.
#'
#' @param models A PaleoCAR batch model, as returned from \code{\link{paleocar_models}}.
#' @param prediction.years The set of years over which to generate reconstruction rasters. Optional.
#' @param ... Further arguments to be passed to other functions.
#' @return A tibble or RasterBrick containing the predictions for each year.
#' @export
predict_paleocar_models <- function(models,
                                    prediction.years = NULL,
                                    ...) {
  
  if (is.null(prediction.years))
    prediction.years <-
      as.numeric(rownames(models$reconstruction.matrix))
  
  if (!all(prediction.years %in% as.numeric(row.names(models$reconstruction.matrix)))) {
    warning('Some specified years not valid for the models provided. Truncating to modeled years.')
    prediction.years <-
      prediction.years[prediction.years %in% 
                         as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  if (class(models$predictands)[[1]] %in% c("RasterBrick", "RasterStack")) {
    null.cells <- which(is.na(raster::getValues(models$predictands[[1]])))
    predictand.matrix <- t(raster::values(models$predictands))
    colnames(predictand.matrix) <-
      as.character(1:raster::ncell(models$predictands))
  } else if (is.vector(models$predictands)) {
    null.cells <- vector(mode = "integer")
    predictand.matrix <- matrix(models$predictands)
    rownames(predictand.matrix) <- names(models$predictands)
    colnames(predictand.matrix) <- 1
  } else if (is.matrix(models$predictands)) {
    null.cells <- vector(mode = "integer")
    predictand.matrix <- models$predictands
    colnames(predictand.matrix) <- 1:ncol(predictand.matrix)
  }
  
  models$models %<>%
    tibble::as_tibble() %>%
    dplyr::select(cell, year, model) %>%
    dplyr::group_by(cell) %>%
    dplyr::mutate(
      endYear = c(year[-1] - 1, 
                  tail(prediction.years, 1)),
      year = ifelse(
        year < head(prediction.years, 1),
        head(prediction.years, 1),
        year
      ),
      endYear = ifelse(
        endYear > tail(prediction.years, 1),
        tail(prediction.years, 1),
        endYear
      ),
      endYear = as.integer(endYear)) %>%
    dplyr::filter(!(year > tail(prediction.years, 1)),
                  !(endYear < head(prediction.years, 1))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cell, year)
  
  my_predict <- function(x){
    
    terms <- models$predictor.matrix[, x$model[[1]], drop = F] %>%
      as.data.frame()
    response <- predictand.matrix[, unique(x$cell), drop = F]
    
    lms <- lm(response ~ .,
              data = terms)
    
    if(!("mlm" %in% class(lms))){
      
      x %>%
        dplyr::select(cell, year, endYear) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(year = list(year:endYear)) %>%
        dplyr::select(cell, year) %>%
        tidyr::unnest(cols = c(year)) %>%
        dplyr::mutate_at(.vars = dplyr::vars(cell, year), as.integer) %>%
        dplyr::left_join(
          predict(lms,
                  newdata = 
                    models$reconstruction.matrix[, x$model[[1]], drop = F] %>% 
                    as.data.frame() %>%
                    na.omit(),
                  interval = "prediction")[, 1:2, drop = FALSE] %>%
            tibble::as_tibble(rownames = "year") %>%
            dplyr::mutate(year = as.integer(year),
                          pi = fit - lwr) %>%
            dplyr::select(year,
                          Prediction = fit,
                          `PI Deviation` = pi),
          by = c("year")
        ) %>%
        dplyr::mutate(
          `Prediction (scaled)` = 
            mean_var_match(calib.vector = lms$model$response[,1], 
                           out.calib.vector = lms$fitted.values, 
                           out.vector = Prediction),
          `PI Deviation (scaled)` = 
            `PI Deviation` * 
            (stats::sd(lms$model$response[,1], na.rm  = T)/
               stats::sd(lms$fitted.values, na.rm = T))
        )
    } else {
      scalar <- 
        apply(lms$model$response, FUN = sd, MARGIN = 2)/
        apply(lms$fitted.values, FUN = sd, MARGIN = 2)
      
      transform <- 
        apply(lms$model$response, FUN = mean, MARGIN = 2) -
        (scalar * apply(lms$fitted.values, FUN = mean, MARGIN = 2))
      
      x %>%
        dplyr::select(cell, year, endYear) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(year = list(year:endYear)) %>%
        dplyr::select(cell, year) %>%
        tidyr::unnest(cols = c(year)) %>%
        dplyr::arrange(cell, year) %>%
        dplyr::left_join(
          predict_mlm(object = lms, 
                      newdata = 
                        models$reconstruction.matrix[, x$model[[1]], drop = F] %>% 
                        as.data.frame() %>%
                        na.omit()
          ),
          by = c("cell","year")) %>%
        # dplyr::mutate(cell = as.character(cell)) %>%
        dplyr::left_join(tibble::tibble(cell = as.integer(names(scalar)), scalar = scalar) %>% dplyr::distinct(),
                         by = "cell") %>%
        dplyr::left_join(tibble::tibble(cell = as.integer(names(transform)), transform = transform),
                         by = "cell") %>%
        dplyr::mutate(`Prediction (scaled)` = (Prediction * scalar) + transform,
                      `PI Deviation (scaled)` = `PI Deviation` * scalar,
                      cell = as.integer(cell)) %>%
        dplyr::select(-scalar,
                      -transform) %>%
        dplyr::arrange(cell, year)
      
    }
    
  }
  
  out <- 
    models$models %>%
    dplyr::arrange(year) %>%
    dplyr::group_by(model) %>%
    dplyr::group_split() %>%
    purrr::map(my_predict) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(cell, year)
  
  
  if (class(models$predictands)[[1]] %in% c("RasterBrick", "RasterStack")) {
    out %<>% 
      dplyr::arrange(year, cell) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(dplyr::across(!c(cell), 
                                     ~ list(
                                       raster::setValues(raster::raster(models$predictands), 
                                                         values = .x, 
                                                         index = cell)
                                     ))) %>%
      dplyr::mutate(dplyr::across(!c(year), ~magrittr::set_names(.x, year))) %>%
      dplyr::summarise(dplyr::across(!c(year), ~list(raster::brick(.x)))) %>%
      unlist(recursive = FALSE) %>% 
      purrr::map(raster::readAll)
  }
  
  return(out)
}
