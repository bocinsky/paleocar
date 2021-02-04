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
      prediction.years[prediction.years %in% as.numeric(rownames(models[['reconstruction.matrix']]))]
  }
  
  if (class(models$predictands) %in% c("RasterBrick", "RasterStack")) {
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
    
    # x <- models$models %>%
    #   split(models$models$model %>% purrr::map_chr(stringr::str_c, collapse = ";")) %>%
    #   magrittr::extract2(c("CO607RA;CO614RA;CO621RA;CO622RA;CO625RA;CO626RA;CO637RA"))
    
    terms <- models$predictor.matrix[, x$model[[1]], drop = F] %>%
      as.data.frame()
    response <- predictand.matrix[, x$cell, drop = F]
    
    lms <- lm(response ~ .,
              data = terms)
    
    if(!("mlm" %in% class(lms))){
      
      x %>%
        dplyr::select(cell, year, endYear) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(year = list(year:endYear)) %>%
        dplyr::select(cell, year) %>%
        tidyr::unnest() %>%
        dplyr::mutate_at(.vars = dplyr::vars(cell, year), as.integer) %>%
        dplyr::left_join(predict(lms, 
                                 newdata = models$reconstruction.matrix[, x$model[[1]], drop = F] %>% 
                                   as.data.frame() %>%
                                   na.omit(), interval = "confidence")[, 1:2, drop = FALSE] %>%
                           as.data.frame() %>%
                           tibble::rownames_to_column("year") %>%
                           dplyr::full_join(predict(lms,
                                                    newdata = models$reconstruction.matrix[, x$model[[1]], drop = F] %>% 
                                                      as.data.frame() %>%
                                                      na.omit(),
                                                    interval = "prediction")[, 2, drop = FALSE] %>%
                                              as.data.frame() %>%
                                              tibble::rownames_to_column("year"),
                                            by = "year") %>%
                           tibble::as_tibble() %>%
                           magrittr::set_names(c(
                             "year",
                             "fit",
                             "ci.lwr",
                             "pi.lwr"
                           )) %>%
                           dplyr::mutate(year = as.integer(year),
                                         ci = fit - ci.lwr,
                                         pi = fit - pi.lwr) %>%
                           dplyr::select(year, fit, ci, pi) %>%
                           dplyr::rename(Prediction = fit,
                                         `CI Deviation` = ci,
                                         `PI Deviation` = pi),
                         by = c("year"))
    } else {
      x %>%
        dplyr::select(cell, year, endYear) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(year = list(year:endYear)) %>%
        dplyr::select(cell, year) %>%
        tidyr::unnest() %>%
        dplyr::arrange(cell, year) %>%
        dplyr::left_join(predict_mlm(object = lms, 
                                     newdata = models$reconstruction.matrix[, x$model[[1]], drop = F] %>% 
                                       as.data.frame() %>%
                                       na.omit()),
                         by = c("cell","year"))
    }
    
  }
  
  out <- 
    models$models %>%
    split(models$models$model %>% 
            purrr::map_chr(stringr::str_c, collapse = ";")) %>%
    purrr::map(my_predict) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(cell, year)
  
  
  
  # 
  # my_predict <- function(start, end, model){
  #   
  #   predict(model, 
  #           newdata = models$reconstruction.matrix[rownames(models$reconstruction.matrix) %in% 
  #                                                    (seq(start, end) %>% 
  #                                                       as.character()),,drop = FALSE] %>%
  #             as.data.frame(),
  #           interval = "confidence")[, 1:2, drop = FALSE] %>%
  #     tibble::as_tibble() %>%
  #     cbind(predict(model,
  #                   newdata = models$reconstruction.matrix[rownames(models$reconstruction.matrix) %in%
  #                                                            (seq(start, end) %>%
  #                                                               as.character()),,drop = FALSE] %>%
  #                     as.data.frame(),
  #                   interval = "prediction")[, 2, drop = FALSE] %>%
  #             tibble::as_tibble()) %>%
  #     magrittr::set_colnames(c(
  #       "fit", 
  #       "ci.lwr", 
  #       "pi.lwr"
  #     )) %>%
  #     tibble::as_tibble() %>%
  #     dplyr::mutate(year = seq(start, end),
  #                   ci = fit - ci.lwr,
  #                   pi = fit - pi.lwr) %>%
  #     dplyr::select(year, fit, ci, pi) %>%
  #     dplyr::rename(Prediction = fit,
  #                   `CI Deviation` = ci,
  #                   `PI Deviation` = pi)
  #   
  # }
  # 
  # out <- 
  #   models$models %>%
  #   tibble::as_tibble() %>%
  #   dplyr::group_by(cell) %>%
  #   dplyr::mutate(
  #     endYear = c(year[-1] - 1, 
  #                 tail(prediction.years, 1)),
  #     year = ifelse(
  #       year < head(prediction.years, 1),
  #       head(prediction.years, 1),
  #       year
  #     ),
  #     endYear = ifelse(
  #       endYear > tail(prediction.years, 1),
  #       tail(prediction.years, 1),
  #       endYear
  #     )) %>%
  #   dplyr::filter(!(year > tail(prediction.years, 1)),
  #                 !(endYear < head(prediction.years, 1))) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::rowwise() %>%
  #   dplyr::mutate(prediction = my_predict(start = year, 
  #                                         end = endYear, 
  #                                         model = model) %>% 
  #                   list()) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select(cell,
  #                 prediction) %>%
  #   tidyr::unnest(prediction)
  # 
  
  if (class(models$predictands) %in% c("RasterBrick", "RasterStack")) {
    out %<>% 
      dplyr::arrange(year, cell) %>%
      tidyr::nest(data = c(cell, 
                           Prediction, 
                           `CI Deviation`, 
                           `PI Deviation`)) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(data = 
                      c("Prediction",
                        "CI Deviation",
                        "PI Deviation") %>%
                      magrittr::set_names(.,.) %>%
                      purrr::map(function(x){
                        out_rast <- raster::raster(models$predictands)
                        out_rast[data[[x]]] <- data$cell
                        out_rast
                      }) %>% 
                      list()) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(data = magrittr::set_names(data, year)) %$% 
      data %>% 
      purrr::transpose() %>% 
      purrr::map(raster::brick) %>% 
      purrr::map(raster::readAll)
  }
  
  # out %>%
  #   purrr::map(raster::mean) %>%
  #   purrr::walk(raster::plot)
  
  return(out)
}
