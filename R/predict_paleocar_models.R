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
  
  my_predict <- function(start, end, model){
    predict(model, 
            newdata = models$reconstruction.matrix[rownames(models$reconstruction.matrix) %in% 
                                                     (seq(start, end) %>% 
                                                        as.character()),,drop = FALSE] %>%
              as.data.frame(),
            interval = "confidence")[, 1:2, drop = FALSE] %>%
      tibble::as_tibble() %>%
      cbind(predict(model,
                    newdata = models$reconstruction.matrix[rownames(models$reconstruction.matrix) %in%
                                                             (seq(start, end) %>%
                                                                as.character()),,drop = FALSE] %>%
                      as.data.frame(),
                    interval = "prediction")[, 2, drop = FALSE] %>%
              tibble::as_tibble()) %>%
      magrittr::set_colnames(c(
        "fit", 
        "ci.lwr", 
        "pi.lwr"
      )) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(year = seq(start, end),
                    ci = fit - ci.lwr,
                    pi = fit - pi.lwr) %>%
      dplyr::select(year, fit, ci, pi) %>%
      dplyr::rename(Prediction = fit,
                    `CI Deviation` = ci,
                    `PI Deviation` = pi)
    
  }
  
  out <- 
    models$models %>%
    tibble::as_tibble() %>%
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
      )) %>%
    dplyr::filter(!(year > tail(prediction.years, 1)),
                  !(endYear < head(prediction.years, 1))) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(prediction = my_predict(start = year, 
                                          end = endYear, 
                                          model = model) %>% 
                    list()) %>%
    dplyr::ungroup() %>%
    dplyr::select(cell,
                  prediction) %>%
    tidyr::unnest(prediction)
  
  
  if (class(models$predictands) %in% c("RasterBrick", "RasterStack")) {
    out %<>%
      dplyr::arrange(year, cell) %>%
      dplyr::select(-cell) %>%
      tidyr::nest(-year) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(data = data %>% 
                      
                      purrr::map(.f = function(x){
                        raster::setValues(raster::raster(models$predictands), 
                                          values = x)
                      }) %>%
                      list()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(data = magrittr::set_names(data, year)) %$%
      data %>%
      purrr::transpose() %>%
      purrr::map(raster::brick)
    
  }
  
  return(out)
}
