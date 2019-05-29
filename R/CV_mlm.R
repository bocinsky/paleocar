#' Computation of LOOCV and AICc for mlm objects
#'
#' Extends the CV method from the forecast package to handle 
#' multi-predictand linear models (objects of class mlm). Specifically, it calculates
#' the leave-one-out cross validataion statistic (also known as the PRESS statistic) and
#' the corrected Akaike's Information Criterion.
#'
#' @param obj An object of class mlm or lm.
#' @return A matrix of of CV and AICc values for each model.
#' @importFrom stats hatvalues
#' @importFrom forecast CV
CV_mlm <- function (obj) 
{
  if(class(obj)[[1]]=="lm"){
    return(forecast::CV(obj)[c("CV", "AICc", "AdjR2")] %>%
      as.matrix() %>%
      t() %>%
      magrittr::set_rownames(colnames((obj$model)[[1]])))
  }
  n <- nrow(obj$residuals)
  aic.raw <- AIC_mlm(obj)
  k <- aic.raw$df - 1
  aic <- aic.raw$AIC + 2
  rm(aic.raw)
  aicc <- aic + 2 * (k + 2) * (k + 3)/(n - k - 3)
#   bic <- aic + (k + 2) * (log(n) - 2)
  cv <- colMeans((obj$residuals/(1 - stats::hatvalues(obj)))^2, na.rm = TRUE)
  adjr2 <- summary(obj) %>%
    purrr::map(`[[`,"adj.r.squared") %>%
    unlist()
  
  out <- do.call(cbind,list(cv, aicc, adjr2))
  colnames(out) <- c("CV", "AICc", "AdjR2")
  return(out)
}