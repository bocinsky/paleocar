#' Tidy function for mlm objects including confidence intervals
#'
#' @param model An object of class mlm.
#' @param conf.int Whether to calculate confidence intervals.
#' @param level Level for confidence interval
#' @return A tibble of tidy model results.
#' @importFrom broom tidy
tidy_mlm <- function (model, conf.int = FALSE, level = 0.95) {
  if(!("mlm" %in% class(model)))
    return(broom::tidy(model, conf.int = conf.int, level = level) %>%
             dplyr::mutate(response = colnames((model$model)[[1]])))
  
  if(!conf.int)
    return(broom::tidy(model))
  
  beta <- coef(model)
  
  Rinv <- with(model$qr, backsolve(qr, diag(rank)))
  ## unscaled standard error
  std_unscaled <- sqrt(rowSums(Rinv ^ 2)[order(model$qr$pivot)])
  ## residual standard error
  sigma <- sqrt(colSums(model$residuals ^ 2) / model$df.residual)
  ## return final standard error
  ## each column corresponds to a model
  se <- "dimnames<-"(outer(std_unscaled, sigma), list = dimnames(model$coefficients))
  
  alpha <- qt((1 - level) / 2, df = model$df.residual)
  
  out <- list(lower = beta + alpha * se, 
              upper = beta - alpha * se) %>%
    purrr::map(tibble::as_tibble, rownames = "term") %>%
    # purrr::map(broom::tidy) %>%
    dplyr::bind_rows(.id = "level") %>%
    # dplyr::rename(term = `.rownames`) %>%
    tidyr::gather("response","estimate",-level,-term) %>%
    tidyr::spread("level","estimate") %>%
    dplyr::arrange(response,term) %>%
    dplyr::select(response, term, lower, upper) %>%
    dplyr::rename(conf.low = lower,
                  conf.high = upper)
  
  return(dplyr::full_join(broom::tidy(model),
                   out,
                   by = c("response","term")))
}