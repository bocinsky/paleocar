#' Computation of AIC for mlm objects
#'
#' Extends the \code{extractAIC} method from the \pkg{stats} package to handle 
#' multi-predictand linear models (objects of class mlm).
#'
#' @param fit An object of class mlm.
#' @param scale The estimate of the error variance.
#' \code{scale = 0} indicates that it is to be 
#' estimated by maximum likelihood.
#' @param k Numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#' @return A list of length 2 giving
#' \itemize{
#'   \item{\code{df}  The ‘equivalent degrees of freedom’ for the fitted model \code{fit}.}
#'   \item{\code{AIC}  A vector of the (generalized) Akaike Information Criterion for the fits.}
#' }
mlmAIC <- function(fit, scale = 0, k = 2, ...){
  n <- nrow(fit$residuals)
  edf <- n - fit$df.residual
  RSS <- deviance(fit)
  dev <- if (scale > 0) 
    RSS/scale - n
  else n * log(RSS/n)
  list(df=edf, AIC=dev + k * edf)
}