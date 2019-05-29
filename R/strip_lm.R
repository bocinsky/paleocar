strip_lm <- function (object) 
{
  op <- object
  op$y <- NULL
  op$model <- NULL
  # op$residuals <- NULL
  op$fitted.values <- NULL
  op$effects <- NULL
  op$weights <- NULL
  op$prior.weights <- NULL
  op$linear.predictors <- NULL
  attr(op$terms, ".Environment") <- NULL
  attr(op$formula, ".Environment") <- NULL

  class(op) <- class(object)

  return(op)
}
