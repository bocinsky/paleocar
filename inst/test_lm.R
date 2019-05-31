test_lm <- lm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris)

head(test_lm$qr$qr)
head(qr(as.matrix(iris[,c("Sepal.Width","Petal.Length")]))$qr)
head(qr(iris[,c("Sepal.Length","Sepal.Width","Petal.Length")])$qr)
head(qr(iris[,c("Sepal.Length","Sepal.Width","Petal.Length")])$qr)


strip_lm <- function (object) 
{
  op <- object
  op$y <- NULL
  op$model <- NULL
  # op$residuals <- NULL
  names(op$residuals) <- NULL
  op$qr$qr <- NULL
  attr(op$qr$qr,"dimnames") <- NULL
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

test_lm_strip <- strip::strip(test_lm, keep = "predict")
test_lm_slim <- strip_lm(test_lm)

list(test_lm = test_lm,
     test_lm_slim = test_lm_slim,
     test_lm_strip = test_lm_strip) %>%
  purrr::map_dfr(function(x){
    x %>% 
      purrr::map(object.size) %>%
      tibble::as_tibble()
  },
  .id = "Type")

identical(predict(test_lm, newdata = iris),
          predict(test_lm_slim, newdata = iris))

identical(predict(test_lm, newdata = iris, interval = "confidence"),
          predict(test_lm_slim, newdata = iris, interval = "confidence"))

identical(predict(test_lm, newdata = iris, interval = "prediction"),
          predict(test_lm_slim, newdata = iris, interval = "prediction"))

identical(predict(test_lm, newdata = iris),
          predict(test_lm_strip, newdata = iris))

identical(predict(test_lm, newdata = iris, interval = "confidence"),
          predict(test_lm_strip, newdata = iris, interval = "confidence"))

identical(predict(test_lm, newdata = iris, interval = "prediction"),
          predict(test_lm_strip, newdata = iris, interval = "prediction"))



