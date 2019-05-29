unlist_mlm <- function(object){
  nmodels <- ncol(object$coefficients)
  out <- list()
  
  # Coefficients
  out$coefficients <- lapply(seq_len(nmodels), 
                             function(i) object$coefficients[, i])
  names(out$coefficients) <- colnames(object$coefficients)
  
  # names(object$residuals) <- rownames(object$model$response)
  out$residuals <- lapply(seq_len(nmodels), 
                          function(i) object$residuals[, i])
  names(out$residuals) <- colnames(object$residuals)
  
  out$effects <- lapply(seq_len(nmodels), 
                        function(i) object$effects[, i])
  names(out$effects) <- colnames(object$effects)
  
  out$rank <- rep(list(object$rank), nmodels)
  
  # names(object$fitted.values) <- rownames(object$model$response)
  out$fitted.values <- lapply(seq_len(nmodels), 
                              function(i) object$fitted.values[, i])
  names(out$fitted.values) <- colnames(object$fitted.values)
  
  out$assign <- rep(list(object$assign), nmodels)
  
  out$qr <- rep(list(object$qr), nmodels)
  out$df.residual <- rep(list(object$df.residual), nmodels)
  out$xlevels <- rep(list(object$xlevels), nmodels)
  out$call <- rep(list(object$call), nmodels)
  out$terms <- object$terms
  attr(out$terms, "dataClasses")["response"] <- "numeric"
  out$terms <- rep(list(out$terms), nmodels)
  
  out$model <- object$model %>% 
    as.list()
  out$model$response <- lapply(seq_len(nmodels), 
                               function(i) {
                                 object$model$response[, i, drop = FALSE]
                                 # dimnames(out) <- list(names(out),colnames(object$model$response)[[i]])
                               })
  # names(out$model$response) <- rep("response", nmodels)
  
  
  setdiff(names(out$model),"response") %>%
    magrittr::set_names(setdiff(names(out$model),"response")) %>%
    purrr::walk(function(x){
      out$model[[x]] <<- rep(list(object$model[[x]]), nmodels)
    })

  out$model %<>%
    purrr::transpose() %>%
    purrr::map(function(x){
      the.dims <- dimnames(x$response)
      x %<>%
        as.data.frame()
      names(x)[[1]] <- "response"
      attr(x,"terms") <- attr(object$model,"terms")
      attr(attr(x,"terms"), "dataClasses")["response"] <- "numeric"
      dim(x$response) <- c(nrow(x),1)
      attr(x$response, "dimnames") <- the.dims
      
      return(x)
    })
   
  out %<>%
    purrr::transpose() %>%
    purrr::map(function(x){
      class(x) <- "lm"
      return(x)
    })
  
  return(out)
}
