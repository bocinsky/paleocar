predict_mlm <- 
  function (object, 
            newdata, 
            na.action = na.pass,
            level = 0.95,
            ...) 
  {
    
    # if (missing(newdata))
    #   stop("the 'se.fit' argument is not yet implemented for \"mlm\" objects")
    if (missing(newdata)) {
      newdata <- as.data.frame(object$model$terms)
      X <- model.matrix(object)
      offset <- object$offset
    } else {
      tt <- terms(object)
      Terms <- delete.response(tt)
      m <- model.frame(Terms, 
                       newdata, 
                       na.action = na.action, 
                       xlev = object$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses"))) 
        .checkMFClasses(cl, m)
      X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
      offset <- if (!is.null(off.num <- attr(tt, "offset"))) 
        eval(attr(tt, "variables")[[off.num + 1]], newdata)
      else if (!is.null(object$offset)) 
        eval(object$call$offset, newdata)
    }
    
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq_len(p) 
    
    piv <- base::qr(object)$pivot[p1]
    # piv <- stats:::qr.lm(object)$pivot[p1]
    pred <- X[, piv, drop = FALSE] %*% object$coefficients[piv, ]
    
    ## model formula
    form <- formula(object)
    ## drop response (LHS)
    form[[2]] <- NULL
    ## prediction matrix
    X <- model.matrix(form, newdata)
    Q <- forwardsolve(t(qr.R(object$qr)), t(X))
    
    ## unscaled prediction standard error
    unscaled.se <- sqrt(colSums(Q ^ 2))
    
    # residual sum of squared errors
    rss <- colSums(residuals(object) ^ 2)
    
    # mean squared error
    mse <- rss / (object$df.residual + 2)
    
    ## residual standard error
    sigma <- sqrt(rss / object$df.residual)
    
    ## scaled standard error of the mean
    se.mean <- tcrossprod(unscaled.se, sigma)
    
    ## scaled standard error of the prediction
    se.pred = sqrt(sweep((se.mean ^ 2), 2, (sigma ^ 2), "+"))
    
    tfrac <- abs(qt((1 - level)/2, object$df.residual))
    
    ci_error <- se.mean * tfrac
    pi_error <- se.pred * tfrac
    
    if (!is.null(offset)) 
      pred <- pred + offset
    
    ci_error <- ci_error[,!duplicated(colnames(pred))]
    pi_error <- pi_error[,!duplicated(colnames(pred))]
    pred <- pred[,!duplicated(colnames(pred))]
    
    return(
      pred %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "year") %>%
        tidyr::gather(cell,Prediction,-year) %>%
        dplyr::bind_cols(    
          tibble::tibble(
            `CI Deviation` = 
              ci_error %>% 
              as.vector(),
            `PI Deviation` = 
              pi_error %>% 
              as.vector()
          )
        ) %>%
        tibble::as_tibble() %>%
        dplyr::mutate_at(.vars = dplyr::vars(cell, year), as.integer) %>%
        dplyr::select(cell, year, dplyr::everything()) %>%
        dplyr::arrange(cell, year)
    ) 
  }

