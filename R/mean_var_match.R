#' @importFrom stats sd
mean_var_match <- function(calib.vector, out.calib.vector, out.vector){
  scalar <- stats::sd(calib.vector, na.rm  = T)/stats::sd(out.calib.vector, na.rm = T)
  transform <- mean(calib.vector, na.rm = T)-(scalar*mean(out.calib.vector, na.rm = T))
  out <- (out.vector*scalar) + transform
  return(out)
}