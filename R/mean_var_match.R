mean_var_match <- function(calib.vector, out.calib.vector, out.vector){
  scalar <- sd(calib.vector)/sd(out.calib.vector)
  transform <- mean(calib.vector)-(scalar*mean(out.calib.vector))
  out <- (out.vector*scalar) + transform
  return(out)
}