#' Make a vector strictly monotonic
#'
#' Given a vector \code{v} of length greater than 1, returns a logical vector
#' the same length as \code{v} indicating whether each element of \code{v} is
#' greater (or less-than, if \code{decreasing==TRUE}) than all previous elements.
#'
#' @param vect A numeric vector.
#' @param decreasing Should the output indicate monotonically decreasing?
#' @return A logical vector.
makeMonotonic <- function(vect,decreasing=T){
  if(length(vect)==0) return(NULL)
  if(length(vect)==1) return(T)
  if(decreasing){
    keep <- c(T,sapply(2:length(vect),function(i){
      ifelse(all(vect[i]<vect[1:i-1]),T,F)
    }))
  }else{
    keep <- c(T,sapply(2:length(vect),function(i){
      ifelse(all(vect[i]>vect[1:i-1]),T,F)
    }))
  }
  
  return(keep)
}