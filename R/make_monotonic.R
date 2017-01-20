#' Make a vector strictly monotonic
#'
#' Given a vector \code{v} of length greater than 1, returns a logical vector
#' the same length as \code{v} indicating whether each element of \code{v} is
#' greater (or less-than, if \code{decreasing==TRUE}) than all previous elements.
#'
#' @param vect A numeric vector.
#' @param decreasing Should the output indicate monotonically decreasing?
#' @return A logical vector.
make_monotonic <- function(vect,decreasing=T){
  if(length(vect)==0) return(NULL)
  if(length(vect)==1) return(T)
  if(decreasing){
    keep <- c(T,sapply(2:length(vect),function(i){
      !any(sapply(vect[1:i-1],all.equal,current = vect[i]) == "TRUE") & all(vect[i] < vect[1:i-1])
    }))
  }else{
    keep <- c(T,sapply(2:length(vect),function(i){
      !any(sapply(vect[1:i-1],all.equal,current = vect[i]) == "TRUE") & all(vect[i] > vect[1:i-1])
    }))
  }
  
  return(keep)
}