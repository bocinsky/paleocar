#' Remove the duplicate entries in a matrix-like object
#'
#' This is a convenient function that removes the duplicated entries (rows)
#' in a matrix, but keeps a map of which of the original indices 
#' match the duplicates.
#'
#' @param mat A matrix or matrix-like object such as a data.frame or data.table.
#' @return A named list containing a matrix of the unique rows of the matrix and 
#' a vector matching the original indices to rows in the unique matrix.
duplicates <- function(mat){
  out <- rep(NA,ncol(mat))
  temp <- apply(mat, 2, function(x) paste(x, collapse="\r"))
  out[duplicated(temp)] <- match(temp[duplicated(temp)], temp)
  out[which(is.na(out))] <- which(is.na(out))
  unique.out <- unique(out)
  unique.mat <- mat[,unique.out]
  colnames(unique.mat) <- NULL
  out <- match(out, unique.out)
  return(list(models=unique.mat,matches=out)) 
}