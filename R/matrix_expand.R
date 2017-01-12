#' Expand a matrix, filling-down rows as needed
#'
#' Given a matrix \code{m} and a vector \code{v} of new row IDs, returns a new matrix
#' with \code{length(v)} rows, and fills the empty rows with data from the nearest preceding row
#' containing data.
#'
#' @param the.matrix A numeric matrix.
#' @param new.rows A character vector of row names, 
#' of which a subset whould be the row names of \code{the.matrix}.
#' @importFrom zoo na.locf
#' @return A numeric matrix.
matrix_expand <- function(the.matrix, new.rows){
  the.matrix <- data.matrix(the.matrix)
  missing.rows <- new.rows[!(new.rows %in% intersect(as.numeric(rownames(the.matrix)),new.rows))]
  new.rows.matrix <- matrix(nrow=length(missing.rows),ncol=ncol(the.matrix),data=NA)
  colnames(new.rows.matrix) <- colnames(the.matrix)
  rownames(new.rows.matrix) <- missing.rows
  matrix_expand <- rbind(the.matrix,new.rows.matrix)
  matrix_expand <- matrix_expand[order(as.numeric(rownames(matrix_expand))),]
  matrix_expand <- zoo::na.locf(matrix_expand)
  matrix_expand[is.na(matrix_expand)] <- 0
  return(matrix_expand)
}