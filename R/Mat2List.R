#' Split a stacked matrix into a list of per-group matrices
#'
#' Splits a vertically stacked matrix into a list of matrices,
#' where the row counts are determined by a reference list of matrices.
#' Useful for converting a single W_RNA matrix (from e.g. nnTensor::NMF)
#' into a per-chromosome list for Machima2's list mode.
#'
#' @param X_list A list of matrices whose nrow() values define the split points
#' @param W A matrix to split (total rows must equal sum of nrow(X_list[[k]]))
#' @return A list of matrices with the same length as X_list
#' @examples
#' X_list <- list(matrix(0, 10, 5), matrix(0, 20, 5))
#' W <- matrix(runif(30*3), 30, 3)
#' W_list <- Mat2List(X_list, W)
#' # nrow(W_list[[1]]) == 10, nrow(W_list[[2]]) == 20
#' @export
Mat2List <- function(X_list, W){
    .Mat2List(X_list, W)
}
