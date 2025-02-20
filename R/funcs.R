#' @useDynLib TESSA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' SVD a matrix using RcppArmadillo
#'
#' @param X A matrix
#' @return The singular values of SVD
#' @export
baseSVD <- function(X) {
  .Call('_TESSA_baseSVD', X)
}

#' Inverse of a matrix using RcppArmadillo
#'
#' @param X A matrix
#' @return The inverse matrix
#' @export
invert <- function(X) {
  .Call('_TESSA_invert', X)
}
