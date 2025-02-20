#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec baseSVD(const arma::mat& X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "standard");
  return S;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat invert(const arma::mat& X) {
  return arma::inv(X);
}
