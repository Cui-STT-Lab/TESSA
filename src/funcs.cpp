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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double TT_cpp(const arma::mat& M1, const arma::mat& M2) {
  if (M1.n_cols != M2.n_rows || M1.n_rows != M2.n_cols) {
    Rcpp::stop("Non-conformable matrices in TT_cpp().");
  }

  double result = 0.0;
  for (arma::uword i = 0; i < M1.n_rows; ++i) {
    for (arma::uword j = 0; j < M1.n_cols; ++j) {
      result += M1(i, j) * M2(j, i);
    }
  }
  return result;
}
