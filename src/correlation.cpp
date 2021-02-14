#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cor_cpp(const arma::mat& mat) {
  arma::mat cor_mat(mat.n_cols, mat.n_cols);

  if (mat.n_rows == 1) {
    cor_mat.fill(NA_REAL); // to mimic R base::cor function
  } else {
    cor_mat = arma::cor(mat); // default is N-1 normalization
  }

  return cor_mat;
}

