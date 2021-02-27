#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat get_ev_from_evd(const arma::mat& mat) {
  arma::mat cor_mat(mat.n_cols, mat.n_cols);
  cor_mat = arma::cor(mat);
  return arma::eig_sym(cor_mat);
}

// [[Rcpp::export]]
arma::mat get_ev_from_svd(const arma::mat& mat) {
  // standard SVD
  const int df = mat.n_rows - 1;
  arma::vec s;
  s = arma::svd(mat);
  return arma::square(s) / df;
}


// arma::mat cor_cpp(const arma::mat& mat) {
//   arma::mat cor_mat(mat.n_cols, mat.n_cols);

//   if (mat.n_rows == 1) {
//     cor_mat.fill(NA_REAL); // to mimic R base::cor function
//   } else {
//     cor_mat = arma::cor(mat); // default is N-1 normalization
//   }

//   return cor_mat;
// }


// arma::mat get_ev_from_svd(const arma::mat& mat) {
//   // economic SVD
//   const int df = mat.n_rows - 1;
//   arma::mat U;
//   arma::vec s;
//   arma::mat V;
//   svd_econ(U, s, V, mat, "left");
//   return arma::square(s) / df;
// }

