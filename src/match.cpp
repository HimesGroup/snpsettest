#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector match_cpp(CharacterVector x, CharacterVector table) {
  return na_omit(match(x, table)); // do not return NA
}
