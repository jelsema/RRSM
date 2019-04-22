// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector pairDiffs( NumericVector A, NumericVector B) {
  int n1 = A.size();
  int n2 = B.size();
  int n3 = n1 * n2;
  
  Rcpp::NumericVector allSum(n3);
  int idx = -1;
  for (int ii=0; ii < n1; ii++){
    for (int jj=0; jj < n2; jj++){
      idx++;
      allSum[idx] = A[ii] - B[jj];
    }
  }
  
  return( allSum );
}
