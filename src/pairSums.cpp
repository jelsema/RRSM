// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
 Rcpp::NumericVector pairSums( Rcpp::NumericVector A, Rcpp::NumericVector B) {
int n1 = A.size();
int n2 = B.size();

Rcpp::NumericVector allSum(n1*n2);
int idx = -1;
for (int ii=0; ii < n1; ii++){
for (int jj=0; jj < n2; jj++){
idx++;
allSum[idx] = A[ii] + B[jj];
}
}

return(allSum);
}
