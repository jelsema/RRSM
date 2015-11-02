#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
RcppExport NumericVector pairSums( NumericVector A, NumericVector B) {
int n1 = A.size();
int n2 = B.size();

NumericVector allSum(n1*n2);
int idx = -1;
for (int ii=0; ii < n1; ii++){
for (int jj=0; jj < n2; jj++){
idx++;
allSum[idx] = A[ii] + B[jj];
}
}

return(allSum);
}
