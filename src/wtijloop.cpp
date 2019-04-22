// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector wtijloop( float cc, NumericVector ahat) {
int n1 = ahat.size();
int m1 = n1*(n1-1)/2;

Rcpp::NumericVector wtij(m1);

long double K;
long double cc2 = cc*cc;

int idx = -1;
for (int ii=0; ii < n1; ii++){
  for (int jj=ii+1; jj < n1; jj++){
    idx++;
    K = cc2 / ( abs(ahat[ii])*abs(ahat[jj]) ) ;
    
    if( K<1 ){
      wtij[idx] = K;
    } else{
      wtij[idx] = 1;
    }
  }
}

return(wtij);
}
