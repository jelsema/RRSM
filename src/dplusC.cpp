// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
RcppExport long double dplusC( NumericVector xx) {

int n1 = xx.size();

// Get the HL estimator
Rcpp::NumericVector pairSum( (n1*(n1+1))/2 );
int idx = -1;
for (int ii=0; ii < n1; ii++){
 for (int jj=ii; jj < n1; jj++){
  idx++;
  pairSum[idx] = (xx[ii] + xx[jj])/2;
 }
}

std::sort( pairSum.begin(), pairSum.end() );
int size = pairSum.size();
long double mid = size/2;

long double HL ;
if( size % 2 == 0 ){
 HL = (pairSum[mid] + pairSum[mid-1]) / 2 ;
}
else{
 HL = pairSum[mid] ;
}

// Subtract HL from values
long double yy1;
Rcpp::NumericVector yy( n1 );
for (int ii=0; ii < n1; ii++){
  // long double  y1     = xx[ii] - HL ;
  yy1     = xx[ii] - HL ;
  yy[ii] = ((yy1>0) - (yy1<0))*yy1 ;
}

std::sort( yy.begin(), yy.end() );

int runStart = 0 ;
int runStop  = 0 ;
long double rankSum  = 0 ;
int gomore   = 0 ;
long double avRank = 0.0;

int ii  = 0 ;
int ii2 = 0;
//int ii3 = 0;

Rcpp::NumericVector ranks(n1);


// Get the midranks
while( ii <= n1 - 1 ){
  
  if( ii==(n1 - 1) ){
    ranks[ii] = ii + 1;
    ii = ii + 1;
  } else{
    
    runStart = ii;
    runStop  = ii;
    rankSum  = 0;
    ii2      = ii;
    gomore   = 0;
    
    while( gomore==0 & ii2<=(n1 - 1) ){
      
      long double diff = (yy[ii] - yy[ii2]);
      long double adiff= ((diff>0)-(diff<0))*diff;
      
      if( adiff == 0 | adiff==-0 ){
        rankSum = rankSum + (ii2+1);
        ii2     = ii2 + 1 ;
      } else{
        gomore = 1;
      }
    }
    runStop = ii2 - 1;
    
    for( int ii3=runStart; ii3<=runStop; ii3++){
	  //long double avRank = rankSum /(long double) (runStop-runStart+1) ;  
      avRank = rankSum /(long double) (runStop-runStart+1) ;  
      ranks[ii3] = avRank;
    }
    
    ii = runStop + 1 ;
  }
  
}

// Finally produce dplus

long double ss1 = sqrt(3);

rankSum = 0;
for( int ii=0; ii<n1; ii++ ){
  rankSum = rankSum + (ranks[ii]*yy[ii]) ;
}

long double dplus = (ss1 /(long double) (n1+1))*rankSum;

return( dplus );

}
