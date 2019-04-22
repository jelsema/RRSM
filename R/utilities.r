


#' Create a grid of points close to a target number
#' 
#' @description
#' Given a target number of points, creates an M x N grid of points, such that |M-N| <= 1
#' 
#' @param nk target number of points
#' @param kmin minimum number of points of the grid in one direction
#' @param kmax maximum number of points of the grid in one direction
#' 
#' @return
#' A matrix of locations
#' 
#' @export
#' 

approx_nk_grid <- function(nk , kmin=3 , kmax=30){
  ylens <- xlens <- sort( rep( kmin:kmax , 2 ) )
  ylens <- ylens[ -length(ylens) ]
  xlens <- xlens[ -1 ]
  knot.numbers <- xlens*ylens
  
  ind1 <- which( abs(knot.numbers-nk) == min(abs(knot.numbers-nk)) )
  xlen <- max( xlens[ind1] )
  ylen <- max( ylens[ind1] )
  c(xlen,ylen)
  
  gen_loca( method="grid", x.num=xlen, y.num=ylen )
  
}




#' Create an identity matrix
#' 
#' @description
#' Port of Matlab function \code{eye}.
#' 
#' @param a dimension of matrix to create.
#' 
#' @return
#' Itentity matrix of size \code{a}
#' 
#' @export
#' 

eye <- function(a){
  diag( a )
}


#' Pairwise sums and differences
#' 
#' @description
#' Computes the sums and differences between each pair of values of two vectors.
#' 
#' @param A a numeric vector.
#' @param B a numeric vector.
#' 
#' @return
#' A list with elements \code{S} and \code{D}
#' 
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @useDynLib RRSM
#' 
#' @export
#' 

pSD <- function(A,B){
  ss <- pairSums(A,B)
  dd <- pairDiffs(A,B)
  return( list(S=ss, D=dd) )
}



#' Find square root of matrix
#' 
#' @description
#' Computes the square root of a matrix by square-rooting the eigenvalues.
#' 
#' @param A dimension of matrix to create.
#' @param ... space for additional arguments (currently unused).
#' 
#' @return
#' Identity matrix of size \code{a}
#' 
#' @export
#' 

sigma12 <- function(A, ...){
  ## REPLACE THIS WITH A "FASTER" VERSION
  ## AND MAYBE USE A RIDGE TO PREVENT SINGULARITY PROBLEMS
  
  Aeig <- eigen(A)
  Avec <- Aeig$vectors
  Aval <- Aeig$values
  ## Corrects for numerical issues that sometimes appear
  if( any(Aval <= 0) ){
    Aval[Aval<=0] <- 0
  }
  B <- Avec %*% ( sqrt(diag(Aval)) ) %*% t(Avec)
  return(B)
}



#' Function used in lifting eigenvalues of 
#' 
#' @description
#' Compares the difference in the sum of eigenvalues
#' 
#' @param aa a multiplier
#' @param lam0 the smallest permissible eigenvalue, eigenvalues below this value will be lifted
#' @param EigVal the vector of observed eigenvalues
#' 
#' @return
#' The difference in the sum of eigenvalues
#' 
#' 
tr_var <- function(aa, lam0, EigVal ){
  sEigVal <- EigVal
  idx1    <- which( EigVal < lam0 )
  sEigVal[idx1] <- lam0 * exp( aa*(sEigVal[idx1]-lam0) )
  return( sum(EigVal) - sum(sEigVal)  )
}


#' Voronii
#' 
#' @description
#' For a vector of values, computes indicator vector for the minimum value. When multiple 
#' values are minimum, randomly selects one of the values to indicate.
#' 
#' @param x vector of values (usually distances between locations and knots).
#' 
#' @return
#' Indicator vector of same length as \code{x} (typically denotes closest knot to given location).
#' 
#' @export
#' 

voronii <- function(x){
  out <- ( which(x == min(x) ) )
  if( length(out)>1 ){
    select <- sample( which( x==min(x) ) , 1 , replace=FALSE )
    out <- 0*x
    out[select] <- 1
  }
  return(out)
}



#' Create a null matrix
#' 
#' @description
#' Port of Matlab function \code{zeros}.
#' 
#' @param a number of rows.
#' @param b number of columns.
#' 
#' @return
#' Null matrix with \code{a} rows and \code{b} columns.
#' 
#' @export
#' 
zeros <- function(a,b=1){
  matrix( 0 , nrow=a , ncol=b )
}


################################################################
##
## FUNCTIONS PLANNED BUT NON-FUNCTIONAL
##

#
#voronii_bins <- function(x, y, dist.args=NULL ){   
#  
#  nx <- nrow(x)
#  
#  for( ii in 1:nx ){
#    
#  }
#  return(out)
#}
#


##
## THESE CAN BE MADE MORE EFFICIENT
##


## More efficient version(s) ?
## http://stackoverflow.com/questions/24314878/compute-all-pairwise-differences-within-a-vector-in-r


if( FALSE ){
  library("Rcpp")
  
  
  #include <Rcpp.h>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
  
  
  
  
  cppFunction(  "RcppExport Rcpp::NumericVector pairDiffs( NumericVector A, NumericVector B) {
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
                
                return(allSum);
}")
  
  
  pairDiffs(1:3, 2:4)
  
  
  
  
  
  
}



























