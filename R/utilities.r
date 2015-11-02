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





