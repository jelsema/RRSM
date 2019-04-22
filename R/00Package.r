
##
## This is the PACKAGE documentation
##

## Run this when adding/changing the C++ functions
#   # Delete  *.o and *.sh from /src
#   library("Rcpp")
#   compileAttributes()

#  LOOK INTO THIS:
#  RcppArmadillo.package.skeleton()



#' Reduced rank spatial models
#' 
#' @docType package
#' @name RRSM
#' 
#' @description
#' Various methods for reduced rank spatial analysis, including Empirical Bayesian and fully Bayesian methods.
#' Important note: This package is still a developmental version. It is only available because it is
#' the easiest way in which to make available the codes for one of the methods implemented (\code{TKS}) as
#' a supplement for a submitted manuscript. 
#' 
#' 
#' Appropriate credit should be given when publishing results obtained using \pkg{RRSM}, or when 
#' developing other programs/packages based off of this one. Use \code{citation(package="RRSM")} for 
#' Bibtex information.
#' 
#' @author blinded <\email{blinded@blinded}>
#' 
#'
#' @import methods
#' @import stats
#'
#' @useDynLib RRSM
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' 
"_PACKAGE"












