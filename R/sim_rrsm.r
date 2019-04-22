#' Simulate large spatial dataset
#' 
#' @description
#' Simulates a large spatial dataset according to a Spatial Random Effects model: \eqn{Y = X*beta + S*eta + epsilon}{Y = X\beta + S\eta + \epsilon}
#' 
#' @param X the design matrix for the large-scale variation (fixed effects).
#' @param S the matrix of basis functions (design matrix for the random effect).
#' @param beta the large-scale variation vector (fixed-effects slopes).
#' @param V the covariance matrix of the reduced rank process.
#' @param ssq the residual variance.
#' @param theta the spatial dependence parameter, for if \code{V} needs to be calculated.
#' @param coords the observed locations, for if \code{S} needs to be calculated.
#' @param knots the knot locations, for if \code{S} needs to be calculated.
#' @param dist the distribution to use for simulating data. Currently only "mvrnorm" is supported. 
#' @param seed optional means to control the RNG seed.
#' @param ... space for additional arguments.
#' 
#' @return
#' A vector containing the simulated values.
#' 
#' @note
#' Several arguments, including \code{S}, \code{V}, and \code{ssq} may be left empty, 
#' if locations and knots are provided to calculate them.
#' 
#' @importFrom MASS mvrnorm
#' 
#' @export
#' 


sim_rrsm <- function( X, S, beta=NULL, V=NULL,
                               ssq=0.001, theta=NULL, coords=NULL, knots=NULL,
                               dist="mvrnorm", seed=NULL, ... ){
  
  if( is.numeric(seed) ){
    if( (seed == ceiling(seed)) & (seed == floor(seed)) ){
      set.seed(seed)
    }
  }
  
  ## Is X provided?
  if( missing(X) ){
    stop("The design matrix 'X' is missing")
  }
  
  ## Is S provided?
  if( missing(S) ){
    stop("Matrix of basis functions 'S' is missing. See function `basis_mat' ")
  }
  
  ## Is V provided?
  ## - Work in alternate covariance structures here when time permits
  ## - e.g. a matern class
  if( is.null(V) ){
    if( is.null(theta) || is.null(knots) ){ stop("Missing one of V or {theta, knots}") }
    warning("Covariance matrix 'V' was not provided. Constructing 'V' based 'knots' and 'theta', but unsupervised construction of 'V' is not recommended.")
    warning("Using exponential covariance function to construct 'V' ")
    
    kk_dist <- pairwise_dist( knots, knots, ... )
    V       <- exp( -1*kk_dist/theta )
  }
  
  ## Simulate the data
  nd <- nrow(S)
  nk <- ncol(S)
  
  eta <- mvrnorm( 1, rep(0, nk), V )
  eps <- rnorm( nd, 0, sqrt(ssq) )
  
  Y <- X%*%beta + S%*%eta + eps
  
  return(Y)
  
}
