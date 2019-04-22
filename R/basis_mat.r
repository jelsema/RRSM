#' Calculate pairwise distances
#' 
#' @description
#' Calculates pairwise distances between two sets of locations.
#' 
#' @param coords1 the coordinates of the first set of locations (dimension \eqn{n1 x p}).
#' @param coords2 the coordinates of the second set of locations (dimension \eqn{n2 x p}). 
#' @param dmethod the type of distance to compute, see Details.
#' @param ... space for additional arguments
#' 
#' @details
#' This function is a wrapper for existing distance functions. Argument \code{dmethod} selects from
#' two types of distances to be computed:
#' \itemize{
#' \item \code{Euclidean distance:}{ Computes distance using \code{SpatialTools::dist2}. 
#'        For this method, \code{identical(dmethod, SpatialTools::dist2)} should return \code{TRUE} (default behavior).}
#' \item \code{Great circle distance:}{ Distance on a sphere/ellipsoid using \pkg{geosphere}. Argument 
#'       \code{dmethod} should match one of the distance functions allowed for argument \code{fun} in \code{geosphere::distm}. }
#' }
#' 
#' @return
#' A vector of matrix containing distances between the points.
#' 
#' @examples
#' coords1 <- cbind( c(1,1,2) , c(1,2,1) )
#' coords2 <- cbind( c(3,3) , c(2,3) )
#' pairwise_dist( coords1, coords2 )
#' 
#' @import geosphere 
#' @importFrom SpatialTools dist2
#' 
#' @export
#' 
pairwise_dist <- function( coords1 , coords2=coords1 , dmethod=SpatialTools::dist2, ... ){
  
  imc1 <- is.matrix(coords1)
  imc2 <- is.matrix(coords2)
  
  
  if( identical(dmethod, SpatialTools::dist2) ){
    ## Euclidean distance
    all_dist <- SpatialTools::dist2( as.matrix(coords1), as.matrix(coords2) )  
  } else{
    ## Distance on the sphere
    # dmethod2 <- get( paste0("geosphere::",dmethod ) )
    all_dist <- geosphere::distm( coords1, coords2, ... ) / 1000
  }
  
  # Return distance
  return( all_dist )
  
}


if( FALSE ){
  pairwise_dist <- function( coords1, coords2=coords1, ... ){
    
    n1 <- nrow(coords1)
    n2 <- nrow(coords2)
    dist_mat <- zeros(n1, n2)
    
    if( n1 >= n2 ){
      for( ii in 1:n2 ){
        dist_mat[,ii] <- vector_dist( coords1, coords2[ii,], ... )
      }
    } else{
      for( ii in 1:n1 ){
        dist_mat[ii,] <- vector_dist( coords2, coords1[ii,], ... )
      }
    }
    
    return( dist_mat )
    
  }
}

#' Construct matrix of basis functions
#' 
#' @description
#' Constructs the matrix of basis functions from a set of locations and a set of knots.
#' 
#' @param coords the coordinates of the data locations (dimension \eqn{n1 x p}).
#' @param knots the coordinates of the knots (dimension \eqn{n2 x p}). 
#' @param mult multiplier for range of each basis function.
#' @param FUN the basis functions, currently bisquare (bisq) and modified bisquare (mod_bisq) are accepted.
#' @param ... space for additional arguments
#' 
#' @details
#' Distances are computed using \code{\link{pairwise_dist}}, see documention of that function
#' for different methods of calculating distance.
#' 
#' @note
#' In future releases, custom basis functions will be permitted.
#' 
#' @return
#' A matrix of size \eqn{n1 \times n2} containing pairwise distances.
#' 
#' @export
#' 
basis_mat <- function( coords, knots, mult=1.3, FUN="bisq", ...){
  
  nd <- nrow(coords)
  nk <- nrow(knots)
  
  S <- knot_ind <- zeros( nd , nk )
  
  kd_dist <- pairwise_dist( coords , knots, ... )
  k_dist  <- pairwise_dist( knots  , knots, ...)    
  
  ## Enable other choices of basis function
  ## Create the S-matrix (bisquare)
  if( FUN=="bisq" ){
    for( ii in 1:nk ){
      a1 <- mult*min( k_dist[k_dist[,ii]>0, ii]  )
      knot_ind[,ii] <- kd_dist[,ii] <= a1
      S[,ii] <- ( (1 - (kd_dist[,ii]/a1)^2 )^2 )*knot_ind[,ii]
    }
  }
  if( FUN=="mod_bisq" ){
    
    for( ii in 1:nk ){
      dd     <- sqrt( ((coords[,1]-knots[ii,1])/mult[1])^2 + ((coords[,2]-knots[ii,2])/mult[2])^2 )
      S[,ii] <- as.numeric(dd <= 2)*(1 - .25*dd)
    }
  }
  
  ## Return the S matrix
  return(S)
  
}


