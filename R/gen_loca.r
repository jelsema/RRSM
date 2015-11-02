#' Generate locations over a spatial domain
#' 
#' @description
#' Generates a pattern of locations over a multiple of the unit square. Various point patterns are permitted.
#' 
#' @param method the pattern of points to generate, allowed values are grid, uniform, radiual, and spots.
#' @param ns the number of locations to generate.
#' @param unit the dimension of the domain, will be \code{unit} x \code{unit}
#' @param sphere logical, if \code{TRUE}, then the locations are projected onto a sphere.
#' @param x.num for \code{method="grid"}, the number of locations in the x-direction.
#' @param y.num for \code{method="grid"}, the number of locations in the y-direction.
#' @param alpha for \code{method="radial"}, the first parameter of the beta distribution.
#' @param beta for \code{method="radial"}, the second parameter of the beta distribution.
#' @param rad for \code{method="spots"}, the approximate (relative) radius of each spot (not on same scale as unit).

#' @param ... space for additional arguments
#' 
#' 
#' @return
#' A matrix of locations.
#' 
#' @examples
#' set1 <- gen_loca( method="grid", x.num=4, y.num=5 )
#' set2 <- gen_loca( method="radial", ns=5000, alpha=1, beta=2.5 )
#' 
#' @export
#' 

gen_loca <- function( method="grid" , ns=NULL , unit=100 , sphere=FALSE, x.num=NULL , y.num=x.num , 
                      alpha=1 , beta=2.5 , rad=0.2 , ... ){
  
  # Set ns if not specified
  if( is.null(ns) ){
    ns <- 10000
  }
  
  # Force ns to be multiple of 100 for spots pattern
  if( (ns %% 100 != 0) & ( method == "spots" ) ){
    ns <- ceiling(ns/100)*100    
  }
  
  ## ##################################################
  ## Generate locations on a regular grid
  ##
  
  if( method=="grid" ){
    
    # Get as a regular grid with as close as possible to the given sample size
    if( is.null(x.num)==T ){
      nsf <- floor( sqrt( ns ) ) - 5
      nsc <- ceiling( sqrt( ns ) ) + 5
            
      ylens <- xlens <- sort( rep( nsf:nsc , 2 ) )
      ylens <- ylens[ -length(ylens) ]
      xlens <- xlens[ -1 ]
      coords.numbers <- xlens*ylens
      ind1  <- which( abs(coords.numbers-ns) == min(abs(coords.numbers-ns)) )
      x.num <- max( xlens[ind1] )
      y.num <- max( ylens[ind1] )
    }
    
    x.sep <- unit/x.num
    y.sep <- unit/y.num
    
    xknot <- seq( 0 , unit-x.sep , x.sep ) + x.sep/2
    yknot <- seq( 0 , unit-y.sep , y.sep ) + y.sep/2
    locas <- expand.grid( xknot , yknot )    
    
  }
  
  ## ##################################################
  ## Generate locations randomly across domain (approximately uniform but non-gridded)
  ##
  
  if( method=="uniform" ){
    locas <- cbind( runif(ns,0,unit) , runif(ns,0,unit) )      
  }
    
  
  ## ##################################################
  ## Radial distributed locations using scaled BETA distribution
  ##
  
  if( method=="radial" ){
    locas <- cbind( rbeta(ns,alpha,beta) , rbeta(ns,alpha,beta) )*unit
  }
  
  
  ## ##################################################
  ## Generate clustered data in five "SPOTS"
  ##
  if( method=="spots" ){
    
    nss <- ns*c( 0.18 , 0.18 , 0.18 , 0.18 , 0.18 , 0.1 )
    
    centers <- cbind( c(0.25,0.25,0.50,0.75,0.75) , c(0.25,0.75,0.5,0.25,0.75) ) * unit
    
    locas  <- c(0,0)
    
    ## Points inside the circles
    for( i in 1:length(centers[,1]) ){
      prop <- cbind( rnorm( nss[1] ) , rnorm( nss[1] ) ) * (unit*rad/3)
      
      prop[,1] <- prop[,1] + centers[i,1]
      prop[,2] <- prop[,2] + centers[i,2]
      
      locas <- rbind( locas , prop )
    }
    
    count <- 0
    
    ## Points outside the circles
    while( count < nss[length(nss)] ){
      prop <- runif(2,0,unit)
      accept <- ( (centers[,1]-prop[1])^2 + (centers[,2]-prop[2])^2 )/unit^2 <= 0.9*(rad^2) 
      if( sum(accept)==0 ){
        locas <- rbind( locas , prop )
        count <- count+1
      }            
    }
    
    locas  <- locas[-1,]    
    
  }

  
  ##
  ## Return the locas
  ##
  
  colnames(locas) <- c("Xdim","Ydim")
  locas <- as.matrix( locas )
  
  return( locas )
  
}
