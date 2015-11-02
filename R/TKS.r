#' Estimate knots for reduced rank spatial process
#' 
#' @description
#' For a given bandwith, estimate knots for a reduced rank spatial models using the centroids of interval sets, iterating through all thresholds.
#' 
#' @param coords the observed locations.
#' @param R the data to be used in estimating knots.
#' @param delta1 the first bandwidth, used to determine the knots for a given threshold interval set, passed to \code{\link{threshold_knots}}.
#' @param delta2 the second bandidth, used to detect if newly estimated knots are too close to existing knots.
#' @param slices the threshold cut-off values to define interval sets.
#' @param augment whether to augment the estimated set of knots with a regular grid of knots (fills in sparse areas). See Notes.
#' @param knots_grid (optional) the grid of knots to use if \code{augment=TRUE}.
#' @param ... space for additional arguments.
#' 
#' @return
#' A matrix of locations.
#' 
#' @note
#' The argument \code{augment} was considered in the course of development of the Threshold Knot 
#' Selection (TKS) algorithm, but was not pursued with much depth. The properties of TKS when 
#' augmenting with a grid have not been studied.
#' 
#' @seealso
#' \code{\link{TKS}}
#' 
#' @export
#' 

empirical_knots <- function( coords , R , delta1 , delta2=delta1 , slices , 
                             augment=FALSE , knots_grid=NULL, ... ){
  
  absR <- abs(R)
  
  # Method 1: Evenly spaced VALUES based on percentiles
  # Us <- quantile( absR , c( 0.98 , 0.3 ) )
  # Us <- seq( Us[1] , Us[2] , length.out=30 )
  
  # Method 2: Evenly spaced PERCENTILES
  # Us <- quantile( absR , slices  )
  
  Us <- quantile( absR , slices  )    
  
  ##
  ## Iteration 1
  ##
  
  knots1 <- threshold_knots( coords ,  R , Us[1] , delta1)
  knots2 <- threshold_knots( coords , -R , Us[1] , delta1)
  knots  <- rbind( knots1 , knots2 )  
  
  ##
  ## Loop through the rest of the percentiles
  ##
  if( length(Us)>=2 ){
    for( ii in 2:length(Us) ){
      
      # Method 1: Not strictly Excursion sets, but "Interval" set
      ExSet <- ( absR <= Us[ii-1]) & (absR >= Us[ii])
      
      # Method 2: Excursion sets (strictly growing clusters, which may cause
      #           problems for large, connected areas)
      # ExSet <- (absR >= Us[ii])
      
      locaii <- coords[ ExSet==TRUE , ]
      Rii    <- R[ ExSet==TRUE ]
      
      # Get the knots for the next percentile
      knew1 <- threshold_knots( locaii ,  Rii , Us[ii] , delta1, ...)
      knew2 <- threshold_knots( locaii , -Rii , Us[ii] , delta1, ...)
      k_new <- rbind( knew1 , knew2 )  
      
      # Check to see if any are too close to existing percentiles
      nk2 <- dim(k_new)[1]
      if( is.null(nk2) ){ nk2 <- 0 }
      if( nk2 > 0 ){
        drop.knot <- vector( length=nk2 )
        
        for( jj in 1:nk2 ){
          dists     <- pairwise_dist( knots , k_new[jj,,drop=FALSE] , ...)
          too.close <- (sum( dists< delta2 ) > 0)*1    
          drop.knot[jj] <- too.close
        }
        
        k_new <- k_new[ drop.knot==0 , , drop = FALSE ]
        
        # Add the 'unique' new knots to the current knots
        if( dim(k_new)[1] > 0 ){
          knots <- rbind( knots , k_new )
        }
      }
    } # END: FOR-ii
  }
  
  knots  <- merge_close_knots( knots , delta2 , ...)
  knots1 <- knots
  
  
  if( augment == TRUE ){
    # Merge in any good knots from a grid of approximately the same size
    if( is.null(knots_grid) ){
      k_new <-  as.matrix( gen_loca( method="grid" , x.num=ceiling(sqrt(dim(knots)[1])) ) )
    } else{
      k_new <- knots_grid
    }
    
    nk2 <- dim(k_new)[1]
    drop.knot <- vector( length=nk2 )
    for( jj in 1:nk2 ){
      dists     <- pairwise_dist( knots , k_new[jj,,drop=FALSE] , ...)
      too.close <- (sum( dists< delta2 ) > 0)*1    
      drop.knot[jj] <- too.close
    }
    k_new <- k_new[ drop.knot==0 , , drop = FALSE]
    if( dim(k_new)[1] > 0 ){
      knots <- rbind( knots , k_new )
    }
  }
  
  ## Return knots
  return(knots)
}




#' Merge close knots
#' 
#' @description
#' Merges (via centroid) locations that are within some bandwidth of each other.
#' 
#' @param knots the observed locations (usually, knot locations)
#' @param delta the bandwidth to determine which locations are 'close'.
#' @param ... space for additional arguments.
#' 
#' @return
#' A matrix containing a set of locations.
#' 
#' @seealso
#' \code{\link{TKS}}
#' 
#' @export
#' 
merge_close_knots <- function( knots , delta , ... ){
  # Merge any knots that are within delta of each other
  nk <- dim(knots)[1]
  k_dist <- pairwise_dist( knots , knots, ...)
  k_dist <- k_dist * upper.tri(k_dist)
  k_dist[ which( k_dist==0 ) ] <- Inf
  
  combine.any <- min( k_dist ) < delta
  
  while( combine.any ){
    
    ind.min <- arrayInd( which.min(k_dist), dim(k_dist) )
    
    k1 <- ind.min[1]
    k2 <- ind.min[2]
    
    knots[k1,] <- (knots[k1,] + knots[k2,])/2
    knots <- knots[-k2,,drop=FALSE]
    
    nk <- dim(knots)[1]
    k_dist <- pairwise_dist( knots , knots, ... )
    k_dist <- k_dist * upper.tri(k_dist)
    k_dist[ which( k_dist==0 ) ] <- Inf
    
    combine.any <- min( k_dist ) < delta
  }
  
  ## Return knots
  return(knots)
  
}



#' Exceedence set knots
#' 
#' @description
#' For a given bandwidth, determines the locations 
#' 
#' @param coords the observed locations.
#' @param R the data to be used to define knots.
#' @param u the current threshold.
#' @param delta the bandwidth to determine which locations are 'close'.
#' @param ... space for additional arguments.
#' 
#' @return
#' A matrix containing a set of locations.
#' 
#' @seealso
#' \code{\link{TKS}}
#' \code{\link{empirical_knots}}
#' 
#' @export
#' 

threshold_knots <- function( coords , R , u , delta, ...){
  
  #
  # Defines the connected components for the excursion set of a spatial process
  #
  # Excursion / Threshold exceedence set: Au := { x in R | x >= u }
  # Connected components: clusters of elements in Au within some bandwidth (delta)
  #
  # Code based off of Matlab function by Madrid (2010)
  #
  # INPUT:
  # coords : Matrix of locations
  # R    : Vector of response data
  # u    : Threshold
  # delta: Bandwidth
  # ...  : Extra arguments passed to other functions
  #
  # OUTPUT:
  # num:  number of connected components
  # ind:  index of elements above the threshold (u)
  # comp: label (number) of component for each element above the threshold
  #
  
  ind  <- which( R >= u)  # Shouldn't this be an interval??
  nu   <- length(ind)
  
  # Only do anything if any values above the threshold
  if( nu>0 ){
    
    # Set up matrix to organize results
    comp     <- matrix( 0 , nrow=nu , ncol=3)
    comp.idx <- 1
    
    comp[1,1] <- ind[1]
    comp.area <- vector()
    
    lake <- ind
    
    
    # Loop through the elements
    while(   sum(comp[,1]!=0) != sum(comp[,3]!=0)   ){
      
      # Locate those locations that have been included, but not explored
      idx <- which( comp[,1]!=0 & comp[,3]==0 )
      t1  <- min(idx)
      t2  <- max(idx)
      
      # Give points component label and mark them as explored
      comp[ idx , 2 ] <- rep( comp.idx , length(idx) )
      comp[ idx , 3 ] <- rep( 1 , length(idx) )
      
      # Remove elements from the lake
      hose <- match( comp[comp[,1]!=0,1] , ind)
      lake <- ind[ -hose ]
      
      for( jj in idx ){
        if( length(lake)>0 ){ 
          
          dists <- pairwise_dist( coords[lake,,drop=FALSE] , 
                                  coords[ comp[jj,1],,drop=FALSE] , ...)
          hose  <- which( dists<=delta )
          
          if( length(hose)>0 ){
            
            t1 <- min( which(comp[,2]==0) )
            t2 <- t1 + length(hose) - 1
            
            comp[ t1:t2 , 1 ] <- lake[hose]
            comp[ t1:t2 , 2 ] <- rep( comp.idx ,t2-t1+1 )
            
            lake <- lake[ -hose ]
          } # End IF-length(hose)>0 #2
          
        }
      } # End FOR-jj
      
      ## Detect if the loop will exit before empying the lake
      if(  (sum(comp[,1]!=0) == sum(comp[,3]!=0)) & length(lake)>0  ){
        
        
        # comp.area
        
        comp.idx   <- comp.idx + 1
        t1         <- min( which(comp[,1]==0) )          
        comp[t1,1] <- lake[1]
        
      }
      
    } # End WHILE-mismatch
    
    ## Look for "large" components and separate them
    num <- max( comp[,2] )
    
    # Obtain the knots
    ExSet <- list( num=num , ind=comp[,1] , comp=comp[,2] )  
    LPK   <- as.numeric( table(comp[,2]))
    
    abth <- coords[ ExSet$ind , ]
    cent <- zeros( ExSet$num , 2 );
    for(i in 1:ExSet$num ){
      if( LPK[i]>1 ){
        cent[i,] = colMeans( abth[ ExSet$comp==i ,1:2] )
      } else{
        if( LPK[i]==1 ){
          cent[i,] <- coords[ comp[ which(comp[,2]==i) , 1 ] , ]
        }
      }
      
    }
    
    # Return the knots
    return( cent )
    
    # No values above the threshold?  
  } else{
    # Need to return something that wont cause an error?
    return(NULL)
  }
  
}





#' Threshold knot selection
#' 
#' @description
#' Uses data to estimate the number and locations of knots for a reduced rank spatial process.
#' 
#' @param Y the observed data. Should already be on log-scale.
#' @param X the design matrix for the large-scale variation (fixed effects).
#' @param coords the observed locations.
#' @param R the residuals to use to select knots (computed if missing).
#' @param test_set the indices of locations to use as the hold-out test set.
#' @param M the number of locations to use for a hold-out test set (ignored if \code{test_set} is non-null).
#' @param penalty the exponent in the penalty for the objective function of bandwidth selection. The MSE is weighted against the number of knots to this power.
#' @param BWs the bandwidths to test.
#' @param num_bw if \code{BWs} is null, this will set the number of bandwidths to test.
#' @param bw_int if \code{BWs} is null, the maximum distance in the coordinates is multiplied by the min and max of \code{bw_int} to set the minimum and maximum bandwidths to test.
#' @param slices the threshold cut-off values to define interval sets, defined as percentiles.
#' @param num_slice if \code{slices} is null, this will set the number of percentile slices.
#' @param slice_int if \code{slices} is null, the slices will be evenly spaces between the min and max of \code{slice_int}.
#' @param vpred method to estimate parameters for the validation stage of TKS.
#' @param nres number of resolutions of knots to seelct (currently only 1 or 2 resolutions can be handled).
#' @param ... space for additional arguments.
#' 
#' @return
#' A matrix containing a set of locations.
#' 
#' @seealso
#' \code{\link{empirical_knots}}
#' \code{\link{threshold_knots}}
#' 
#' @importFrom MASS ginv
#' @importFrom fields cover.design
#' 
#' @export
#' 

TKS <- function( Y, X, coords, R, test_set=NULL, M=round(0.1*nrow(X)), penalty=0.5, 
                 BWs=NULL, num_bw=NULL, bw_int=c(0.03,0.12),
                 slices=NULL, num_slice=NULL, slice_int=c(0.98,0.2), vpred="EM", nres=1, ... ){
  
  nd <- nrow(X)
  
  if( missing(R) ){
    R  <- Y - X %*% (ginv(t(X)%*%X) %*% t(X)%*%Y)
  }
  
  if( is.null(test_set) ){ test_set <- sample( 1:nd , M )  }
  
  if( is.null(slices) ){
    if( is.null(num_slice) ){ num_slice <- 15  }
    slices <- seq( max(slice_int) , min(slice_int) , length.out=num_slice )
  }
  
  
  if( is.null(BWs) ){
    if( is.null(num_bw) ){ num_bw <- 6 }
    
    xlen <- max( coords[,1] ) - min( coords[,1] )
    ylen <- max( coords[,2] ) - min( coords[,2] )
    
    dlen <- round(sqrt( xlen^2 + ylen^2 ))
    
    sbw <- dlen*min(bw_int)
    lbw <- dlen*max(bw_int)
    
    BWs <- seq( sbw , lbw , length.out= num_bw )
    
  }
    
  
  Y_test  <- Y[  test_set  ]
  X_test  <- X[  test_set ,]
  coords_test <- coords[ test_set , ]
  
  Y_train <- Y[ -test_set  ]
  X_train <- X[ -test_set ,]
  R_train <- R[ -test_set  ]
  coords_train <- coords[ -test_set , ]
    
  BandSel <- matrix( 0 , nrow=length(BWs) , ncol=3 )
  min_SS  <- Inf
  all_knot_sets <- list()
  
  ## Loop through the bandwidths
  for( ii in 1:length(BWs) ){
    
    h1 <- BWs[ii]
    #h2 <- 1.5*BWs[ii]
    h2 <- h1
    
    knots <- empirical_knots( coords=coords_train , R=R_train , delta1=h1 , delta2=h2 ,
                              slices=slices , ... )
    
    all_knot_sets[[ii]] <- knots
    
    # dist.args=dist.args)
    
    nk_bw  <- nrow(knots)
    
    
    ## Predictions / MSE using FRK type of model
    S      <- basis_mat( coords , knots, ... )
    S_test  <- S[  test_set , ]
    S_train <- S[ -test_set , ]
    
    if( vpred=="MOM"){
      # Method 1: MOM Estimation + Kriging
      ests_mom <- rr_est( method="MOM", empirical="cj", fitting="frobenius", 
                          Y=Y_train, X=X_train, S=S_train, coords=coords_train, ... )
      preds <- rr_universal_krige( Y=Y_train, X=X_train, S=S_train, coords=coords_train, 
                                   Xpred=X_test, Spred=S_test, pgrid=coords_test ,
                                   V = ests_mom$V, ssq=ests_mom$ssq, ... )$Preds
      Yhat  <- preds[,3]
    } else if( vpred == "EM"){
      # Method 2: EM Estimation + Kriging
      ests_mom <- rr_est( method="EM", Y=Y_train, X=X_train, S=S_train, coords=coords_train, ... )
      preds <- rr_universal_krige( Y=Y_train, X=X_train, S=S_train, coords=coords_train, 
                                   Xpred=X_test, Spred=S_test, pgrid=coords_test ,
                                   V = ests_mom$V, ssq=ests_mom$ssq, ... )$Preds
      Yhat  <- preds[,3]
    } else{
      # Method 3: Approximation
      Psi    <- cbind( X , S )    
      Yhat   <- c( lm( Y ~  Psi - 1 )$fitted.values )[test_set]
    }
    
    SSE   <- mean( (Y_test - Yhat)^2 )
    
    
    ## Predictions / MSE using PredProc type of model
    # --> Need to get the parameter estimates for a given covariance model
    
    ## Evaluate the objective function
    if( SSE*sqrt(nk_bw) < min_SS ){
      min_SS   <- SSE * (nk_bw^penalty)
      knots_r1 <- knots
      bw_selk  <- h1
    }
    
    BandSel[ii,1] <- h1
    BandSel[ii,2] <- nk_bw
    BandSel[ii,3] <- SSE
    
  }
  
  
  if( nres==2 ){
    # Select a second resolution of knots
    nk_sel          <- nrow(knots_r1)
    have_second_res <- FALSE
    slices_r2       <- slices
    
    iter <- 0
    
    if( nk_sel < 120 ){
      # Look for larger set of knots within those already found
      if( any( BandSel[,2]>= (nk_sel*3) ) ){
        idx             <- which.max( BandSel[,2] >= (nk_sel*3)  )
        have_second_res <- TRUE
        knots_r2        <- all_knot_sets[[idx]]
      }
      
      # Otherwise, need to do some more searching
      #slices_r2 <- seq( max(slice_int) , min(slice_int) , length.out= ceiling(length(slices)*1.25) )
      while( have_second_res==FALSE ){
        iter <- iter+1
        h1_r2 <- bw_selk*0.9
        knots_r2 <- empirical_knots( coords=coords_train , R=R_train , delta1=h1_r2 , delta2=h1_r2 ,
                                     slices=slices_r2 , ... )
        if( nrow(knots_r2) >= nk_sel*3 | iter>30 ){
          have_second_res <- TRUE
          if( nrow(knots_r2) < nk_sel*3 ){
            # knots_r2 <- gen_loca( method="grid" , ns=nk_sel * 3 )
            knots_r2 <- cover.design( coords_train[ duplicated(coords_train)==FALSE, ] , nd=(nk_sel*3) )
          }
        }
      }
    }  else{
      # Look for smaller set of knots within those already found
      if( any( BandSel[,2]<= ceiling(nk_sel/3) ) ){
        idx             <- which.min( BandSel[,2] <= ceiling(nk_sel/3) )
        have_second_res <- TRUE
        knots_r2        <- all_knot_sets[[idx]]
      }
      
      # Otherwise, need to do some more searching
      #slices_r2 <- seq( max(slice_int) , min(slice_int) , length.out= ceiling(length(slices)*0.75) )
      while( have_second_res==FALSE ){
        iter <- iter+1
        h1_r2 <- bw_selk*1.1
        knots_r2 <- empirical_knots( coords=coords_train , R=R_train , delta1=h1_r2 , delta2=h1_r2*1.5 ,
                                     slices=slices_r2 , ... )
        if( nrow(knots_r2) <= ceiling(nk_sel/3) | iter>30 ){
          have_second_res <- TRUE
          if( nrow(knots_r2) > ceiling(nk_sel/3) ){
            # knots_r2 <- gen_loca( method="grid" , ns=ceiling(nk_sel/3) )
            knots_r2 <- cover.design( coords_train[ duplicated(coords_train)==FALSE, ] , nd=ceiling(nk_sel/3) )
          }
        }
      }
    }
  } else{
    knots_r2 <- NULL
  }
  
  # plot( BandSel[,1] , BandSel[,2]*BandSel[,3] , type='b' )
  # cbind( BandSel[,1] , BandSel[,2]*BandSel[,3] )
  # BandSel
  # easy.map( cbind(coords,Y) , knots=knots_r1 )
  
  ## Return Objects
  colnames(BandSel) <- c("h1" , "r" , "SSE")
  Return.Obj <- list( knots=knots_r1, knots2=knots_r2 , BandSel=BandSel , h1=bw_selk )
  
  return( Return.Obj )
  
}




