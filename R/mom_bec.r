#' Binned Empirical Covariance Matrix
#' 
#' @description
#' Computes the binned empirical covariance matrix for Method of Moments (MOM) estimation 
#' in a reduced rank spatial model.
#' 
#' @param empirical string identifying which empirical estimate to used for the empirical binned covariance matrix. See details.
#' @param coords the matrix of observed locations.
#' @param Y the observed data.
#' @param X the fixed-effects design matrix.
#' @param S the matrix of basis functions.
#' @param bin_centers matrix of coordinates for the bin centers for the MOM estimation.
#' @param ... space for additional arguments to be passed to the estimation function.
#' 
#' @details
#' The Method of Moments (MOM, \code{method="MOM"}) estimation is a two-step procedure. First, an empirical binned
#' covariance matrix is computed, and then the covariance matrix is fitted (see \code{\link{mom_fit}}). Different options for the empirical 
#' estimate can be selected using \code{empirical}. Currently there are four choices:
#' \itemize{
#' \item \code{"cj"}{The Cressie-Johanneson empirical binned covariance matrix.}
#' \item \code{"robust"}{A robust empirical binned estimator based on the median absolute deviation. }
#' \item \code{"dispersion"}{A robust empirical binned estimator based on a more efficient robust variance. }
#' \item \code{"wt_dispersion"}{A weighted version of the dispersion estimate. }
#' }
#' 
#' @note
#' Alone, the output of this function is not useful, it should be used in tandem with \code{\link{mom_fit}}
#' 
#' @return
#' A matrix.
#' 
#' @seealso
#' \code{\link{rr_est}}
#' \code{\link{mom_fit}}
#' 
#' 
#' @importFrom MASS ginv
#' @importFrom fields cover.design
#' @importFrom Rfit pairup
#' 
#' @export
#' 


mom_bec <- function( empirical="cj", coords, Y, X, S, bin_centers=NULL, ...){
  
  cc <- match.call( expand.dots=TRUE )  
  
  nd <- nrow(S)
  nk <- ncol(S)
    
  bhat <- ginv( t(X)%*%X ) %*% ( t(X)%*%Y )
  R    <- Y - X%*%bhat
    
  ##########################################
  
  ## Make the boxes
  
  ## Are there bin_centers provided?
  if( !is.null(bin_centers) ){
    return_bins <- FALSE
    if( nk >= nrow(bin_centers) ){
      stop("Number of specified bins is less than number of knots")
    }
  } else{    
    ## Generate a plausible number of bin centers, hopefully approx. equadistant
    return_bins <- TRUE
    bin_centers <- cover.design( R=coords, nd=2*nk, ... )$design
  }
  
  
  
  ## boxes[,8] === NUMBER of box
  
  ## Check to make sure there are at least 15 obs / box  
  ## Use Voronii tesslation based on bin_centers
  nbox    <- nrow(bin_centers)  
  boxFreq <- cbind( 1:nrow(bin_centers), rep(0,nbox) )
  #kd_dist <- pairwise_dist( coords, bin_centers, ... ) # <----------
  kd_dist <- pairwise_dist( coords, bin_centers  )
  vBins   <- apply( kd_dist, 1 , 'which.min' ) 
  for( ii in 1:nbox){
    jj <- boxFreq[ii,1]
    boxFreq[ii,2] <- sum( vBins==jj)
  }
  binMin <- min( boxFreq[,2] )
  
  while( binMin < 15 ){
    id_bin <- which.min( boxFreq[,2] )
    boxFreq[id_bin,]
    
    ## Just remove the bin with n_in_bin < 15
    #bin_centers <- bin_centers[-id_bin,]        
    
    ## Merge bin with next closest bin    
    # kk_dist  <- pairwise_dist( bin_centers[-id_bin,,drop=FALSE], bin_centers[id_bin, ,drop=FALSE], ... )
    # id_merge <- which.min(kk_dist)
    # new_bin  <- colMeans( bin_centers[c(id_bin,id_merge),] )
    # bin_centers <- rbind(bin_centers[-c(id_bin,id_merge),], new_bin)
    
    # dis_bin  = bin to discard
    # rem_bins = bins remaining after removing discarded bin
    dis_bin  <- bin_centers[ id_bin, ,drop=FALSE]
    rem_bins <-  bin_centers[-id_bin,,drop=FALSE]
    kk_dist  <- pairwise_dist( rem_bins, dis_bin )
    id_merge <- which.min(kk_dist)
    new_bin  <- colMeans( rbind(rem_bins[id_merge,], dis_bin ) )
    bin_centers <- rbind( rem_bins[-id_merge,], new_bin)
    
    
    nbox    <- nrow(bin_centers)  
    #kd_dist <- pairwise_dist( coords, bin_centers, ... ) # <----------
    kd_dist <- pairwise_dist( coords, bin_centers )
    vBins   <- apply( kd_dist, 1 , 'which.min' ) 
    boxFreq <- cbind( 1:nrow(bin_centers), rep(0,nbox) )
    for( ii in 1:nbox){
      jj <- boxFreq[ii,1]
      boxFreq[ii,2] <- sum( vBins==jj)
    }
    binMin <- min( boxFreq[,2] )
    
  }  
  
  if( nbox<=nk ){
    stop("After ensuring n>15 locations / bin, there are fewer bins than knots")
  }
  rownames(bin_centers) <- NULL
  
  
  ### Set up the empirical binned covariance matrix
  # off_diag  is to calculate the off-diag element
  # on_diag is to calculate the diagonal element
  off_diag <- rep(0,nbox)
  on_diag  <- rep(0,nbox)
  
  if( empirical=="cj" ){
    ## The CRESSIE-JOHANNESON estimate
    ## Diagonals: mean{r_i^2}
    ## Off-diagonals: mean{r_i}*mean{r_i}    
    for(ii in 1:nbox){
      off_diag[ii] <- mean( R[ vBins==boxFreq[ii] ] )
      on_diag[ii]  <- mean( R[ vBins==boxFreq[ii] ] * R[ vBins==boxFreq[ii] ] )
    }
    SigM <- off_diag %*% t(off_diag)
    diag(SigM) <- on_diag
    
  } else if( empirical=="robust"){
    ## The ROBUST estimate
    ## Diagonals: mad^2{r_i}
    ## Off-diagonals: (mad^2{r_i + r_j} - mad^2{r_i - r_j}) / 4
    SigM <- matrix( 0 , nbox, nbox )
    for(ii in 1:nbox){
      on_diag[ii] <- mad( R[ vBins==boxFreq[ii] ] )^2
    }
    for(ii in 1:(nbox-1)){
      for( jj in (ii+1):nbox){
        # SD <- pSD( R[ vBins==boxFreq[ii] ] , R[ vBins==boxFreq[jj] ] , ... )  # <----------
        SD <- pSD( R[ vBins==boxFreq[ii] ] , R[ vBins==boxFreq[jj] ]  )
        PS11 <- SD$S
        PD11 <- SD$D
        SigM[ii,jj] <- ( mad(PS11)^2 - mad(PD11)^2 )/4
      }
    }
    
    SigM <- SigM + t(SigM)
    diag(SigM) <- on_diag
    
  } else if( empirical=="dispersion"){
    ## The GENERALIZED ROBUST estimate
    ## Diagonals:  V(D+(r_i))^2 
    ## Off-diagonals: (mad^2{r_i + r_j} - mad^2{r_i - r_j}) / 4
    SigM <- matrix( 0 , nbox, nbox )
    for(ii in 1:nbox){
      on_diag[ii] <- sig2dplus( R[ vBins==boxFreq[ii] ] )
    }
    for(ii in 1:(nbox-1)){
      for( jj in (ii+1):nbox){
        SD <- pSD( R[ vBins==boxFreq[ii] ] , R[ vBins==boxFreq[jj] ] , ... )
        PS11 <- SD$S
        PD11 <- SD$D
        SigM[ii,jj]  <- ( sig2dplus(PS11) - sig2dplus(PD11) ) / 4
      }
    }
    
    SigM <- SigM + t(SigM)
    diag(SigM) <- on_diag
    
  } else if( empirical=="wt_dispersion"){
    ## The GENERALIZED ROBUST estimate
    ## Diagonals:  V(D+(r_i))^2 
    ## Off-diagonals: (mad^2{r_i + r_j} - mad^2{r_i - r_j}) / 4
    SigM <- matrix( 0 , nbox, nbox )
    for(ii in 1:nbox){
      on_diag[ii] <- sig2dpluswts( R[ vBins==boxFreq[ii] ] )
    }
    
    for(ii in 1:(nbox-1)){
      for( jj in (ii+1):nbox){
        SD <- pSD( R[ vBins==boxFreq[ii] ] , R[ vBins==boxFreq[jj] ] , ... )
        PS11 <- SD$S
        PD11 <- SD$D
        SigM[ii,jj]  <- ( sig2dpluswts(PS11) - sig2dpluswts(PD11) ) / 4
      }
    }
    
    SigM <- SigM + t(SigM)
    diag(SigM) <- on_diag
    
  }
  
  
  ## Calculate the binned S matrix
  # M    <- ncol( SigM )
  Sbar <- matrix( 0 , nbox , nk)
  for( ii in 1:nbox){
    #Sbar[ii,] <- apply( S[ vBins==boxFreq[i] , ] , 2 , mean )
    Sbar[ii,] <- colMeans( S[ vBins==boxFreq[ii] , ] )
  }
  
  ## Make the return object
  return_obj <- list( SigM=SigM, Sbar=Sbar )
  if( return_bins ){ return_obj$bin_centers <- bin_centers }
  
  
  return( return_obj )
  
}





