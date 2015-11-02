#' Estimate the parameters of a spatial mixed effects model
#' 
#' @description
#' Computes estimates of parameters for spatial mixed effects models using an EM algorithm.
#' 
#' @rdname rr_em
#' 
#' @param Y the observed data.
#' @param X the fixed-effects design matrix.
#' @param S the matrix of basis functions.
#' @param epsilon the stopping condition for EM algorithms.
#' @param max_iter for EM algorithms, the maximum number of iterations.
#' @param ... space for additional arguments to be passed to the estimation function.
#' 
#' @note
#' Some of the estimation functions do not estimate the fixed effects parameter (e.g., 
#' \code{rr_em_nobeta}). In this case, \code{bhat} in the return object is \code{NULL}. 
#' In the estimation functions of \pkg{RRSM}, if \code{bhat} is \code{NULL}, then a 
#' GLS estimate will be computed automatically.
#' 
#' @seealso
#' \code{\link{rr_est}}
#' \code{\link{mom_bec}}
#' \code{\link{mom_fit}}
#' \code{\link{rr_universal_krige}}
#' 
#' @export
#' 


rr_em <- function( Y, X, S, epsilon=0.00001 , max_iter=500, ... ){
  
  nd <- nrow(S)
  nk <- ncol(S)
  
  ## Big outside computations
  SS <- t(S)%*%S
  XX <- t(X)%*%X
  SX <- t(S)%*%X
  SY <- t(S)%*%Y
  XY <- t(X)%*%Y
  
  XXi <- ginv(XX)
  
  ## Initial values
  V    <- diag(nk)
  bhat <- ginv(XX) %*% XY
  Rt   <- (Y - X%*%bhat)
  ssq  <- c(var(Rt))
  SdvS <- ginv( ( ginv(V) + SS/ssq ) )
  
  DONE_EST <- iter <- 0
  
  while( DONE_EST==0 ){
    iter  <- iter+1
    
    V1    <- V
    bhat1 <- bhat
    ssq1  <- ssq
    
    Rt <- (Y - X%*%bhat)
    RR <- sum( Rt^2 )
    SR <- t(S)%*%Rt
    
    ## Estimate ssq
    mean_nu <- SdvS %*% (SR/ssq)    
    ssq     <- (1/nd)*c(RR - 2*t(SR)%*%mean_nu + sum(diag(SS%*%SdvS)) + t(mean_nu)%*%SS%*%mean_nu)  
    
    ## Estimate V    
    V    <- SdvS + mean_nu %*% t(mean_nu)
    SdvS <- ginv( ginv(V) + (SS/ssq) )
    
    ## Estimate beta
    XPX  <- XX/ssq - t(SX)%*%SdvS%*%SX
    bhat <- (XXi %*% t(SX) %*% mean_nu) + (XXi %*% XY)
    
    ## Assess convergence    
    v_conv <- (sqrt(sum(sum((V-V1)^2)))/(nk^2) < epsilon)
    b_conv <- (sqrt(sum(sum((bhat-bhat1)^2)))/(nk^2) < epsilon)
    
    if( ((v_conv & b_conv) | (iter >= max_iter)) ){ DONE_EST <- 1 }
  }
  
  
  ## Create the rr_est object
  results <- list( bhat=bhat, V=V, ssq=ssq, niter=iter )
    
  return( results )
  
}
