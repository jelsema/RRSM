#' @rdname rr_est
#' 
#' @inheritParams rr_em
#' 
#' @export
#' 


## Estimate V and sigma^2 --- EM Algorithm 
## REWRITE THIS FUNCTION WITH MY UPDATED NOTATION TO BE CONSISTENT WITH RR_CEST()

rr_em_nobeta <- function( Y, X, S, epsilon=0.00001, max_iter=500, ... ){
  
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
  SR   <- t(S)%*%Rt
  ssq  <- c(var(Rt))
  SdvS <- ginv( ( ginv(V) + SS/ssq ) )
  RR   <- sum( Rt^2 )
  
  DONE_EST <- iter <- 0
  
  while( DONE_EST==0 ){
    iter <- iter+1
    
    V1   <- V
    ssq1 <- ssq
    
    ## Estimate ssq
    mean_nu <- SdvS %*% (SR/ssq)    
    ssq     <- (1/nd)*c(RR - 2*t(SR)%*%mean_nu + sum(diag(SS%*%SdvS)) + t(mean_nu)%*%SS%*%mean_nu)
    
    ## Estimate V
    V    <- SdvS + mean_nu %*% t(mean_nu)
    SdvS <- ginv(( ginv(V) + SS/ssq ))
    
    ## Assess convergence
    v_conv <- (sqrt(sum(sum((V-V1)^2)))/(nk^2) < epsilon)
    
    if( (v_conv  | (iter >= max_iter)) ){ DONE_EST <- 1 }
    
  }  # End EM-algorithm
  

  ## Create the rr_est object
  results <- list( bhat=NULL, V=V , ssq=ssq , niter=iter )
    
  return( results )
  
}
