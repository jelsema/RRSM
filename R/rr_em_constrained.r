#' @rdname rr_em
#' 
#' @inheritParams rr_em
#' @param X2 for constrained EM algorithm, optional design matrix for unconstrained fixed-effects.
#' @param cone for constrained EM algorithm, matrix of -1, 0, and 1's defining the order restriction.
#' @param mySolver the solver to be used in the order restriction (passed to \code{isotone::activeSet}.
#' 
#' 
#' @importFrom isotone activeSet
#' @importFrom isotone gpava
#' 
#' @export
#' 


rr_em_constrained <- function( Y, X, X2=NULL, S , cone=NULL, epsilon=0.00001 , max_iter=500, mySolver="LS", ... ){
  
  X1 <- X
  P1 <- ncol(X1)
  X  <- as.matrix(cbind(X1, X2))
  nd <- nrow(S)
  nk <- ncol(S)
  pp <- ncol(X)
  
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
  
  
  if( FALSE ){
    ## Put in a check here. Look at theta = c * bhat
    
    # If very negative, give warning, and provide an option for continuing 
    # anyway (by rerunning the model with an extra argument to ignore the error)
    
    # E.g.
    # if( theta<0 & (ignore_warn()==FALSE) ){ error(...) }
    # if( theta<0 & (ignore_warn()==TRUE) ){ warning(...  EM algorithm may be unstable.) }
    
  }
  
  
  
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
    # GLS estimaton of beta (UNSTABLE, ONLY USE FOR THE WEIGHTS)
    XPX  <- XX/ssq - t(SX)%*%SdvS%*%SX
    bhat <- (XXi %*% t(SX) %*% mean_nu) + (XXi %*% XY)
    
    ## Apply order constraints / isotonization
    if( !is.null(cone) & iter>5 ){
      cov_bhat <- ginv( XPX )
      if( mySolver=="gpava" ){
        wts        <- diag(ginv(cov_bhat))[1:P1]
        bhat[1:P1] <- gpava( z= seq(1,P1), y=bhat[1:P1], weights=wts )$x
      } else{
        if( mySolver=="GLS"){
          wts <- ginv( cov_bhat )[1:P1, 1:P1, drop=FALSE]
        } else{
          wts <- diag(ginv(cov_bhat))[1:P1]
        }
        bhat[1:P1] <- activeSet( cone, y=bhat[1:P1], weights=wts, mySolver=mySolver )$x
      }
    }
    
    ## Assess convergence    
    v_conv <- (sqrt(sum(sum((V-V1)^2)))/(nk^2) < epsilon)
    b_conv <- (sqrt(sum(sum((bhat-bhat1)^2)))/(nk^2) < epsilon)
    
    if( ((v_conv & b_conv) | (iter >= max_iter)) & iter>6 ){ DONE_EST <- 1 }
  }
  
  
  ## Create the rr_est object
  results <- list( bhat=bhat, V=V, ssq=ssq, niter=iter )
    
  return( results )
  
}
