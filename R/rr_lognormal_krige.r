#' Spatial predictions for lognormal process
#' 
#' @description
#' Obtains spatial predictions at selected locations for a lognormal process.
#' 
#' @inheritParams rr_predict
#' @param Y the observed data. Should already be on log-scale.
#' @param X the design matrix for the large-scale variation (fixed effects).
#' @param S the matrix of basis functions (design matrix for the random effect).
#' @param coords the matrix of observed locations.
#' @param pgrid the prediction grid.
#' @param Xpred the design matrix for the predicted locations.
#' @param Spred the matrix of basis functions for the predicted locations.
#' @param V   the estimated covariance matrix.
#' @param ssq the estimated residual variance.
#' @param tsq the estimated variance due to dimension reduction.
#' @param ncore if \code{ncore>1}, number of cores to use with \pkg{doParallel} backend.
#' @param ... space for additional arguments.
#' 
#' @return
#' A list containing [[1]] a matrix with the coordinates of the predicted locations, the 
#' predicted value, and measures of uncertainty, and [[2]] the estimate of large-scale variation.
#' 
#' @note
#' The parameter estimates may be left empty, and \code{...} may be used to pass arguments to estimation functions.
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export
#' 

## Need to clean up the parameters, make it more efficient to run the code

rr_lognormal_krige <- function( Y, X, S=NULL, coords, pgrid=coords, Xpred=X , 
                                Spred=S, V=NULL, ssq, tsq=0, ncore=1, ... ){
  
  cc <- match.call( expand.dots=TRUE )  
  
  # If S is not provided then it must be calculated
  if( is.null(S) ){
    
    ## Extract knots from the "..."
    
    if( is.null(coords) || is.null(knots) ){
      stop("Must provide either S-matrix or {locations , knots}") 
    } else{
      S     <- basis_mat( coords , knots , a.mult=1.3)$S    
      Spred <- basis_mat( pgrid , knots , a.mult=1.3)$S  
    }
  }
    
  # If V is not provided then it must be estimated
  if( is.null(V) ){
    
    method <- cc$method
    # epsilon
    # max_iter
    
    if( is.null(method) ){
      #warning("No estimation method provided, attempting EM Algorithm") 
    }
    
    v_est <- rr_est( method="EM", em_type="default", Y=Y, X=X, S=S )
    V     <- v_est$V
    ssq   <- v_est$ssq    
  }
  
  # Write logic to properly set up sigma^2 vs tau^2
  
  ngrid  <- nrow(pgrid)
  nd     <- nrow(S)
  nk     <- ncol(S)
  Preds  <- cbind( pgrid , rep(0,ngrid) , rep(0,ngrid) )
  
  ## Calculate constants (especially large computations)
  SS <- t(S)%*%S
  XX <- t(X)%*%X
  SX <- t(S)%*%X
  SY <- t(S)%*%Y
  XY <- t(X)%*%Y
  
  SdvS <- ginv( ginv(V) + (SS/ssq) )
  
  XPX  <- XX/ssq - (t(SX)/ssq) %*% SdvS %*% (SX/ssq)
  XPY  <- XY/ssq - (t(SX)/ssq) %*% SdvS %*% (SY/ssq)
  bhat <- ginv(XPX)%*%XPY

  R    <- Y - X%*%bhat
  SR   <- t(S)%*%R  
  
  XB   <- X %*% bhat
  SXB  <- t(S)%*%XB
  
  XPXB <- t(X)%*%XB/ssq - (t(SX)/ssq) %*% SdvS %*% (SXB/ssq)
  
  
  ## START OF KRIGING
  if( round(ncore)==ncore & ncore > 1 ){
    ## SIMULTANEOUS PREDICTION (PARALLEL PROCESSING)
    registerDoParallel(ncore)
    
    Preds1 <- foreach( krg = 1:ngrid, .combine='rbind', .packages="foreach") %dopar% {
      
      SkSo <- Spred[krg,]
      
      # Get the proper Cwk(So)
      Ck <- SkSo %*% V %*% t(S)
      
      # Get the Xso vector
      XSo1 <- Xpred[krg,]
      
      
      CSoSo <- c( (SkSo %*% V %*% SkSo) + tsq )
      SdC   <- ( t(S)%*%t(Ck) )/ssq 
      CPC   <- (Ck%*%t(Ck))/ssq - (t(SdC)%*%(SdvS %*% SdC))    
      XPC   <- (t(X)%*%t(Ck))/ssq - (t(SX/ssq)%*%(SdvS%*%SdC))
      YPC   <- (t(Y)%*%t(Ck))/ssq - (t(SY/ssq)%*%(SdvS%*%SdC))
      
      AkY <- c( ((XSo1 - t(XPC)) %*% bhat) + t(YPC) )
      Bk  <- c( 0.5*(CSoSo - CPC) )
      
      pred <- exp( AkY + Bk )
      
      xb   <- c( XSo1 %*% bhat )
      XBPC <- (t(XB)%*%t(Ck))/ssq - (t(SXB/ssq)%*%(SdvS%*%SdC))
      AkMU <- c( (XSo1 - t(XPC)) %*% ginv(XPX)  %*% XPXB + t(XBPC) )
      
      AkPAk <- c( ((XSo1 - t(XPC)) %*% ginv(XPX) %*% t(XSo1 - t(XPC))) + 
                    (2*(XSo1 - t(XPC)) %*% ginv(XPX) %*% XPC) + CPC )
      
      AC    <- c(((XSo1 - t(XPC)) %*% ginv(XPX)  %*% XPC) + CPC)
      
      
      mspe <- (exp( 2*(xb + CSoSo) ) + exp( 2*(AkMU + AkPAk + Bk)  ) - 
                 2*exp( xb + 0.5*CSoSo + AkMU + 0.5*AkPAk + AC + Bk ))
      
      c(pred, mspe)
      
    }
    
    Preds[,3:4] <- Preds1
    
  } else{
    ## SEQUENTIAL PREDICTIONS
    
    for( krg in 1:ngrid ){
      
      SkSo <- Spred[krg,]
      
      # Get the proper Cwk(So)
      Ck <- SkSo %*% V %*% t(S)
      
      # Get the Xso vector
      XSo1 <- Xpred[krg,]
      
      
      CSoSo <- c( (SkSo %*% V %*% SkSo) + tsq )
      SdC   <- ( t(S)%*%t(Ck) )/ssq 
      CPC   <- (Ck%*%t(Ck))/ssq - (t(SdC)%*%(SdvS %*% SdC))    
      XPC   <- (t(X)%*%t(Ck))/ssq - (t(SX/ssq)%*%(SdvS%*%SdC))
      YPC   <- (t(Y)%*%t(Ck))/ssq - (t(SY/ssq)%*%(SdvS%*%SdC))
      
      AkY <- c( ((XSo1 - t(XPC)) %*% bhat) + t(YPC) )
      Bk  <- c( 0.5*(CSoSo - CPC) )
      
      pred <- exp( AkY + Bk )
      
      xb   <- c( XSo1 %*% bhat )
      XBPC <- (t(XB)%*%t(Ck))/ssq - (t(SXB/ssq)%*%(SdvS%*%SdC))
      AkMU <- c( (XSo1 - t(XPC)) %*% ginv(XPX)  %*% XPXB + t(XBPC) )
      
      AkPAk <- c( ((XSo1 - t(XPC)) %*% ginv(XPX) %*% t(XSo1 - t(XPC))) + 
                    (2*(XSo1 - t(XPC)) %*% ginv(XPX) %*% XPC) + CPC )
      
      AC    <- c(((XSo1 - t(XPC)) %*% ginv(XPX)  %*% XPC) + CPC)
      
      
      mspe <- (exp( 2*(xb + CSoSo) ) + exp( 2*(AkMU + AkPAk + Bk)  ) - 
                 2*exp( xb + 0.5*CSoSo + AkMU + 0.5*AkPAk + AC + Bk ))
      
      Preds[krg,3] <- pred
      Preds[krg,4] <- mspe
      
    }
  }  ## END OF KRIGING
  
  colnames(Preds) <- c( "Xdim" , "Ydim" , "Y" , "MSPE" )
  rownames(Preds) <- NULL    
  
  ## Return the results  
  return_obj <- list( Preds=Preds , bhat=bhat )
  return( return_obj )
}







