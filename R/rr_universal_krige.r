#' Universal kriging predictions for spatial process
#' 
#' @description
#' Obtains spatial predictions at selected locations.
#' 
#' @inheritParams rr_predict
#' @inheritParams rr_lognormal_krige
#' 
#' @rdname rr_krige
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export
#' 
#' 


rr_universal_krige <- function( Y, X, S=NULL, coords, pgrid=coords, Xpred=X , 
                                Spred=S, V=NULL, ssq, tsq=0, ncore=1, ... ){
  
  cc <- match.call( expand.dots=TRUE )  
  
  # If S is not provided then it must be calculated
  if( is.null(S) ){
    
    ## Extract knots from the "..."
    
    if( is.null(coords) || is.null(knots) ){
      stop("Must provide either S-matrix or {locations , knots}") 
    } else{
      S     <- basis_mat( coords , knots , mult=1.3)$S    
      Spred <- basis_mat( pgrid , knots , mult=1.3)$S  
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
  
  
  
  ## Write logic to properly set up sigma^2 vs tau^2
  ngrid  <- dim(pgrid)[1]
  nk     <- dim(S)[2]
  Preds  <- cbind( pgrid , rep(0,ngrid) , rep(0,ngrid) )
  
  ## Big outside computations
  SS <- t(S)%*%S
  XX <- t(X)%*%X
  SX <- t(S)%*%X
  SY <- t(S)%*%Y
  XY <- t(X)%*%Y
  
  SdvS <- ginv( ginv(V) + (SS/ssq) )
  
  XPX  <- XX/ssq - (t(SX)/ssq) %*% SdvS %*% (SX/ssq)
  XPY  <- XY/ssq - (t(SX)/ssq) %*% SdvS %*% (SY/ssq)
  bhat <- ginv(XPX)%*%XPY

  R  <- Y - X%*%bhat
  SR <- t(S)%*%R  
  
  ## START OF KRIGING
  if( round(ncore)==ncore & ncore > 1 ){
    ## SIMULTANEOUS PREDICTION (PARALLEL PROCESSING)
    registerDoParallel(ncore)
    
    Preds1 <- foreach( krg = 1:ngrid, .combine='rbind', .packages="foreach" ) %dopar% {
      
      SkSo <- Spred[krg,]
      
      # Get the proper Cwk(So)
      Ck <- SkSo %*% V %*% t(S);
      
      # Get the Xso vector (using Legendre polynomials)
      XSo1 <- Xpred[krg,]
      
      # Prediction: Preds[krg,3]
      pred <- XSo1%*%bhat + (Ck%*%R)/ssq - (Ck%*%S/ssq)%*%( SdvS %*% (SR/ssq) ) ;
      
      # Calculate MSPE
      CSoSo <- (SkSo %*% V %*% SkSo) + tsq;
      
      SdC <- ( t(S)%*%t(Ck) )/ssq ;
      CdC <- (Ck%*%t(Ck))/ssq - (t(SdC)%*%(SdvS %*% SdC))
      
      XdC <- (t(X)%*%t(Ck))/ssq - (t(SX/ssq)%*%(SdvS%*%SdC))
      R1  <- XSo1 - XdC
      
      mspe <- sqrt(c( CSoSo - CdC + t(R1)%*% ginv(XPX) %*% R1 ))
      
      # MSPE: Preds[krg,4]
      c(pred, mspe)
      
    }
    Preds[,3:4] <- Preds1
    
  } else{
    ## SEQUENTIAL PREDICTIONS
    
    for( krg in 1:ngrid ){
      
      SkSo <- Spred[krg,]
      
      # Get the proper Cwk(So)
      Ck <- SkSo %*% V %*% t(S);
      
      # Get the Xso vector (using Legendre polynomials)
      XSo1 <- Xpred[krg,]
      
      # Prediction: Preds[krg,3]
      pred <- XSo1%*%bhat + (Ck%*%R)/ssq - (Ck%*%S/ssq)%*%( SdvS %*% (SR/ssq) ) ;
      
      # Calculate MSPE
      CSoSo <- (SkSo %*% V %*% SkSo) + tsq;
      
      SdC <- ( t(S)%*%t(Ck) )/ssq ;
      CdC <- (Ck%*%t(Ck))/ssq - (t(SdC)%*%(SdvS %*% SdC))
      
      XdC <- (t(X)%*%t(Ck))/ssq - (t(SX/ssq)%*%(SdvS%*%SdC))
      R1  <- XSo1 - XdC
      
      mspe <- sqrt( c( CSoSo - CdC + t(R1)%*% ginv(XPX) %*% R1 ) )
      
      # MSPE: Preds[krg,4]
      Preds[krg,3] <- pred
      Preds[krg,4] <- mspe
      
    }    
    
  } ## END OF KRIGING
  colnames(Preds) <- c( "Xdim" , "Ydim" , "Y" , "RMSPE" )
  rownames(Preds) <- NULL
    
  
  ## Return the results  
  return_obj <- list( Preds=Preds , bhat=bhat )
  return(return_obj)
}
