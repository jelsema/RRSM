#' Fit the reduced dimension covariance matrix
#' 
#' @description
#' Uses the binned empirical covariance matrix for Method of Moments (MOM) 
#' estimation to fit the reduced dimension covariance matrix of a reduced rank spatial model.
#' 
#' @param fitting string describing which fitting method to be used to fit the covariance matrix. See details.
#' @param SigM the binned empirical covariance matrix.
#' @param Sbar a similarly binned version of the matrix of basis functions.
#' @param ... space for additional arguments to be passed to the fitting function.
#' 
#' @details
#' The Method of Moments (MOM, \code{method="MOM"}) estimation is a two-step procedure. 
#' First, an empirical binned covariance matrix is computed (see \code{\link{mom_bec}}), 
#' and then the covariance matrix is fitted. Different options for the fitting stage
#' can be selected using \code{fitting}. Currently there are two choices:
#' \itemize{
#' \item \code{"frobenius"}{The covariance is fitted by minimizing the Frobenius norm (REF).}
#' \item \code{"robust"}{The covariance is fitted by minimizing a robust norm (using quantile regression) (REF). }
#' }
#' 
#' 
#' @note
#' The values of \code{SigM} and \code{Sbar} are intended to be the output values of 
#' 
#' @return
#' A matrix.
#' 
#' @seealso
#' \code{\link{rr_est}}
#' \code{\link{mom_bec}}
#'  
#' @importFrom quantreg rq.fit.br
#' 
#' @export
#' 



mom_fit <- function( fitting="frobenius", SigM, Sbar , ...){
    
  nk <- ncol(Sbar)  # Number of knots
  nb <- ncol(SigM)  # Number of bins
  
  if( fitting=="frobenius" ){
    ## The LEAST-SQUARES / FROBENIUS Fitting method
    
    # Get the QR decomposition of Sbar
    Q <- qr.Q( qr(Sbar))
    R <- qr.R( qr(Sbar))
    
    E   <- Q%*%t(Q)
    L   <- SigM - E%*%SigM%*%E
    L1  <- L[1:nb, 1:nb]
    LL1 <- c(L1)
    
    Pk  <- diag(nb) - E
    Pk2 <- c(Pk)
    
    ssq <- ginv(t(Pk2)%*%Pk2) %*% (t(Pk2)%*%LL1)
    ssq <- ssq[1,1]
    
    # R1in <- solve(R) %*% diag(nk)
    # RRin <- R1in
    RRin <- solve(R) %*% diag(nk)
    
    QQ <- Q
    
    Dbar <- ssq * diag(nb)
    # Dbarhin <- 1/ss1*diag(nb)
    Vhat <- RRin%*%t(QQ)%*%(SigM - Dbar)%*%QQ%*%t(RRin)
        
  }
  
  if( fitting=="robust" ){
    ## The QUANTILE REGRESSION / ROBUST Fitting method
    
    SbSinv  <- solve( t(Sbar) %*% Sbar )
    SigProj <- diag(nb) - Sbar %*% SbSinv %*% t(Sbar)

    # Estimate sigma^2
    sigY <- c( SigM %*% SigProj )
    sigX <- c( SigProj )
    sig  <- rq.fit.br( sigX , sigY , tau=0.5, alpha=0.1, ci=FALSE )
    ssq  <- sig$coef
    
    SigSSbS <- (SigM - ssq*diag(nb)) %*% Sbar %*% SbSinv
    
    Vhat    <- matrix( 0 , nrow=nk , ncol=nk )
    for( ii in 1:nk ){
      test <- rq.fit.br( Sbar , SigSSbS[,ii] , tau=0.5, alpha=0.1, ci=FALSE )
      Vhat[,ii] <- test$coef
    }
    
    Vhat <- 0.5*( Vhat + t(Vhat) )
    
  }
  
  ## LIFT THE EIGENVALUES
  ## Look at the automatic way to do this described elsewhere (Kang? Katzfuss?)
  SVD    <- eigen(Vhat)
  EigVec <- SVD$vectors
  EigVal <- SVD$values
  
  if( min(EigVal) < 0 ){
    EigValSum <- sum(EigVal)                              # Sum of Eignevalues, to be partitioned
    EigVal <- EigVal - min(EigVal)                        # Shift Eigenvalues if (probably) necessary
    NewEigVals <- EigValSum * (EigVal/sum(EigVal))        # Compute new eignevalues by weighting shifted eigenvalues
    Vhat <- EigVec %*% diag(NewEigVals) %*% solve(EigVec) # Reconstitute Vhat
  }
  
  ## Make the return object
  return_obj <- list( bhat=NULL, V=Vhat, ssq=ssq )
  return( return_obj )
  
}




