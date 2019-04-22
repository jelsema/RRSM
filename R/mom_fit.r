#' Fit the reduced dimension covariance matrix
#' 
#' @description
#' Uses the binned empirical covariance matrix for Method of Moments (MOM) 
#' estimation to fit the reduced dimension covariance matrix of a reduced rank spatial model.
#' 
#' @param fitting string describing which fitting method to be used to fit the covariance matrix. See details.
#' @param SigM the binned empirical covariance matrix.
#' @param Sbar a similarly binned version of the matrix of basis functions.
#' @param mcn if positive, winsorizes eigenvalues of covariance matrix so that the condition number is at max \code{mcn}.
#' @param vlift if \code{TRUE}, lifts eigenvalues of robust covariance matrix ("proportional").
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
#' The values of \code{SigM} and \code{Sbar} are intended to be the output values of \code{\link{mom_bec}}
#' 
#' @return
#' A matrix.
#' 
#' @seealso
#' \code{\link{rr_est}}
#' \code{\link{mom_bec}}
#'  
#' @importFrom quantreg rq.fit.br
#' @importFrom pracma fzero
#' 
#' @export
#' 



mom_fit <- function( fitting="frobenius", SigM, Sbar , mcn=1000, vlift="none",  ...){
    
  nk <- ncol(Sbar)  # Number of knots
  nb <- ncol(SigM)  # Number of bins
  
  if( fitting=="frobenius" ){
    ## The LEAST-SQUARES / FROBENIUS Fitting method
    
    ## LIFT THE EIGENVALUES OF [SigM] IF NECESSARY
    SVD    <- eigen(SigM)
    EigVec <- SVD$vectors
    EigVal <- SVD$values
    
    if( min(EigVal) < 0 ){
      MM <- length(EigVal)
      qq <- ncol(Sbar)
      aa <- -1
      while( aa < 0 ){
        if( qq >= 0 ){
          lam0 <- max( c( quantile(EigVal, (MM - qq)/MM), 0) )
        } else{
          lam0 <- EigVal[1] - qq
        }
        # aa <- fzero( tr_var, x=c(0,10), lam0=lam0, EigVal=EigVal )$x
        chk <- tryCatch( 
          expr    = fzero( tr_var, x=c(0,10), lam0=lam0, EigVal=EigVal )$x,
          warning = function(c) "bad_Result" ,
          error   = function(c) "bad_Result"
        )
        if( is.numeric(chk) ){
          aa <- chk
        } else{
          qq <- qq - 1
        }
      }
      
      sEigVal <- EigVal
      idx1    <- which( EigVal < lam0 )
      sEigVal[idx1] <- lam0 * exp( aa*(sEigVal[idx1]-lam0) )
      SigM <- EigVec %*% diag(sEigVal) %*% t(EigVec)
    } ## End of lifting eigenvalues
    
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
    
    ## LIFT THE EIGENVALUES OF [Vhat] IF NECESSARY
    ## MAYBE WINSORIZE THE EIGENVALUES HERE INSTEAD OF THIS SCHEME?

    if( vlift == "proportional" ){
      SVD    <- eigen(Vhat)
      EigVec <- SVD$vectors
      EigVal <- SVD$values
      
      if( min(EigVal) < 0 ){
        EigValSum <- sum(EigVal)                              # Sum of Eignevalues, to be partitioned
        EigVal <- EigVal - min(EigVal)                        # Shift Eigenvalues
        NewEigVals <- EigValSum * (EigVal/sum(EigVal))        # Compute new eignevalues by weighting shifted eigenvalues
        Vhat <- EigVec %*% diag(NewEigVals) %*% solve(EigVec) # Reconstitute Vhat
      } 
      
    }
    ## End of lifting eigenvalues
    
  }
  
  if( mcn > 0 ){
    SVD    <- eigen(Vhat)
    EigVec <- SVD$vectors
    EigVal <- SVD$values
    idx_lift <- max( which( EigVal[1] / EigVal < mcn ) )
    
    NewEigVals <- EigVal
    NewEigVals[ idx_lift:length(EigVal) ] <- EigVal[idx_lift]
    Vhat <- EigVec %*% diag(NewEigVals) %*% solve(EigVec)
  }
  
  
  ## Make the return object
  return_obj <- list( bhat=NULL, V=Vhat, ssq=ssq )
  return( return_obj )
  
}




