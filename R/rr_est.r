#' Estimate the parameters of a spatial mixed effects model
#' 
#' @description
#' Computes estimates of parameters for spatial mixed effects models using either EM algorithm or
#' a Method of Moments estimation.
#' 
#' @param method string describing the type of parameter estimation; currently 
#'        supported are Expectation-Maximization ("EM") and Method of Moments ("MOM").
#' @param em_type string describing further specification for the EM algorithm. See details.
#' @param empirical string identifying which empirical estimate to used for the empirical binned covariance matrix. See details.
#' @param fitting   string identifying which fitting method to use for MOM estimation. See details.
#' @param Y the observed data.
#' @param X the fixed-effects design matrix.
#' @param S the matrix of basis functions.
#' @param X2 for constrained EM algorithm, optional design matrix for unconstrained fixed-effects.
#' @param cone for constrained EM algorithm, matrix of -1, 0, and 1's defining the order restriction.
#' @param coords the matrix of observed locations.
#' @param knots the matrix of knot locations.
#' @param epsilon the stopping condition for EM algorithms.
#' @param max_iter for EM algorithms, the maximum number of iterations.
#' @param ... space for additional arguments to be passed to the estimation function (see especially \code{rr_em_constrained}).
#' 
#' @details
#' \code{method} is used to define a general type of estimation method.
#' For \code{method="EM"}, the supported options for \code{em_type} are:
#' \itemize{
#' \item \code{"default"}{A standard EM algorithm to estimate all model parameters.}
#' \item \code{"nobeta"}{Similar to the default method, but the fixed-effects parameter is treated not estimated.}
#' \item \code{"constrained"}{An EM algorithm in which the fixed effects parameter is subject to linear order constraints.}
#' }
#' The Method of Moments (MOM, \code{method="MOM"}) estimation is a two-step procedure. 
#' First, an empirical binned covariance matrix is computed (see \code{\link{mom_bec}}), 
#' and then the covariance matrix is fitted (see \code{\link{mom_fit}}. Different options 
#' for the empirical estimate can be selected using \code{empirical}. Currently the options are:
#' \itemize{
#' \item \code{"cj"}{The Cressie-Johanneson empirical binned covariance matrix.}
#' \item \code{"robust"}{A robust empirical binned estimator based on the median absolute deviation. }
#' \item \code{"dispersion"}{A robust empirical binned estimator based on a more efficient robust variance (computationally intensive). }
#' \item \code{"wt_dispersion"}{A weighted version of the dispersion estimate (computationally intensive). }
#' }
#' Additionally, there are multiple options for the fitting method (\code{fitting}):
#' \itemize{
#' \item \code{"frobenius"}{The covariance is fitted by minimizing the Frobenius norm.}
#' \item \code{"robust"}{The covariance is fitted by minimizing a robust norm (using quantile regression). }
#' }
#'  
#' @note
#' Some of the estimation functions do not estimate the fixed effects parameter (e.g., 
#' \code{rr_em_nobeta}, and several of the MOM estimation functions). In this case, \code{bhat} 
#' in the return object is \code{NULL}. In the estimation functions of \pkg{RRSM}, if \code{bhat}
#' is \code{NULL}, then a GLS estimate will be computed automatically.
#' 
#' @return
#' A list with elements \code{bhat} (estimate of beta), \code{V}, \code{ssq}, \code{niter}, 
#' and \code{method}.
#' 
#' @seealso
#' \code{\link{rr_em}}
#' \code{\link{mom_bec}}
#' \code{\link{mom_fit}}
#' \code{\link{rr_universal_krige}}
#' 
#' @export
#' 


##
################# ADD ALL THE ARGUMENTS NECESSARY -- DATA, ETC
##

rr_est <- function( method="EM", em_type="default", empirical="cj", fitting="frobenius",
                    Y, X, S, X2=NULL, cone=NULL, coords=NULL, knots=NULL, ... ){
  
  cc <- match.call( expand.dots=TRUE )  
  
  ## Some error-catching
  if( missing(Y) ){
    stop("The data vector 'Y' is missing")
  }
  if( missing(X) ){
    stop("The design matrix 'X' is missing")
  }
  if( missing(S) ){
    if( is.null(knots) || is.null(coords) ){
      stop("Either matrix of basis functions 'S', or both 'coords' and 'knots' are required")
    }
    warning("Matrix of basis functions 'S' was not provided. Constructing 'S' based on 'coords' and 'knots', but unsupervised construction of 'S' is not recommended.")
    S <- basis_mat( coords , knots , ... )    
  }
  
  
  ## Capitalize the text parameters
  method   <- toupper(method)
  matchMth <- charmatch( method, c("EM","MOM") , nomatch = 0)
  if( matchMth==0 ){
    stop( "Argument 'method' does not match a supported value ('EM' or 'MOM')" )
  }
  
  ## Start the estimation procedures
  if( method=="EM" ){
    ## Some error-fixing for the function calls
    em_type <- tolower(em_type)
    all_types <- c("default", "nobeta", "constrained")
    matchType <- charmatch( em_type, all_types, nomatch = 0)
    
    if( matchType > 0 ){ em_type <- all_types[matchType]
    } else{ stop("Fitting method cannot be determined") }
    
    if( em_type=="default" ){
      return_obj <- rr_em(Y=Y, X=X, S=S, ... )
    } else if( em_type=="nobeta" ){
      return_obj <- rr_em_nobeta(Y=Y, X=X, S=S, ... )
    } else if( em_type=="constrained" ){
      return_obj <- rr_em_constrained(Y=Y, X=X, X2=X2, S=S, cone=cone, ...)
    } else{
      stop( "argument 'em_type' does not match a supported value (default, nobeta, constrained)" )
    }
    
    full_method         <- list()
    full_method$type    <- "Estimation type: EM"
    full_method$em_type <- paste0( "EM method: ", em_type )  
    full_method$m_iter  <- paste0( "EM algorithm finished after ", return_obj$niter, " iterations" )
    
  } else if( method=="MOM" ){
    
    ## Some error-fixing for the function calls
    all_emps <- c("cj", "robust", "dispersion", "wt_dispersion" )
    empirical <- tolower(empirical)
    matchEmp <- charmatch( empirical, all_emps, nomatch = 0)
    if( matchEmp > 0 ){ empirical <- all_emps[matchEmp]
    } else{ stop("Fitting method cannot be determined") }
    
    all_fits <- c("frobenius", "robust")
    fitting   <- tolower(fitting)
    matchFit <- charmatch( fitting, all_fits, nomatch = 0)
    if( matchFit > 0 ){ fitting <- all_fits[matchFit]
    } else{ stop("Fitting method cannot be determined") }
    
    ## Get the Empirical binned covariance matrix
    if( empirical %in% all_emps ){
      emp_cov <- mom_bec( empirical=empirical, coords=coords, Y=Y, X=X, S=S, ... )  
    }
    
    ## Fit the covariance matrix
    return_obj <- mom_fit( fitting=fitting, SigM=emp_cov$SigM, Sbar=emp_cov$Sbar, ... )
        
    full_method           <- list()
    full_method$type      <- "Estimation type: MOM"
    full_method$fitting   <- paste0( "Fitting method: ", fitting )
    full_method$empirical <- paste0( "Empirical binned covariance matrix: ", empirical )
    
  }
  
  ## Estimate tau vs. sigma
  
  
  
  ## Add some values to the output object
  return_obj$full_method <- full_method
  if( is.null(return_obj$niter) ){ return_obj$niter <- NULL }
  
  class(return_obj) <- "rr_est"
  
  return( return_obj )
  
}
