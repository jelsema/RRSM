#' Calculate despersion, \eqn{D^{+}{D^+}}
#' 
#' @description
#' Calculate the value of the dispersion function, \eqn{D^{+}{D^+}}
#' 
#' @param x a vector of numeric values.
#' 
#' @return
#' Numeric value
#' 
#' @importFrom ICSNP hl.loc
#' @export
#' 

dplus <- function(x){
  return( dplusC(x) )
}


#' Calculate variance using dispersion
#' 
#' @description
#' Estimates the variance of a sample using the dispersion function.
#' 
#' @param x a vector of numeric values.
#' 
#' @return
#' Numeric value
#' 
#' @export
#' 
sig2dplus <- function(x){
  return( sig2dplusC(x) )
}


#####################################################################3


#' Weight function for robust variance
#' 
#' @description
#' Computes weights for the weighted dispersion based variance estimate.
#' 
#' @param xin a vector of numeric values.
#' 
#' @return
#' Numeric
#' 
#' @export
#' 
wtii <- function(xin){
  x      <- xin - median(xin)
  sighat <- mad(x)
  ahat   <- x/sighat
  n      <- length(x)
  cc     <- abs(median(ahat) + 3*mad(ahat))
  wtii   <- cc/abs(ahat)
  wtii   <- 1*(wtii>=1) + wtii*(wtii<1)
  return( wtii )
}



#' Weight function for robust variance
#' 
#' @description
#' Computes weights for the weighted dispersion based variance estimate.
#' 
#' @param xin a vector of numeric values.
#' 
#' @return
#' Numeric
#' 
#' @export
#' 
wtij <- function(xin){
  x      <- xin - median(xin)
  sighat <- mad(x)
  ahat   <- x/sighat
  n      <- length(x)
  cc     <- abs(median(ahat) + 3*mad(ahat))
  xx     <- pairup(xin)
  wtij0  <- wtijloop( cc, ahat )
  wtij   <- cbind(xx, wtij0 )
  return( wtij )
}



#' Weight function for robust variance
#' 
#' @description
#' Computes weights for the weighted dispersion based variance estimate.
#' 
#' @param yin a vector of numeric values.
#' 
#' @return
#' Numeric
#' 
#' @importFrom ICSNP hl.loc
#' 
#' @export
#' 
normij <- function(yin){
  #fit <- wilcox.test(yin,conf.int=TRUE)$est
  fit   <- hl.loc(yin)
  n   <- length(yin)
  y   <- yin - fit
  
  bii   <- wtii(y)
  part1 <- sum(bii*abs(y))
  
  bij <- wtij(y)
  c1  <- bij[,1]
  c2  <- bij[,2]
  c3  <- bij[,3]
  
  part2 <- sum(c3*abs((c1+c2)/2))
  part3 <- sum(c3*abs((c1-c2)/2))
  
  normij <- (sqrt(3)/(n+1))*(part1 + part2 + part3)
  return(normij)
}


#' Weighted variance using dispersion
#' 
#' @description
#' Estimates the variance of a sample using the weighted dispersion function.
#' 
#' @param x a vector of numeric values.
#' 
#' @return
#' Numeric value
#' 
#' @export
#' 
sig2dpluswts <- function(x){
  
  dp <- normij(x)
  n  <- length(x)
  c1 <- sqrt(pi/3)
  c2 <- (n+1)/(n*(n-1))
  
  sig2dpluswts <- (c1*c2*dp)^2
  
  return(sig2dpluswts)
}











