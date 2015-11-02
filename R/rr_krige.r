#' Make spatial predictions
#' 
#' @description
#' Obtains spatial predictions at selected locations
#' 
#' @rdname rr_predict
#' 
#' @param method the prediction method to use. Currently supported are "universal" and "lognormal".
#' @param ... space for additional arguments (passed to relelvent prediction function).
#' 
#' @return
#' A matrix with the coordinates of the predicted locations, the predicted value, and measures of uncertainty.
#' 
#' @export
#' 

rr_predict <- function( method="universal", ... ){
  
  
  if( method=="universal" ){
    kriging <- rr_universal_krige( ... )
  }  
  
  
  if( method=="lognormal" ){
    kriging <- rr_lognormal_krige( ... )
  }  
  
  
  
  
  return(kriging)
  
  
}



