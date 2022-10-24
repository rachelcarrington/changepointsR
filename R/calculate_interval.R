#' Calculate interval
#'
#' @description ...
#'
#' @param y vector of data (ignored if cs & nuTy are supplied)
#' @param method Character string; \code{"bs"} for binary segmentation; \code{"wbs"} for wild binary segmentation;
#' \code{"not"} for narrowest over threshold.
#' @param results ...
#' @param nu ...
#' @param nu2 value of ||nu||_2^2
#' @param nuTy value of nu^T y
#' @param threshold minimum changepoint threshold used in binary segmentation algorithm
#' @param n.cp max. number of changepoints to detect of binary segmentation algorithm
#' @param cps2 ?
#' @param cs ?
#'
#' @return A 2-dimensional vector
#' @export
#'
#' @examples
#' x <- 0
#'
calculate_interval <- function(y, method, results, nu, nu2=NULL, nuTy=NULL, threshold=NULL, n.cp=NULL, cps2=NULL, cs=NULL ){

  if ( method=="bs" ){

    b <- results$results$b[ results$results$cp==1 ]
    d <- results$results$d[ results$results$cp==1 ]
    interval <- calculate_interval_bs(y, nu, b, d, cs=cs, nu2=nu2, nuTy=nuTy, threshold=threshold, n.cp=n.cp)

  } else if ( method=="wbs" ){
    
    interval <- calculate_interval_wbs(y, nu, results, nu2=nu2, nuTy=nuTy, threshold=threshold, n.cp=n.cp)
    
  } else if ( method=="not" ){
    
    interval <- calculate_interval_not(y, nu, results, nu2=nu2, nuTy=nuTy, cps2=cps2)
    
  }
  
  return(interval)

}