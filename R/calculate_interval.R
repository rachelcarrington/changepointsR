#' Calculate interval
#'
#' @description Wrapper function for \code{calculate_interval_bs}, \code{calculate_interval_wbs}, and \code{calculate_interval_not}.
#'
#' @param y Numeric vector of data
#' @param method Character string; \code{"bs"} for binary segmentation; \code{"wbs"} for wild binary segmentation;
#' \code{"not"} for narrowest over threshold.
#' @param results Output of changepoint algorithm (\code{binary_segmentation}, \code{wild_binary_segmentation}, or \code{narrowest_over_threshold}).
#' @param nu Numeric vector.
#' @param nu2 Value of \eqn{||\nu||_2^2}.
#' @param nuTy Value of \eqn{\nu^T y}.
#' @param n.cp Maximum number of changepoints to detect.
#'
#' @return A 2-dimensional vector.
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(100) + c(rep(1,50), rep(-1,50))
#' results <- binary_segmentation(y, threshold=4)
#' b <- results$results$b[ results$results$cp==1 ]
#' h <- 10
#' nu <- c(rep(0, b[1]-h), rep(1/h, h), rep(-1/h, h), rep(0, length(y)-b[1]-h))
#' calculate_interval(y, "bs", results, nu)
#'
calculate_interval <- function(y, method, results, nu, nu2=NULL, nuTy=NULL, n.cp=NULL ){

  if ( method=="bs" ){

    b <- results$results$b[ results$results$cp==1 ]
    d <- results$results$d[ results$results$cp==1 ]
    interval <- calculate_interval_bs(y, nu, b, d, threshold=results$threshold, nu2=nu2, nuTy=nuTy, n.cp=n.cp)

  } else if ( method=="wbs" ){
    
    interval <- calculate_interval_wbs(y, nu, results, nu2=nu2, nuTy=nuTy, n.cp=n.cp)
    
  } else if ( method=="not" ){
    
    interval <- calculate_interval_not(y, nu, results, nu2=nu2, nuTy=nuTy)
    
  }
  
  return(interval)

}