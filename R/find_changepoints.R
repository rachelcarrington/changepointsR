#' Find changepoints.
#'
#' @description Detect changepoints in the mean using binary segmentation, wild binary segmentation, or narrowest over threshold.
#'
#' @param y Numeric vector of data.
#' @param method Character string: \code{"bs"} for binary segmentation; \code{"wbs"} for wild binary segmentation;
#' \code{"not"} for narrowest over threshold.
#' @param threshold Numeric; changepoint detection threshold for CUSUM statistic.
#' @param maxiter Integer; maximum number of changepoints to detect.
#' @param num_rand_ints Integer; number of random intervals. Ignored if \code{rand_ints} is specified, or if
#' \code{method = "bs"}.
#' @param rand_ints Matrix containing random intervals for wild binary segmentation or narrowest over threshold.
#' Ignored if \code{method = "bs"}. Optional.
#' @param seeded Logical; if \code{TRUE} and \code{method = "wbs"}, then seeded binary segmentation is implemented.
#' @param decay Decay parameter for seeded binary segmentation; only used if \code{method = "wbs"} and \code{seeded = TRUE}.
#'
#' @return A list:
#' \itemize{
#' \item \code{results} Dataframe containing results
#' \item \code{changepoints} Vector of changepoints detected 
#' \item \code{rand_ints} (Not if \code{method = "bs"}) Matrix containing random intervals used
#' \item \code{threshold} Value of \code{threshold}
#' \item \code{maxiter} Value of \code{maxiter}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(100) + c(rep(1,45), rep(-1,10), rep(1,45))
#' results_bs <- find_cps(y, "bs", threshold=4)
#' print(results_bs$results)
#'
#' results_wbs <- find_cps(y, "wbs", threshold=4, num_rand_ints=100)
#' print(results_wbs$results)
#'
#' results_not <- find_cps(y, "not", threshold=4, num_rand_ints=100)
#' print(results_not$results)
#'
find_changepoints <- function(y, method, threshold=NULL, maxiter=NULL, num_rand_ints=NULL, rand_ints=NULL, seeded=FALSE, decay=NULL){

  stopifnot( method == "bs" || method == "wbs" || method == "not" )

  if ( method == "bs" ){

    results <- binary_segmentation(y, threshold, maxiter)

  } else if ( method == "wbs" ){

    if ( !is.null(rand_ints) ){
      num_rand_ints <- nrow(rand_ints)
    }
    results <- wild_binary_segmentation(y, num_rand_samples=num_rand_ints, random_samples=rand_ints, threshold=threshold,
                                          maxiter=maxiter, seeded=seeded, decay=decay)

  } else if ( method == "not" ){

    stopifnot( !is.null(threshold) )
    results <- narrowest_over_threshold(y, threshold, N=num_rand_ints, rand_ints=rand_ints, max_cps=maxiter)

  }

  return(results)

}
