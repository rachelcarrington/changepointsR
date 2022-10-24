#' Find changepoints.
#'
#' @description
#' Implements binary segmentation, wild binary segmentation, narrowest over threshold for change in mean model.
#'
#' @param y Numeric vector of data.
#' @param method Character string; \code{"bs"} for binary segmentation; \code{"wbs"} for wild binary segmentation;
#' \code{"not"} for narrowest over threshold.
#' @param threshold Numeric; changepoint detection threshold for CUSUM statistic.
#' @param maxiter Integer; maximum number of changepoints to detect.
#' @param num_rand_ints Integer; number of random intervals. Ignored if \code{rand_ints} is specified, or if
#' \code{method = "bs"}.
#' @param rand_ints Matrix containing random intervals for wild binary segmentation or narrowest over threshold.
#' Ignored if \code{method = "bs"}.
#' @param seeded Logical; if \code{TRUE} and \code{method = "wbs"}, then seeded binary segmentation is implemented.
#' @param decay Decay parameter for seeded binary segmentation; only used if \code{method = "wbs"} and \code{seeded = TRUE}.
#'
#' @return ...
#' @export
#'
#' @examples
#' # ...
#'
find_cps <- function(y, method, threshold=NULL, maxiter=NULL, num_rand_ints=NULL, rand_ints=NULL, seeded=FALSE, decay=NULL){

  if ( method=="bs" ){

    results <- binary_segmentation(y, threshold, maxiter)

  } else if ( method=="wbs" ){

    if ( !is.null(rand_ints) ){
      num_rand_ints <- nrow(rand_ints)
    }
    results <- wild_binary_segmentation(y, num_rand_samples=num_rand_ints, random_samples=rand_ints, threshold=threshold,
                                          maxiter=maxiter, seeded=seeded, decay=decay)

  } else if ( method=="not" ){

    results <- narrowest_over_threshold(y, threshold, N=num_rand_ints, rand_ints=rand_ints, max_cps=maxiter)

  }

  return(results)

}
