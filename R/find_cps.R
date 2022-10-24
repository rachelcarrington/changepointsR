## dependent functions: binary_segmentation, wild_binary_segmentation, narrowest_over_threshold, cusum, generate_random_intervals

#' Title
#'
#' @description
#' Implements binary segmentation, wild binary segmentation, narrowest over threshold for change in mean model.
#'
#' @param y vector of data
#' @param method "bs" = binary segmenation; "wbs" = wild binary segmentation; "not" = narrowest over threshold
#' @param threshold changepoint detection threshold for CUSUM statistic
#' @param maxiter max. number of changepoints to detect
#' @param num_rand_ints number of random intervals (ignored if rand_ints specified, or if method = "bs")
#' @param rand_ints random intervals for wild binary segmentation or narrowest over threshold (ignored if method = "bs")
#' @param seeded if TRUE and method = "wbs", then seeded binary segmentation is implemented
#' @param decay decay parameter for seeded binary segmentation; only used if method = "wbs" and seeded = TRUE
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
