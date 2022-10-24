#' Wild Binary Segmentation
#'
#' @description Implements wild binary segmentation algorithm for changepoint detection.
#'
#' @param y A vector of data
#' @param num_rand_samples Number of random intervals to use (ignored if random_samples supplied)
#' @param random_samples N x 2 matrix of random intervals
#' @param threshold Minimum threshold for determining changepoint candidates
#' @param maxiter Maximum number of changepoints to find
#' @param seeded Logical; if TRUE, seeded binary segmentation is used
#' @param decay Decay parameter for seeded binary segmentation
#'
#' @return A list
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(100) + c(rep(1,40), rep(0,10), rep(1,50))
#' results <- wild_binary_segmentation(y, num_rand_samples=500)
#' print(results$results)
#'
wild_binary_segmentation <- function( y, num_rand_samples=1000, random_samples=NULL, threshold=NULL, maxiter=NULL, seeded=FALSE,
                                        decay=sqrt(2) ){

  ## Wild / Seeded Binary Segmentation
  #### random_samples : N x 2 matrix of random intervals. If this is supplied, these are the intervals used
  #### seeded : whether to implement seeded binary segmentation. If TRUE (and random_samples not supplied), seeded intervals are
  ###### selected (num_rand_samples is ignored)
  #### num_rand_samples : number of random samples to draw (if seeded is FALSE and random_samples not supplied; otherwise it is ignored)
  #### threshold : threshold for detecting changepoints, a default value is used if unspecified
  #### maxiter : maximum number of changepoints to detect, defaults to number of data points if not specified

  n <- length(y)

  if ( is.null(threshold) ){
    threshold <- sd(y) * sqrt(2 * log(n))
  }

  if ( is.null(maxiter) ){
    maxiter <- n
  }

  ### Draw random/seeded samples
  if ( is.null(random_samples) ){
    if ( seeded ){
      random_samples <- seeded.intervals(n, decay)
      num_rand_samples <- nrow(random_samples)
    } else {
      random_samples <- generate_random_intervals(n, num_rand_samples)
    }
  } else {
    num_rand_samples <- nrow(random_samples)
  }

  ### Calculate estimate for each interval
  results0 <- numeric(0)
  for ( ind in 1:num_rand_samples ){
    cusum_stats <- cusum(y, random_samples[ind,1], random_samples[ind,2], return_full=TRUE)
    b_hat <- which.max( abs(cusum_stats) )
    d_hat <- ifelse( cusum_stats[b_hat] > 0, -1, 1 )
    cs <- cusum_stats[b_hat]
    cp <- ifelse( abs(cs) > threshold, 1, 0 )
    results0 <- rbind( results0, c(ind, random_samples[ind,1], random_samples[ind,2], b_hat, d_hat, cs, cp) )
  }

  results0 <- data.frame(results0)
  colnames(results0) <- c("index", "s", "e", "b", "d", "cs", "cp")
  results0 <- results0[results0$cp==1,] ## delete intervals below threshold
  results0 <- results0[order(abs(results0$cs), decreasing=TRUE),]

  results <- numeric(0)
  iter <- 1
  while( iter <= maxiter ){
    results <- rbind(results, results0[1,])
    b_hat <- results0$b[1]
    cp <- results0$cp[1]

    ### Remove intervals which contain the changepoint
    z <- (b_hat >= results0$s) & (b_hat < results0$e)
    results0 <- results0[!(z),]

    ### Check whether there are any intervals left
    if ( nrow(results0)==0 ){
      iter <- maxiter + 1
    } else {
      iter <- iter + 1
    }
  }

  return( list( results=results, threshold=threshold, maxiter=maxiter, rand_ints=random_samples ) )
}
