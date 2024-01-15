#' Wild Binary Segmentation
#'
#' @description Implements wild binary segmentation algorithm for changepoint detection.
#'
#' @param y Numeric vector of data.
#' @param num_rand_samples Number of random intervals to use. Ignored if \code{random_samples} is supplied.
#' @param random_samples N x 2 matrix of random intervals.
#' @param threshold Minimum threshold for determining changepoint candidates.
#' @param maxiter Maximum number of changepoints to find.
#' @param seeded Logical; if \code{TRUE}, seeded binary segmentation is used.
#' @param decay Decay parameter for seeded binary segmentation. Only used if \code{seeded = TRUE}.
#'
#' @return A list.
#' \itemize{
#' \item \code{results} Dataframe containing results
#' \item \code{changepoints} Vector of changepoints detected
#' \item \code{rand_ints} N x 2 matrix containing random intervals used in the algorithm.
#' \item \code{threshold} Value of \code{threshold}
#' \item \code{maxiter} Value of \code{maxiter}
#' }
#'
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

  n <- length(y)

  if ( is.null(threshold) ){
    threshold <- sd(y) * sqrt(2 * log(n))
  }

  if ( is.null(maxiter) ){
    maxiter <- n
  }

  # Draw random/seeded samples
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

  # Calculate changepoint estimate for each interval
  results0 <- numeric(0)
  for ( ind in 1:num_rand_samples ){
    cusum_stats <- cusum(y, random_samples[ind,1], random_samples[ind,2], return_full=TRUE)
    b_hat <- which.max(abs(cusum_stats))
    d_hat <- ifelse(cusum_stats[b_hat] > 0, -1, 1)
    cs <- cusum_stats[b_hat]
    cp <- ifelse(abs(cs) > threshold, 1, 0)
    results0 <- rbind(results0, c(ind, random_samples[ind,1], random_samples[ind,2], b_hat, d_hat, cs, cp))
  }

  results0 <- data.frame(results0)
  colnames(results0) <- c("index", "s", "e", "b", "d", "cs", "cp")

  # Keep only intervals where the CUSUM statistic is above the threshold
  results0 <- results0[results0$cp == 1, ]
  if( nrow(results0) == 0 ){
    stop("No changepoints detected.")
  }

  # Order by absolute value of CUSUM statistics
  results0 <- results0[order(abs(results0$cs), decreasing=TRUE), ]

  results <- numeric(0)
  iter <- 1
  while( iter <= maxiter ){
    results <- rbind(results, results0[1,])
    b_hat <- results0$b[1]
    cp <- results0$cp[1]

    # Remove intervals which contain the changepoint
    contains_changepoint <- (b_hat >= results0$s) & (b_hat < results0$e)
    results0 <- results0[!(contains_changepoint), ]

    # Check whether there are any intervals left; if not, end loop
    iter <- ifelse(nrow(results0) == 0, maxiter + 1, iter + 1)
  }

  return( list(results=results, changepoints=results$b, rand_ints=random_samples, threshold=threshold, maxiter=maxiter) )
}
