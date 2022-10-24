#' Binary segmentation
#'
#' @description
#' Binary segmentation algorithm for detecting changepoints, in the change in mean model.
#' The algorithm terminates when either (i) the maximum number of changepoints have been found; or
#' (ii) there are no points remaining where the CUSUM statistic is above the threshold.
#'
#' @param y A numeric vector.
#' @param threshold Threshold for determining changepoint candidates; defaults to sqrt(2*log(n)*var(y)).
#' @param maxiter Maximum number of changepoints to find; defaults to length(y) - 1
#'
#' @return A list.
#' @export
#'
#' @examples
#' set.seed(10)
#' y <- rnorm(100) + c(rep(1,20), rep(-1,30), rep(1,50))
#' binary_segmentation(y)
#'
binary_segmentation <- function( y, threshold=NULL, maxiter=NULL ){

  n <- length(y)

  ## Set threshold, if unspecified
  if ( is.null(threshold) ){
    threshold <- sd(y) * sqrt(2 * log(n))
  }

  ## Set maximum number of iterations (default to n-1, if unspecified)
  if ( is.null(maxiter) ){
    maxiter <- n - 1
  } else if ( maxiter > n-1 ){
    maxiter <- n - 1
  }

  ## Find first changepoint
  iter <- 1
  s <- 1
  e <- n
  cusum_stats <- cusum(y)
  b <- which.max( abs(cusum_stats) ) ## changepoint candidate
  lrs <- max( abs(cusum_stats) ) ## CUSUM statistic at this point
  d <- ifelse( cusum_stats[b] > 0, -1, 1 ) ## direction of change
  cp <- ( lrs > threshold ) ## = 1 if CUSUM statistic is above threshold, 0 otherwise
  results <- data.frame( iter, s, e, b, d, lrs, cp )

  if ( cp & maxiter > 1 ){

    ## Split dataset at changepoints and repeat
    while ( iter < maxiter ) {

      if ( results$cp[length(results$cp)]==1 ){

        ### If a CP was detected last time
        iter <- iter + 1

        ### Take results for which the changepoint threshold was reached
        cps <- sort( c(results$b[ results$cp==1 ]) )
        s <- c(1, cps + 1)
        e <- c(cps, n)

        cusum_stats <- rep(0, n - 1)
        for ( k in 1:length(s) ){
          if ( e[k] > s[k] ){
            cusum_stats[ s[k]:(e[k] - 1) ] <- cusum(y[s[k]:e[k]])
          }
        }

        b <- which.max( abs(cusum_stats) ) ## changepoint candidate
        lrs <- max( abs(cusum_stats) ) ## CUSUM statistic at this point
        d <- ifelse( cusum_stats[b] > 0, -1, 1 ) ## direction of change
        cp <- ( lrs > threshold ) ## = 1 if CUSUM statistic is above threshold, 0 otherwise
        results <- rbind( results, c(iter, max(s[s <= b]), min(e[e > b]), b, d, lrs, cp) )

      } else {

        ### If no CP was detected, end while loop
        iter <- maxiter

      }

    }

  }

  changepoints <- results$b[ results$cp==1 ]

  return( list( results=results, changepoints=changepoints, threshold=threshold, maxiter=maxiter ) )

}
