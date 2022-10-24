#' Narrowest over threshold for change in mean.
#'
#' @param y A numeric vector of data.
#' @param lambda Threshold; a scalar.
#' @param rand_ints Matrix containing random intervals; optional.
#' @param N Number of random intervals; defaults to 1000.
#' @param max_cps Maximum number of changepoints to return.
#'
#' @return ...
#' @export
#'
#' @examples
#' y <- rnorm(200) + c(rep(1,90), rep(-1,20), rep(1,90))
#' lambda <- 4
#' results <- narrowest_over_threshold(y, lambda, N=500)
#' print(results$results)
#'
narrowest_over_threshold <- function( y, lambda, rand_ints=NULL, N=1000, max_cps=NULL ){

  ## y : vector of data
  ## lambda : threshold
  ## rand_ints : N x 2 matrix of random intervals
  ## N : number of random intervals; ignored if rand_ints is not NULL
  ## max_cps : max. number of changepoints to find, defaults to length(y) - 1

  n <- length(y)
  if ( is.null( max_cps ) ){
    max_cps <- n - 1
  }

  if ( is.null( rand_ints ) ){
    rand_ints <- generate_random_intervals( n, N, 2 )
  }

  ### Sort intervals by width
  interval_widths <- rand_ints[,2] - rand_ints[,1] + 1
  rand_ints <- rand_ints[ order( interval_widths ), ]
  interval_widths <- sort( interval_widths )

  ### For each interval, calculate CUSUM statistic
  changepoint_candidates <- data.frame( index=0, s=0, e=0, b=0, d=0, lrs=0, width=0, cp=0 )
  changepoint_candidates <- changepoint_candidates[ changepoint_candidates$index > 0, ]
  for ( interval_index in 1:nrow(rand_ints) ){
    s <- rand_ints[interval_index,1]
    e <- rand_ints[interval_index,2]
    cusum_stats <- rep( 0, n-1 )
    cusum_stats[s:(e-1)] <- cusum( y[s:e] )
    if ( max(abs(cusum_stats)) > lambda ){
      b1 <- which.max( abs( cusum_stats ) )
      cusum_max <- cusum_stats[b1]
      d1 <- ifelse( cusum_max < 0, 1, -1 )
      changepoint_candidates <- rbind( changepoint_candidates, c( interval_index, s, e, b1, d1, cusum_max, (e-s+1), 1 ) )
    }
  }
  colnames(changepoint_candidates) <- c("index", "s", "e", "b", "d", "lrs", "width", "cp")

  ### Within each width, sort by CUSUM statistic
  widths <- unique(changepoint_candidates$width)
  for ( width in widths ){
    z2 <- changepoint_candidates[ changepoint_candidates$width==width, ]
    if ( nrow(z2) > 1 ){
      changepoint_candidates[ changepoint_candidates$width==width, ] <- z2[ order( abs( z2$lrs ), decreasing=TRUE ), ]
    }
  }

  ### Create empty dataframe
  results <- changepoint_candidates[ changepoint_candidates$index<0, ]

  ### Get results
  while( nrow(changepoint_candidates) > 0 ){
    results <- rbind( results, changepoint_candidates[1,] )
    changepoint_candidates <- changepoint_candidates[ !( changepoint_candidates$s <= changepoint_candidates$b[1] & changepoint_candidates$e > changepoint_candidates$b[1] ), ]
    if ( nrow(results) == max_cps ){
      break
    }
  }

  if ( nrow(results) > max_cps ){
    results <- results[1:max_cps,]
  }

  return( list( results=results, changepoints=results$results$b, rand_ints=rand_ints, threshold=lambda, maxiter=max_cps ) )
}
