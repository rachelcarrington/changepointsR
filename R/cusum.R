#' Calculate CUSUM statistics for a vector.
#'
#' @param y A numeric vector.
#' @param s Starting index.
#' @param e Ending index.
#' @param return_full Logical. If \code{TRUE}, returns vector of same length as \code{y}, with indices outside of \code{(s, e - 1)} set to
#' \code{NA}. If \code{FALSE}, a vector of length \code{e - s + 1} is returned.
#' @param cumsums Vector of cumulative sums of either \code{y} or \code{y[s:e]}. If not supplied, this will be calculated within the function.
#'
#' @description
#' Calculate CUSUM statistics for a vector \code{y}.
#'
#' \code{cumsums} can be supplied to reduce computation time if \code{cusum} is applied to (subsets of) the same dataset multiple times. Otherwise
#' it can be left as \code{NULL} and calculated within the function.
#'
#' @return A numeric vector.
#'
#' @export
#'
#' @examples
#' y <- c(1, 3, 2, 5, 3)
#' cusum(y)
#'
cusum <- function( y, s=1, e=length(y), return_full=FALSE, cumsums=NULL ){

  stopifnot( s >= 1, s < length(y), e > 1, e <= length(y) )

  if ( is.null(cumsums) ){
    cumsums <- cumsum(y[s:e])
  } else if ( length(cumsums) == length(y) ){
    if ( s > 1 ){
      cumsums <- cumsums[s:e] - cumsums[s - 1]
    } else {
      cumsums <- cumsums[s:e]
    }
  } else { # (incorrectly specified cumsums)
    cumsums <- cumsum(y[s:e])
  }

  n <- e - s + 1
  z <- 1:(n - 1)
  a <- sqrt(z * rev(z) / n)

  cusum_vec <- a * (cumsums[-n] / z - (cumsums[n] - cumsums[-n]) / rev(z))
  if ( return_full ){
    cusum_vec <- c(rep(NA, s - 1), cusum_vec, rep(NA, length(y) - e))
  }

  return(cusum_vec)

}
