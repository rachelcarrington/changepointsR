#' Calculate CUSUM statistics for a vector
#'
#' @param y A numeric vector.
#' @param s Starting index.
#' @param e Ending index.
#' @param return_full Logical. If \code{TRUE}, returns vector of same length as \code{y}, with indices outside of \code{(s, e - 1)} set to
#' \code{NA}. If \code{FALSE}, a vector of length \code{e - s + 1} is returned.
#'
#' @description
#' Calculate CUSUM statistics for a vector \code{y}.
#'
#' @return A numeric vector.
#'
#' @export
#'
#' @examples
#' y <- c(1, 3, 2, 5, 3)
#' cusum(y)
#'
cusum <- function( y, s=1, e=length(y), return_full=FALSE ){

  stopifnot( s >= 1, s < length(y), e > 1, e <= length(y) )

  ## Calculates CUSUM statistics for a vector y
  ## (s,e) : start and end indices, to calculate CUSUM statistics for just a subset of y
  ## note that CUSUM statistics are calculated for y_s, ..., y_{e-1}
  ## return_full : if TRUE, returns vector of same length as y, with values outside (s,e-1) replaced by NA
    ## else, returns just vector of CUSUM statistics for (s,e-1)

  n <- e - s + 1
  x <- 1:(n - 1)
  a <- sqrt(x * rev(x) / n)

  cusum_vec <- a * (cumsum(y[s:(e - 1)]) / x - rev(cumsum(y[e:(s + 1)]) / x))
  if ( return_full ){
    cusum_vec <- c(rep(NA, s - 1), cusum_vec, rep(NA, length(y) - e))
  }

  return(cusum_vec)

}
