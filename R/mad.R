#' Estimate variance by median absolute deviation
#'
#' @description Estimate variance of a vector using median absolute deviation, under the assumption that the data follows a
#' Normal distribution.
#'
#' @param y Numeric vector of data.
#'
#' @return Scalar: estimated variance of \code{y}.
#' @export
#'
#' @examples
#' y <- rnorm(200, sd=2) + c(rep(1, 100), rep(-1, 100))
#' mad(y)
#'
mad <- function(y){

  difs <- y[-1] - y[-length(y)]
  mad_est <- median(abs(difs - median(difs))) / qnorm(0.75)
  mad_est <- mad_est^2 / 2

  return(mad_est)
}
