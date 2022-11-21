#' Calculate y'(phi)
#'
#' @description Calculate \eqn{y'(\phi) = y - (||\nu||_2^2)^{-1} \nu \nu^T y + (||\nu||_2^2)^{-1} \nu \phi}.
#'
#' @param y Numeric data vector.
#' @param nu Numeric vector of same length as \code{y}.
#' @param phi Numeric.
#' @param nu2 Value of \eqn{||\nu||_2^2}; optional.
#' @param nuTy Value of \eqn{\nu^T y}; optional.
#'
#' @return Numeric vector.
#' @export
#'
#' @examples
#' y <- rnorm(20) + c(rep(1, 10), rep(-1, 10))
#' nu <- c(rep(0, 5), rep(1/5, 5), rep(-1/5, 5), rep(0, 5))
#' phi <- 2
#' y_phi(y, nu, phi)
#'
y_phi <- function( y, nu, phi, nu2=NULL, nuTy=NULL ){

  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }
  if ( is.null(nuTy) ){
    nuTy <- as.numeric(t(nu) %*% y)
  }
  y2 <- y + nu2^(-1) * (phi - nuTy) * nu

  return(y2)

}
