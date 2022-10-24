#' Calculate y'(phi)
#'
#' @description ...
#'
#' @param y Numeric data vector
#' @param nu Numeric data vector of same length as y
#' @param phi Scalar
#' @param nu2 Scalar; squared 2-norm of nu
#' @param nuTy Scalar; value of nu^T y
#'
#' @return Numeric vector
#' @export
#'
#' @examples
#' y <- rnorm(20) + c(rep(1,10), rep(-1,10))
#' nu <- c(rep(0,5), rep(1/5,5), rep(-1/5,5), rep(0,5))
#' phi <- 2
#' y_phi(y, nu, phi)
#'
y_phi <- function( y, nu, phi, nu2=NULL, nuTy=NULL ){

  ## Calculate y'(phi), given y, nu and phi

  if ( is.null(nu2) ){
    nu2 <- sum( nu^2 )
  }
  if ( is.null(nuTy) ){
    nuTy <- as.numeric( t(nu) %*% y )
  }
  y2 <- y + nu2^(-1) * (phi - nuTy) * nu

  return(y2)

}
