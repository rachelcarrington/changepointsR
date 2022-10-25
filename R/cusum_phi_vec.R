#' Calculate vector of CUSUM statistics for data y, in terms of phi
#'
#' @description Used inside \code{binary_segmentation_psi} (etc.)
#' 
#' @param y Numeric vector of data.
#' @param nu Numeric vector.
#' @param nu2 Numeric; value of the squared 2-norm of \code{nu}.
#' @param nuTy Value of nu^T y.
#' @param s Starting point for calculating CUSUM statistics, defaults to 1.
#' @param e Ending point for calculating CUSUM statistics, degaults to \code{length(y)}.
#'
#' @return Returns \code{length(y)} x 2 matrix such that 
#' \code{cusum_phi_vec(y, nu)^T %*% c(1, phi) = cusum(y_phi(y, nu, phi))}.
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(100) + c(rep(1,40), rep(-1,60))
#' results <- binary_segmentation(y)
#' b <- results$changepoints[1]
#' nu <- c(rep(1/b, b), rep(-1/(100-b), 100-b))
#' cusum_phi_vec(y, nu)
#'
cusum_phi_vec <- function( y, nu, nu2=NULL, nuTy=NULL, s=1, e=length(y) ){

  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }

  if ( is.null(nuTy) ){
    nuTy <- as.numeric(t(nu) %*% y)
  }

  cusum_nu <- cusum(nu, s, e)
  cst <- cusum(y, s, e) - nuTy/nu2 * cusum_nu
  coef <- 1/nu2 * cusum_nu

  return( cbind(cst, coef) )

}
