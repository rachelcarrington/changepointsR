#' Calculate vector of CUSUM statistics for data x, in terms of phi
#'
#' @description Used inside \code{calculate_pvals}.
#' 
#' @param x Numeric vector of data.
#' @param nu Numeric vector.
#' @param nu2 Numeric; value of the squared 2-norm of \code{nu}.
#' @param nuTx Value of nu^T x.
#' @param s Starting point for calculating CUSUM statistics, defaults to 1.
#' @param e Ending point for calculating CUSUM statistics, defaults to \code{length(x)}.
#'
#' @return Returns \code{length(x)} x 2 matrix such that 
#' \code{cusum_phi_vec(x, nu)^T %*% c(1, phi) = cusum(x_phi(x, nu, phi))}.
#'
#' @export
#'
#' @examples
#' x <- rnorm(100) + c(rep(1,40), rep(-1,60))
#' n <- length(x)
#' results <- binary_segmentation(x)
#' b <- results$changepoints[1]
#' nu <- c(rep(1/b, b), rep(-1/(n - b), n - b))
#' cusum_phi_vec(x, nu)
#'
cusum_phi_vec <- function( x, nu, nu2=NULL, nuTy=NULL, s=1, e=length(x), cumsums=NULL, nu_cumsums=NULL ){

  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }

  if ( is.null(nuTx) ){
    nuTx <- as.numeric(t(nu) %*% x)
  }

  cusum_nu <- cusum(nu, s, e, cumsums=nu_cumsums)
  cst <- cusum(x, s, e, cumsums=cumsums) - nuTx/nu2 * cusum_nu
  coef <- cusum_nu / nu2

  return( cbind(cst, coef) )

}
