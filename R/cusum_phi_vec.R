#' Calculate vector of CUSUM statistics for data y, in terms of phi
#'
#' @description Used inside binary_segmentation_psi (etc.)
#'
#' @param y Vector of data
#' @param nu nu
#' @param nu2 Value of ||nu||_2^2
#' @param nuTy Value of nu^T y
#' @param s Starting point for calculating CUSUM statistics, defaults to 1
#' @param e Ending point for calculating CUSUM statistics, degaults to n
#'
#' @return Returns n x 2 matrix s.t. cusum_phi_vec(y)^T (1,phi) = cusm(y'(phi))
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

  if ( is.null( nu2 ) ){
    nu2 <- sum( nu^2 )
  }
  if ( is.null( nuTy ) ){
    nuTy <- as.numeric( t(nu) %*% y )
  }
  cusum_nu <- cusum(nu, s, e)
  cst <- cusum(y, s, e) - nuTy/nu2 * cusum_nu
  coef <- 1/nu2 * cusum_nu

  return( cbind( cst, coef ) )

}
