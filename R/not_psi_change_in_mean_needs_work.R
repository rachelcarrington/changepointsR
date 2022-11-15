#' Post-selection inference for narrowest over threshold change in mean model.
#'
#' @description For the change in mean model, using narrowest over threshold algorithm to detect changepoints.
#' Calculate p-value for first detected changepoint.
#'
#' @param y Numeric vector of data.
#' @param lambda Numeric; threshold for determining changepoint candidates.
#' @param results Output from applying \code{narrowest_over_threshold} to \code{y} (optional).
#' @param max_cps Maximum number of changepoints for the algorithm to detect; defaults to 1. Ignored if \code{results}
#' specified.
#' @param N Number of random intervals; defaults to 1000. Ignored if \code{results} or \code{rand_ints} specified.
#' @param rand_ints Matrix of random intervals for narrowest over threshold algorithm. Ignored if \code{results} specified.
#'
#' @return List:
#' \itemize{
#' \item \code{b} Vector of changepoints.
#' \item \code{d} Vector containing directions of change for each changepoint (\code{1} for a positive change, \code{-1} for a 
#' negative change.
#' \item \code{p_value} P-value.
#' \item \code{S} Dataframe containing intervals for \eqn{\phi} and changepoints obtained when \eqn{\phi} is in each interval.
#' }
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(100) + c(rep(1,20), rep(-1,20), rep(1,20), rep(-1,20), rep(1,20))
#' results <- narrowest_over_threshold(y, 4, N=50)
#' not_psi_change_in_mean(y, 4, results)
#'
not_psi_change_in_mean <- function( y, lambda, results=NULL, max_cps=1, N=1000, rand_ints=NULL ){

  ### Post-selection inference following Jewell et al., for narrowest over threshold algorithm
  ### Tests significance of first CP only.

  n <- length(y)

  ### Apply NOT algorithm to find changepoints
  if ( is.null(results) ){
    results <- narrowest_over_threshold(y, lambda, max_cps=max_cps, N=N, rand_ints=rand_ints)
  }

  b <- results$results$b
  d <- results$results$d
  s <- results$results$s
  e <- results$results$e

  if ( length(b) > 0 ){
    cps <- c(0, sort(b, decreasing=FALSE), n)
    j <- (1:length(cps))[cps==b[1]]
    nu <- rep(0, n)
    nu[ (cps[j-1]+1):cps[j] ] <- 1/(cps[j] - cps[j-1])
    nu[ (cps[j]+1):cps[j+1] ] <- -1/(cps[j+1] - cps[j])

    nu2 <- 1/(cps[j] - cps[j-1]) + 1/(cps[j+1] - cps[j])
    nuTy <- as.numeric(t(nu) %*% y)

    rand_ints <- results$rand_ints

    S <- calculate_S(y, nu, results, method="not")

    ### Calculate probability

    ### Take intervals which are for the correct values of b
    S2 <- S[ S[,3]==0, ] ## empty data.frame
    S2 <- S2[ !is.na( S2[,1] ), ]
    b_sorted <- sort(b)
    for ( i in 1:nrow(S) ){
      if ( !is.na(sum(S[i,3:(2+length(b))])) & sum(abs(sort(S[i,3:(2+length(b))]) - b_sorted)) < 10^(-10)  ){
        S2 <- rbind(S2, S[i,])
      }
    }

    ### Calculate P(phi \in S)
    P_phi_in_S <- sum( pnorm(S2[,2] / sqrt(nu2)) - pnorm(S2[,1] / sqrt(nu2)) )

    ### Calculate P(phi > |nuTy| & phi \in S)
    P_both <- 0
    for ( i in 1:nrow(S2) ){
      if ( abs(nuTy) >= S2[i,1] & abs(nuTy) <= S2[i,2] ){
        P_both <- P_both + pnorm(S2[i,2] / sqrt(nu2)) - pnorm(abs(nuTy) / sqrt(nu2))
      }

    if ( (-1)*abs(nuTy) >= S2[i,1] & (-1)*abs(nuTy) <= S2[i,2] ){
      P_both <- P_both + pnorm((-1)*abs(nuTy) / sqrt(nu2)) - pnorm(S2[i,1] / sqrt(nu2))
      }
    }

    p_value <- P_both / P_phi_in_S

    S <- S[ order(S[,1]), ]

    return( list(b=b, d=d, p_value=p_value, S=S) )

  } else {

    print( "No changepoints detected by narrowest-over-threshold algorithm." )
    return( list(b=NULL, d=NULL, s=NULL, e=NULL, p_value=NULL, S=NULL) )
  }

}
