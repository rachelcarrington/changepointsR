#' Post-selection inference for narrowest over threshold change in mean model.
#'
#' @description Calculates p-value for first detected changepoint, when the narrowest over threshold algorithm
#' is applied to the change in mean model.
#'
#' @param y Numeric vector of data.
#' @param lambda Threshold for determining changepoint candidates.
#' @param results Output of narrowest over threshold algorithm (optional).
#' @param max_cps Maximum number of changepoints for the algorithm to detect; defaults to 1; ignored if results
#' specified.
#' @param N Number of random intervals; defaults to 1000; ignored if results or rand_ints specified.
#' @param rand_ints Random intervals for narrowest over threshold algorithm; ignored if results specified.
#'
#' @return List of results ...
#' @export
#'
#' @examples
#' # ...
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
    nuTy <- as.numeric( t(nu) %*% y )

    cps2 <- cps[(j-1):(j+1)]
      ## will be useful later; nu is constant on any intervals which don't contain one of these points

    rand_ints <- results$rand_ints

    ### Calculate interval s.t. the given values of b & d are obtained
    interval <- calculate_interval_not(y, nu, results, cps2=cps2)
      ## need to rewrite this to use vectorised cs
    S <- data.frame( matrix( c(interval, as.matrix(c(b, d, s, e))), nrow=1 ) )
    colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:length(b)), paste0("d", 1:length(d)),
                     paste0("s", 1:length(s)), paste0("e", 1:length(e)))

    ### Find other intervals
    #### !!!! Need to deal with the problem that the # of CPs may vary between rows of S
    eps <- 0.01
    while ( max(S$upper_lim) < Inf ){
      phi <- max(S$upper_lim) + eps
      y2 <- y_phi(y, nu, phi, nu2=nu2, nuTy=nuTy)
      r2 <- narrowest_over_threshold(y2, lambda, max_cps=max_cps, rand_ints=results$rand_ints)
      if ( nrow(r2$results) >= 1 ){
        b2 <- r2$results$b
        d2 <- r2$results$d
        s2 <- r2$results$s
        e2 <- r2$results$e
      } else {
        b2 <- d2 <- s2 <- e2 <- NA
      }

      interval <- calculate_interval_not_multiple_cps_2(y, nu, r2, cps2)

      ### Check this interval is the next one
      if ( abs(interval[1] - max(S$upper_lim)) < 10^(-10) ){
        S <- rbind( S, c(interval, b2, d2, s2, e2) )
        eps <- 0.01
      } else {
        eps <- eps/2
      }

    }

    eps <- 0.01
    while ( min(S$lower_lim) > -Inf ){
      phi <- min(S$lower_lim) - eps
      y2 <- y_phi(y, nu, phi, nu2=nu2, nuTy=nuTy)
      r2 <- narrowest_over_threshold(y2, lambda, max_cps=max_cps, rand_ints=results$rand_ints)
      if ( nrow(r2$results) >= 1 ){
        b2 <- r2$results$b
        d2 <- r2$results$d
        s2 <- r2$results$s
        e2 <- r2$results$e
      } else {
        b2 <- d2 <- s2 <- e2 <- NA
      }

      interval <- calculate_interval_not(y, nu, r2, cps2)

      ### Check this interval is the next one
      if ( abs(interval[2] - min(S$lower_lim)) < 10^(-10) ){
        S <- rbind( S, c(interval, b2, d2, s2, e2) )
        eps <- 0.01
      } else {
        eps <- eps/2
      }

    }

    ### Calculate probability

    ### Take intervals which are for the correct values of b
    S2 <- S[ S[,3]==0, ] ## empty data.frame
    S2 <- S2[ !is.na( S2[,1] ), ]
    b_sorted <- sort(b)
    for ( i in 1:nrow(S) ){
      if ( !is.na(sum(S[i,3:(2+length(b))])) & sum(abs(sort( S[i,3:(2+length(b))]) - b_sorted)) < 10^(-10)  ){
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

    return( list(b=b, d=d, p_value=p_value, S=S) )

  } else {

    print( "No changepoints detected by narrowest-over-threshold algorithm." )
    return( list(b=NULL, d=NULL, s=NULL, e=NULL, p_value=NULL, S=NULL) )
  }

}
