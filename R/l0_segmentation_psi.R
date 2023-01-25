#' Post-selection inference for L0 segmentation
#'
#' @description Post-selection inference for L0 segmentation. This uses the functions \code{changepoint_estimates} and
#' \code{changepoint_inference} from the package \code{ChangepointInference}.
#'
#' @param y Numeric vector of data
#' @param lambda Threshold parameter
#' @param N Number of samples of psi to take
#' @param h Window size
#' @param sigma2 Variance of \code{y}
#' @param sig Tuning parameter
#' @param include_original Logical; whether to include observed value as psi in place as one of the random samples;
#' defaults to \code{TRUE}
#' @param num_pvals Maximum number of p-values to calculate
#'
#' @return A list:
#' \itemize{
#' \item \code{b} Vector of changepoints
#' \item \code{p_value} Vector of p-values
#' \item \code{p_value_orig} Vector of p-values obtained using fixed \eqn{\psi = \psi_{obs}}
#' \item \code{nuTy} Value of \eqn{\nu^T y}
#' \item \code{P_both} Matrix containing values of \eqn{Pr(|\phi| > |\phi_{obs}| \& \phi \in S)}
#' \item \code{P_phi_in_S} Matrix containing values of \eqn{Pr(\phi \in S)}
#' \item \code{P_both_orig} Vector containing values of \eqn{Pr(|\phi| > |\phi_{obs}| \& \phi \in S)} for fixed \eqn{\psi = \psi_{obs}}
#' \item \code{P_phi_in_S_orig} Vector containing values of \eqn{Pr(\phi \in S)} for fixed \eqn{\psi = \psi_{obs}}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(100) + c(rep(1,40), rep(-1,20), rep(1,40))
#' l0_segmentation_psi(y, 4, 10, 10)
#'
l0_segmentation_psi <- function(y, lambda, N, h, sigma2=1, sig=4, include_original=TRUE, num_pvals=NULL){

#  library(ChangepointInference)

  n <- length(y)

  fit <- ChangepointInference::changepoint_estimates(y, "L0", lambda)
  b <- fit$change_pts
  if ( is.null( num_pvals ) ){
    num_pvals <- length(fit$change_pts)
  }

  if ( include_original ){
    ### Original p-values
    P_both_orig <- P_phi_in_S_orig <- rep(0, num_pvals)
    pvals_orig <- rep(NA, num_pvals)
    z <- ChangepointInference::changepoint_inference(y, "L0-fixed", tuning_parameter=lambda, window_size=h, sig=sig, return_conditioning_set=TRUE)

    for ( jj in 1:num_pvals ){
      nu <- c(rep(0,b[jj]-h), rep(1/h,h), rep(-1/h,h), rep(0,n-b[jj]-h))
      nu2 <- 2/h
      nuTy <- as.numeric(t(nu)%*%y)

      S <- z$conditioning_sets[[jj]]
      S2 <- S[ S$contained==1, ]

      ### Calculate P(phi \in S)
      P_phi_in_S_orig[jj] <- sum( pnorm( S2$max_mean / sqrt(nu2*sigma2) ) - pnorm( S2$min_mean / sqrt(nu2*sigma2) ) )

      ### Calculate P(phi > |nuTy| & phi \in S)
      for ( i in 1:nrow(S2) ){
        if ( S2$max_mean[i] > abs(nuTy) ){
          P_both_orig[jj] <- P_both_orig[jj] +
            pnorm( S2$max_mean[i] / sqrt(nu2*sigma2) ) - pnorm( max( c( abs(nuTy), S2$min_mean[i] ) ) / sqrt(nu2*sigma2) )
        }
        if ( S2$min_mean[i] < (-1)*abs(nuTy) ){
          P_both_orig[jj] <- P_both_orig[jj] +
            pnorm( min( c( -abs(nuTy), S2$max_mean[i] ) ) / sqrt(nu2*sigma2) ) - pnorm( S2$min_mean[i] / sqrt(nu2*sigma2) )
        }
      }
      pvals_orig[jj] <- P_both_orig[jj] / P_phi_in_S_orig[jj]
    }

    N2 <- N - 1

  } else {

    pvals_orig <- P_both_orig <- P_phi_in_S_orig <- NA
    N2 <- N

  }


  ### New p-values

  if ( length(b) >= 1 & N2 >= 1 ){

    alpha <- rep(1, 2*h)
    alpha2 <- 2*h
    Z <- diag(2*h) - 1/alpha2 * crossprod(t(alpha)) - 1/(2/h) * crossprod(t(c(rep(1/h,h), rep(-1/h,h))))
    U <- svd(Z)$u[,1:(2*h-2)]

    if ( is.null(num_pvals) ){
      num_pvals <- length(b)
    } else {
      num_pvals <- min( c(num_pvals, length(b)) )
    }

    p_val <- rep(NA, num_pvals)
    P_both <- P_phi_in_S <- matrix(NA, nrow=num_pvals, ncol=N2)

    for ( j in 1:num_pvals ){

      if ( b[j] >= h & b[j] <= (n-h-1) ){

        nu <- c(rep(0,b[j]-h), rep(1/h,h), rep(-1/h,h), rep(0,n-b[j]-h))
        nu2 <- 2/h
        nuTy <- as.numeric(t(nu)%*%y)

        Psi <- matrix( rnorm( N2*(2*h-2), sd=sqrt(sigma2) ), nrow=N2 )
        for ( iter in 1:N2 ){
          y_new <- y
          y_new[(b[j]-h+1):(b[j]+h)] <- y_new[(b[j]-h+1):(b[j]+h)] + U %*% Psi[iter,]

          fit_new <- ChangepointInference::changepoint_estimates(y_new, "L0", lambda)
          b_new <- fit_new$change_pts

          ### I think this works (we need to make sure that b[j] %in% b_new)
          phi <- 10
          while ( !( b[j] %in% b_new ) ){
            y_new <- y_phi(y_new, nu, phi)
            phi <- phi + 10
            b_new <- ChangepointInference::changepoint_estimates(y_new, "L0", lambda)$change_pts
          }

          z <- ChangepointInference::changepoint_inference(y_new, 'L0-fixed', lambda, window_size = h, sig = sig, return_conditioning_sets = TRUE)

          ### Make sure we select the right changepoint!
          S <- z$conditioning_sets[[ (1:length(z$change_pts))[ z$change_pts==b[j] ] ]]

          ### Keep rows of S for which b[j] %in% b ("contained" tells us whether b_new[1] %in% b)
          S2 <- S[ S$contained==1, ]

          ### Calculate P(phi \in S)
          P_phi_in_S[j,iter] <- sum( pnorm( S2$max_mean / sqrt(nu2*sigma2) ) - pnorm( S2$min_mean / sqrt(nu2*sigma2) ) )

          ### Calculate P(phi > |nuTy| & phi \in S)
          P_both[j,iter] <- 0
          for ( i in 1:nrow(S2) ){
            if ( S2$max_mean[i] > abs(nuTy) ){
              P_both[j,iter] <- P_both[j,iter] +
                pnorm( S2$max_mean[i] / sqrt(nu2*sigma2) ) - pnorm( max( c( abs(nuTy), S2$min_mean[i] ) ) / sqrt(nu2*sigma2) )
            }
            if ( S2$min_mean[i] < (-1)*abs(nuTy) ){
              P_both[j,iter] <- P_both[j,iter] +
                pnorm( min( c( -abs(nuTy), S2$max_mean[i] ) ) / sqrt(nu2*sigma2) ) - pnorm( S2$min_mean[i] / sqrt(nu2*sigma2) )
            }

          }

        }

        ### Calculate p-value estimates
        if ( include_original ){
          p_val[j] <- ( sum(P_both[j,]) + P_both_orig[j] ) / ( sum(P_phi_in_S[j,]) + P_phi_in_S_orig[j] )
        } else {
          p_val[j] <- sum( P_both[j,] ) / sum( P_phi_in_S[j,] )
        }

      }

    }


  } else {

    p_val <- numeric(0)
    print("No changepoints detected.")

  }

  return(list(b=b, p_value=p_val, p_value_orig=pvals_orig, nuTy=nuTy, P_both=P_both, P_phi_in_S=P_phi_in_S, P_both_orig=P_both_orig,
                P_phi_in_S_orig=P_phi_in_S_orig))
}

