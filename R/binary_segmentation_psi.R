#' Post-selection inference for binary segmentation
#'
#' @description Calculate p-values associated with detected changepoints for binary segmentation;
#' uses method of Jewell et al.
#'
#' @details
#' Given a changepoint of interest \eqn{\tau_j}, there are two options for the null hypothesis:
#' \itemize{
#' \item There are no changepoints within a window size \eqn{h} of \eqn{\tau_j}. If \code{nus = NULL} but \code{h} is supplied,
#' then the value of \eqn{\nu} for each changepoint will be calculated using this assumption.
#' \item There are no other changepoints between \eqn{\tau_{j-1}} and \eqn{\tau_{j+1}}. If \eqn{\tau_j} is the first
#' changepoint, then \eqn{\tau_{j-1}} is taken to be 0; if it is the last changepoint, \eqn{\tau_{j+1}} is taken to be 
#' \code{length(y)}. If \code{nus = NULL} and \code{h = NULL}, then the value of \eqn{\nu} for each changepoint will be calculated 
#' using this assumption.
#' }
#' Alternatively \code{nus} can be specified manually.
#'
#' @param y Numeric vector of data.
#' @param results Output of binary_segmentation function; can be NULL.
#' @param nus List of nu's for each detected changepoint. See details.
#' @param threshold Numeric; minimum threshold for CUSUM statistic for changepoint detection. Ignored if results specified.
#' @param maxiter Number of changepoints to find. Ignored if results specified; otherwise defaults to \code{length(y) - 1}.
#' @param eps0 Hyperparameter for calculating S. (Details.)
#' @param sigma2 Variance of \code{y}. If unknown, it can be estimated within the function.
#' @param h Window size. See details.
#' @param first_cp_only Logical. If \code{TRUE}, condition on the fact that the changepoint of interest is in the model; 
#' if \code{FALSE}, condition on all changepoints. Defaults to \code{TRUE} if \code{h} is supplied, and \code{FALSE} otherwise.
#' @param num_pvals Integer. Maximum number of p-values to calculate; defaults to \code{length(y) - 1}. If the number of changepoints
#' detected by the binary segmentation algorithm is less than or equal to \code{num_pvals}, then p-values will be calculated for all 
#' changepoints. If \code{num_pvals} is less than the number of changepoints detected, p-values will only be calculated for the first
#' \code{num_pvals} changepoints, in order of detection.
#'
#' @return A list containing:
#' \itemize{
#' \item \code{b} Numeric vector of changepoints
#' \item \code{p_value} Numeric vector of p-values
#' }
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(200) + c(rep(1,50), rep(0,100), rep(1,50))
#' results <- binary_segmentation(y)
#' pvals <- binary_segmentation_psi(y, results, sigma2=1)
#' print(pvals)
#'
binary_segmentation_psi <- function( y, results=NULL, nus=NULL, threshold=NULL, maxiter=NULL, eps0=0.01, sigma2=NULL,
                                       h=NULL, first_cp_only=FALSE, num_pvals=NULL ){

  n <- length(y)

  ## Implement binary segmentation
  if ( !is.null(results) ){
    threshold <- results$threshold
    maxiter <- results$maxiter
  } else {
    if ( is.null(threshold) ){
      threshold <- sd(y) * sqrt(2 * log(n))
    }
    if ( is.null(maxiter) ){
      maxiter <- n - 1
    }
    results <- binary_segmentation(y, threshold=threshold, maxiter=maxiter)
  }

  b <- results$results$b[ results$results$cp==1 ]
  d <- results$results$d[ results$results$cp==1 ]

  if ( length(b)==0 ){
    stop("No changepoints detected.")
  } else if ( is.null(num_pvals)  || length(b) < num_pvals ){
    num_pvals <- length(b)
  }

  if ( is.null(nus) ){
    nus <- as.list( rep(NA, num_pvals) )
  }

  if ( is.null(sigma2) ){
    sigma2 <- 0
    b_sorted <- c(0, sort(b), n)
    for ( k in 1:(length(b) + 1) ){
      mean_est <- mean(y[(b_sorted[k] + 1):b_sorted[k+1]])
      sigma2 <- sigma2 + sum( (y[(b_sorted[k] + 1):b_sorted[k+1]] - mean_est)^2 )
    }
    sigma2 <- sigma2 / (n - 1)
  }

  if ( !is.null(h) & !(first_cp_only) ){
    first_cp_only <- TRUE
  }

  p_value <- rep(NA, num_pvals)

  for ( jj in 1:num_pvals ){

    if ( is.na(nus[[jj]]) ){
      if ( !is.null(h) ){
        if ( b[jj] - h < 0 ){
          nu <- c( rep(1/b[jj], b[jj]), rep(-1/h, h), rep(0, n - b[jj] - h) )
          nu2 <- 1/b[jj] + 1/h
        } else if ( b[jj] + h > n ){
          nu <- c( rep(0, b[jj] - h), rep(1/h, h), rep(-1/(n - b[jj]),n - b[jj]) )
          nu2 <- 1/h + 1/(n - b[jj])
        } else {
          nu <- c( rep(0, b[jj] - h), rep(1/h, h), rep(-1/h, h), rep(0, n - b[jj] - h) )
          nu2 <- 2/h
        }
      } else {
        cps <- c( 0, sort(b, decreasing=FALSE), n )
        j <- (1:length(cps))[cps==b[jj]]
        nu <- rep(0, n)
        nu[ (cps[j-1] + 1):cps[j] ] <- 1/(cps[j] - cps[j-1])
        nu[ (cps[j] + 1):cps[j+1] ] <- -1/(cps[j+1] - cps[j])
        nu2 <- 1/(cps[j] - cps[j-1]) + 1/(cps[j+1] - cps[j])
      }
    } else {
      nu <- nus[[jj]]
      nu2 <- sum(nu^2)
    }

    nuTy <- as.numeric(t(nu) %*% y)

    if ( jj == 1 ){
      S <- calculate_S(y, results=results, nu=nu, threshold=threshold, maxiter=maxiter, method="bs",
                        nuTy=nuTy, first_cp_only=first_cp_only)
    } else {
      S <- calculate_S(y, results=results, nu=nu, threshold=threshold, maxiter=maxiter, method="bs",
                        nuTy=nuTy)
    }


    ####### Calculate p-values

    ### Take intervals which are for the correct values of b
    ncp_max <- (ncol(S) - 2)/2
    if ( first_cp_only ){
      S2 <- S[ S[,3] == b[jj], ]
      if ( ncp_max > 1 ){
        for ( j in 4:(2 + ncp_max) ){
          S2 <- rbind( S2, S[ S[,j] == b[jj], ] )
        }
      }
      S2 <- S2[ !is.na( S2[,3] ), ]
    } else {
      if ( ncp_max > length(b) ){
        S2 <- S[ is.na( S[,2+(length(b)+1)] ), ] ## remove rows which contain too many changepoints
      } else {
        S2 <- S
      }
      for ( j in 1:length(b) ){
        S2 <- S2[ S2[,j+2] %in% b, ] ## check each found CP is in b
      }
    }


    ### Calculate P(phi \in S)
    P_phi_in_S <- sum( pnorm(S2[,2] / sqrt(nu2*sigma2)) - pnorm(S2[,1] / sqrt(nu2*sigma2)) )

    ### Calculate P(phi > |nuTy| & phi \in S)
    P_both <- 0

    for ( i in 1:nrow(S2) ){
      if ( S2[i,2] > abs(nuTy) ){
        P_both <- P_both + pnorm(S2[i,2] / sqrt(nu2*sigma2)) - pnorm(max(c(abs(nuTy), S2[i,1])) / sqrt(nu2*sigma2))
      }
      if ( S2[i,1] < (-1)*abs(nuTy) ){
        P_both <- P_both + pnorm(min(c(-abs(nuTy), S2[i,2])) / sqrt(nu2*sigma2)) - pnorm(S2[i,1] / sqrt(nu2*sigma2))
      }

    }

    p_value[jj] <- P_both / P_phi_in_S

    ncp_max <- (ncol(S) - 2)/2
    colnames(S) <- colnames(S2) <- c( "lower_lim", "upper_lim", paste0("b", 1:ncp_max), paste0("d", 1:ncp_max) )

  }

  return( list( b=b, p_value=p_value ) )
}
