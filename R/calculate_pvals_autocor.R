#' Calculate p-values for change in mean model, with autocorrelation
#'
#' @param x Numeric vector of data.
#' @param results Output of changepoint algorithm (either \code{binary_segmentation}, \code{wild_binary_segmentation}, or
#' \code{narrowest_over_threshold}).
#' @param Sigma Covariance matrix of \code{x}; not needed if \code{rho} and \code{sigma2} are supplied.
#' @param rho Autocorrelation parameter for AR(1) model; ignored if \code{Sigma} is supplied.
#' @param sigma2 Noise variance for AR(1) model; ignored if \code{Sigma} is supplied.
#' @param h Positive integer (>=2) or vector of length 2 or \code{NULL}. Window size. If \code{NULL} then the changepoints either side of 
#' the changepoint of interest are used to define the window.
#' @param eps0 Hyperparameter for calculating S; defaults to 0.01. (This will only affect the time taken to run the algorithm, not the results.)
#' @param num_pvals Number of p-values to compute. If \code{NULL}, the function will calculate p-values for all detected changepoints.
#' 
#' @details 
#' The model is \eqn{X_t = \mu_t + Z_t}, where \eqn{\mu_t} is constant except at changepoints and \eqn{Z_t \sim N(0, \Sigma)}.
#' Either the covariance matrix \eqn{\Sigma} can be supplied to the function or, if \eqn{Z_t} follows an AR(1) model of the form
#' \eqn{Z_t = \rho Z_{t-1} + \sigma \epsilon_t} with eqn{\epsilon_t \sim N(0, 1)}, then \eqn{\rho} and \eqn{\sigma^2} can be supplied instead.
#'
#' @return x
#' @export
#'
#' @examples x
#' 
calculate_pvals_autocor <- function( x, results, Sigma=NULL, rho=NULL, sigma2=1, h=NULL, eps0=0.01, num_pvals=NULL ){
  
  n <- length(x)
  
  threshold <- results$threshold
  maxiter <- results$maxiter
  b <- results$results$b[results$results$cp == 1]
  d <- results$results$d[results$results$cp == 1]
  
  if ( length(b) == 0 ){
    stop("No changepoints detected")
  } else if ( is.null(num_pvals) ){
    num_pvals <- length(b)
  } else if ( length(b) < num_pvals ){
    print(paste0("Only ", length(b), " changepoints were detected."))
    num_pvals <- length(b)
  }
  
  p_value <- rep(NA, num_pvals)
  P_both <- P_phi_in_S <- rep(0, num_pvals)
  
  if ( is.null(Sigma) ){
    Sigma <- diag((1 - rho^seq(2, 2*n, by=2))/(1 - rho^2))
    for ( i in 1:(n - 1) ){
      Sigma[(i + 1):n, i] <- Sigma[i, (i + 1):n] <- rho^seq(1:(n - i)) * Sigma[i, i]
    }
    Sigma <- sigma2 * Sigma
  }

  for ( jj in 1:num_pvals ){
    
    # Calculate nu
    if ( is.null(h) ){ # use changepoint on either side to determine h
      
      cps <- c(0, sort(b, decreasing=FALSE), n)
      j <- (1:length(cps))[cps==b[jj]]
      h1 <- ceiling((cps[j] - cps[j-1]))
      h2 <- ceiling((cps[j+1] - cps[j]))
      
    } else if ( length(h) >= 2 ) {
      
      cps <- c(0, sort(b[-jj]), n)
      if ( sum(cps %in% (b[jj] - h[1] + 1):b[jj]) >= 1 ){
        h1 <- b[jj] - max(cps[cps %in% (b[jj] - h[1] + 1):b[jj]])
      } else {
        h1 <- h[1]
      }
      if ( sum(cps %in% (b[jj] + 1):(b[jj] + h[2])) >= 1 ){
        h2 <- max(cps[cps %in% (b[jj] + 1):(b[jj] + h[2])]) - b[jj]
      } else {
        h2 <- h[2]
      }
      
    } else {
      
      stopifnot( h >= 2 )
      cps <- c(0, sort(b[-jj]), n)
      if ( sum(cps %in% (b[jj] - h + 1):b[jj]) >= 1 ){
        h1 <- b[jj] - max(cps[cps %in% (b[jj] - h + 1):b[jj]])
      } else {
        h1 <- h
      }
      if ( sum(cps %in% (b[jj] + 1):(b[jj] + h)) >= 1 ){
        h2 <- max(cps[cps %in% (b[jj] + 1):(b[jj] + h)]) - b[jj]
      } else {
        h2 <- h
      }
      
    }
    
    nu <- c(rep(0, b[jj] - h1), rep(1/h1, h1), rep(-1/h2, h2), rep(0, n - b[jj] - h2))
    nh <- h1 + h2    
    nuTx <- as.numeric(t(nu) %*% x)

    # Calculate S
    if ( length(b) >= 1 ){
      if ( !is.null(h) & jj == 1 ){
        S_all <- calculate_S_autocor(x, nu=nu, Sigma=Sigma, results=results, threshold=threshold, maxiter=maxiter, method="bs", first_cp_only=TRUE, nuTx=nuTx)
      } else {
        S_all <- calculate_S_autocor(x, results=results, Sigma=Sigma, nu=nu, threshold=threshold, maxiter=maxiter, method="bs", nuTx=nuTx)
      }
    } else {
      S_all <- calculate_S_autocor(x, results=results, Sigma=Sigma, nu=nu, threshold=threshold, maxiter=maxiter, method="bs", nuTx=nuTx)
    }
    
    
    # Calculate S
    if ( is.null(h) ){
      
      # Find which intervals contain the correct combination of changepoints
      ## First remove intervals which contain too many changepoints
      max_cps_found <- (ncol(S_all) - 2)/2
      if ( max_cps_found > length(results$changepoints) ){
        S <- S_all[is.na(S_all[, (2 + length(results$changepoints) + 1)]), ]
      } else {
        S <- S_all
      }
      ## Remove intervals which contain too few changepoints
      S <- S[!is.na(S[, (2 + length(results$changepoints))]), ]
      ## Check we have the correct combination of changepoints
      z <- rep(0, nrow(S))
      cps <- sort(results$changepoints)
      for ( j in 1:nrow(S) ){
        if ( sum(abs(sort(S[j, 3:(length(cps) + 2)]) - cps)) < 10^(-10) ){
          z[j] <- 1
        }
      }
      S <- S[z == 1,]
      
    } else {
      
      # Find intervals which contain b[jj]
      max_cps_found <- (ncol(S_all) - 2)/2
      S <- S_all[S_all[,3] == b[jj], ]
      if ( max_cps_found > 1 ){
        for ( j in 2:max_cps_found ){
          S <- rbind(S, S_all[S_all[,j + 2] == b[jj], ])
        }
      }
      S <- S[!is.na(S$b1), ]
    }
    
    # Calculate variance of nuTx
    nuTx_var <- sqrt(c(t(nu) %*% Sigma %*% nu))
    
    # Calculate P(phi \in S)
    P_phi_in_S[jj] <- sum(pnorm(S[,2] / nuTx_var) - pnorm(S[,1] / nuTx_var))
    
    # Calculate P(phi > |nuTx| & phi \in S)
    P_both[jj] <- 0
    if ( nrow(S2) >= 1 ){
      for ( i in 1:nrow(S) ){
        if ( S$upper_lim[i] > abs(nuTx) ){
          P_both[jj] <- P_both[jj] + pnorm(S$upper_lim[i] / nuTx_var) - pnorm(max( c(abs(nuTx), S$lower_lim[i]) ) / nuTx_var)
        }
        if ( S$lower_lim[i] < (-1)*abs(nuTx) ){
          P_both[jj] <- P_both[jj] + pnorm(min(c(-abs(nuTx), S$upper_lim[i])) / nuTx_var) - pnorm(S$lower_lim[i] / nuTx_var)
        }
      }
      
    }
    
    # Calculate p-value
    p_value[jj] <- P_both[jj] / P_phi_in_S[jj]
    
  }
  
  return( list(p_value=p_value, P_both=P_both, P_phi_in_S=P_phi_in_S) )
  
}