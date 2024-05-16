#' Calculate S for autocorrelated data.
#'
#' @param x Numeric vector of data.
#' @param nu Numeric vector.
#' @param Sigma 
#' @param results Output of changepoint algorithm (either \code{binary_segmentation}, \code{wild_binary_segmentation}, or 
#' \code{narrowest_over_threshold}.
#' @param b Vector of changepoints; ignored if \code{results} specified.
#' @param d Vector containing directions of changepoints (+1 or -1); ignored if \code{results} specified.
#' @param threshold Threshold for detecting changepoints; ignored if \code{results} specified.
#' @param maxiter Maximum number of changepoints to detect; ignored if \code{results} specified.
#' @param nuTx Value of \eqn{nu^T x}; optional.
#' @param eps0 Hyperparameter.
#' @param first_cp_only Logical. If \code{TRUE}, condition on the fact that the changepoint of interest is in the model; 
#' if \code{FALSE}, condition on all changepoints. Defaults to \code{FALSE}.
#' @param method Character. One of \code{"bs"} (binary segmentation), \code{"wbs"} (wild binary segmentation), or \code{"not"} (narrowest
#' over threshold). Defaults to \code{"bs"}.
#' @param rand_ints Matrix containing random intervals for changepoint algorithm. Ignored if \code{results} specified
#' or \code{method = "bs"}.
#' @param seeded Logical. For \code{method = "wbs"} only, whether to use seeded binary segmentation.
#' @param decay Decay parameter for seeded binary segmentation. Only used if \code{method = "wbs"} and \code{seeded = TRUE}.
#'
#' @return A dataframe containing intervals with the changepoints obtained when \eqn{\phi} is in each interval.
#' @export
#'
#' @examples x
calculate_S_autocor <- function( x, nu, Sigma, results=NULL, b=NULL, d=NULL, threshold=NULL, maxiter=NULL, nuTx=NULL, eps0=0.01, first_cp_only=FALSE,
                                 method="bs", rand_ints=NULL, seeded=FALSE, decay=NULL ){
  
  ## Calculate S given x, b, d
  
  if ( is.null(results) ){
    if ( is.null(b) || is.null(d) || is.null(maxiter) || is.null(threshold) ){
      stop("If results is not specified, then b, d, maxiter, and threshold must all be given.")
    }
  }
  
  n <- length(x)
  
  if ( is.null(nuTx) ){
    nuTx <- c(t(nu) %*% x)
  }
  
  if ( is.null(b) ){
    b <- results$results$b[results$results$cp == 1]
  }
  if ( is.null(d) ){
    d <- results$results$d[results$results$cp == 1]
  }
  if ( is.null(maxiter) ){
    maxiter <- results$maxiter
  }
  if ( is.null(threshold) ){
    threshold <- results$threshold
  }
  
  if ( method != "bs" ){
    if ( is.null(rand_ints) ){
      rand_ints <- results$rand_ints
    }
    num_rand_ints <- nrow(rand_ints)
  }
  
  if ( first_cp_only ){
    interval <- calculate_interval_autocor(x, nu, b, d, threshold=threshold, Sigma=Sigma, nuTx=nuTx, n.cp=1)
    S <- matrix(c(interval, as.matrix(c(b[1], d[1]))), nrow=1)
    colnames(S) <- c("lower_lim", "upper_lim", "b1", "d1")
    
    ncp_max <- 1
    eps <- eps0
    while ( max(S[,"upper_lim"]) < Inf ){
      phi <- max(S[,"upper_lim"]) + eps
      x2 <- x_phi_autocor(x, nu, phi, Sigma=Sigma, nuTx=nuTx)
      r2 <- find_changepoints(x2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_ints, rand_ints=rand_ints,
                              seeded=seeded, decay=decay)
      b2 <- r2$results$b[r2$results$cp == 1]
      d2 <- r2$results$d[r2$results$cp == 1]
      
      if ( b[1] %in% b2 ){ # then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2 == b[1]]
      } else {
        n.cp <- maxiter
      }
      
      interval <- calculate_interval_autocor(x, nu, b2, d2, threshold=threshold, Sigma=Sigma, nuTx=nuTx, n.cp=n.cp)
      
      # Check this interval is the next one
      ncps_found <- min(n.cp, length(b2))
      if ( abs(interval[1] - max(S[,"upper_lim"])) < 10^(-10) ){
        if ( ncps_found == ncp_max ){
          S <- rbind(S, c(interval, b2[1:ncps_found], d2[1:ncps_found]))
        } else if ( ncps_found <= ncp_max ){
          S <- rbind(S, c(interval, b2[1:ncps_found], rep(NA, ncp_max - ncps_found), d2[1:ncps_found], rep(NA, ncp_max - ncps_found)))
        } else {
          S <- cbind( S[,1:(2 + ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max) )
          S <- rbind(S, c(interval, b2[1:ncps_found], d2[1:ncps_found]))
          ncp_max <- ncps_found
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
    eps <- eps0
    while ( min(S[,"lower_lim"]) > -Inf ){
      phi <- min(S[,"lower_lim"]) - eps
      x2 <- x_phi_autocor(x, nu, phi, Sigma=Sigma, nuTx=nuTx)
      r2 <- find_changepoints(x2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_samples, rand_ints=rand_ints,
                              seeded=seeded, decay=decay)
      b2 <- r2$results$b[r2$results$cp == 1]
      d2 <- r2$results$d[r2$results$cp == 1]
      
      if ( b[1] %in% b2 ){ # then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2 == b[1]]
      } else {
        n.cp <- maxiter
      }
      
      interval <- calculate_interval_autocor(x, nu, b2, d2, threshold=threshold, Sigma=Sigma, nuTx=nuTx, n.cp=n.cp)
      
      # Check this interval is the next one
      ncps_found <- min(n.cp, sum(!is.na(b2)))
      if ( abs(interval[2] - min(S[,"lower_lim"])) < 10^(-10) ){
        if ( ncps_found == ncp_max ){
          S <- rbind( S, c(interval, b2[1:ncps_found], d2[1:ncps_found]) )
        } else if ( ncps_found <= ncp_max ){
          S <- rbind( S, c(interval, b2[1:ncps_found], rep(NA, ncp_max - ncps_found), d2[1:ncps_found], rep(NA, ncp_max - ncps_found)) )
        } else {
          S <- cbind( S[,1:(2 + ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max) )
          S <- rbind(S, c(interval, b2[1:ncps_found], d2[1:ncps_found]))
          ncp_max <- ncps_found
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
    
  } else {
    
    # Calculate interval s.t. the given values of b & d are obtained
    interval <- calculate_interval_autocor(x, nu, b, d, threshold=threshold, Sigma=Sigma, nuTx=nuTx, n.cp=maxiter)
    if ( length(b) >= 1 ){
      S <- matrix(c(interval, as.matrix(c(b, d))), nrow=1)
      colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:length(b)), paste0("d", 1:length(d)))
    } else {
      S <- matrix(c(interval, as.matrix(c(NA, NA))), nrow=1)
      colnames(S) <- c("lower_lim", "upper_lim", "b1", "d1")
    }
    
    ncp_max <- max(length(b), 1)
    
    # Find other intervals
    eps <- eps0
    while ( max(S[,"upper_lim"]) < Inf ){
      phi <- max(S[,"upper_lim"]) + eps
      x2 <- x_phi_autocor(x, nu, phi, Sigma=Sigma, nuTx=nuTx)
      r2 <- find_changepoints(x2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_samples, rand_ints=rand_ints,
                              seeded=seeded, decay=decay)
      b2 <- r2$results$b[r2$results$cp == 1]
      d2 <- r2$results$d[r2$results$cp == 1]
      
      interval <- calculate_interval_autocor(x, nu, b2, d2, threshold=threshold, Sigma=Sigma, nuTx=nuTx, n.cp=maxiter)
      
      # Check this interval is the next one
      if ( abs(interval[1] - max(S[,"upper_lim"])) < 10^(-10) ){
        if ( length(b2) == ncp_max ){
          S <- rbind(S, c(interval, b2, d2))
        } else if ( length(b2) <= ncp_max ){
          S <- rbind(S, c(interval, b2, rep(NA, ncp_max - length(b2)), d2, rep(NA, ncp_max - length(b2))))
        } else {
          S <- cbind( S[,1:(2 + ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max) )
          S <- rbind(S, c(interval, b2, d2))
          ncp_max <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
    eps <- eps0
    while ( min(S[,"lower_lim"]) > -Inf ){
      phi <- min(S[,"lower_lim"]) - eps
      x2 <- x_phi_autocor(x, nu, phi, Sigma=Sigma, nuTx=nuTx)
      r2 <- find_changepoints(x2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_samples, rand_ints=rand_ints,
                              seeded=seeded, decay=decay)
      b2 <- r2$results$b[r2$results$cp == 1]
      d2 <- r2$results$d[r2$results$cp == 1]
      
      interval <- calculate_interval_autocor(x, nu, b2, d2, threshold=threshold, Sigma=Sigma, nuTx=nuTx, n.cp=maxiter)
      
      # Check this interval is the next one
      if ( abs( interval[2] - min(S[,"lower_lim"]) ) < 10^(-10) ){
        if ( length(b2) == ncp_max ){
          S <- rbind(S, c(interval, b2, d2))
        } else if ( length(b2) <= ncp_max ){
          S <- rbind(S, c(interval, b2, rep(NA, ncp_max - length(b2)), d2, rep(NA, ncp_max - length(b2))))
        } else {
          S <- cbind( S[,1:(2+ncp_max),drop=FALSE], matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max), S[,-(1:(2 + ncp_max)),drop=FALSE],
                      matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max) )
          S <- rbind(S, c(interval, b2, d2))
          ncp_max <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }
      
    }
    
  }
  
  S <- data.frame(S)
  colnames(S) <- c("lower_lim", "upper_lim", paste0("b", 1:ncp_max), paste0("d", 1:ncp_max))
  
  return(S)
}