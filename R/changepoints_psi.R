#' Post-selection inference for changepoints
#'
#' @description
#' Apply a changepoint algorithm (binary segmentation, wild binary segmentation, seeded binary segmentation, or narrowest
#' over threshold) to a vector of data, where the model is piecewise constant mean. Compute the p-value associated with the
#' first changepoint using method of Jewell et al.
#'
#' @param y Numeric vector of data.
#' @param results Output of changepoint algorithm (either \code{binary_segmentation}, \code{wild_binary_segmentation}, or
#' \code{narrowest_over_threshold}).
#' @param nu nu
#' @param threshold Changepoint detection threshold.
#' @param maxiter Integer. Maximum number of changepoints for algorithm to detect.
#' @param eps0 ...
#' @param sigma2 Variance of y.
#' @param h Window size.
#' @param first_cp_only ...
#' @param method Character string; one of \code{"bs"} (binary segmentation), \code{"wbs"} (wild binary segmentation),
#' or \code{"not"} (narrowest over threshold).
#' @param num_rand_ints Number of random intervals for changepoint algorithm (not for \code{method = "bs"}).
#' @param rand_ints Rand
#'
#' @details
#' Either results or threshold and maxiter should be given.
#'
#'
#' @return ...
#'
#'
#' @export
#'
#' @examples
#' # to do
#'
changepoints_psi <- function( y, results=NULL, nu=NULL, threshold=NULL, maxiter=NULL, eps0=0.01, sigma2=NULL,
                                     h=NULL, first_cp_only=FALSE, method="bs", num_rand_ints=NULL, rand_ints=NULL ){

  ### Calculate p-values for binary segmentation, wild/seeded binary segmentation, or narrowest over threshold

  n <- length(y)

  ## Implement changepoint algorithm
  if (!is.null(results)){

    threshold <- results$threshold
    maxiter <- results$maxiter

  } else {

    if (is.null(threshold)){
      threshold <- sd(y) * sqrt( 2 * log( n ) )
    }

    if (is.null(maxiter)){
      maxiter <- n - 1
    }

    if ( method!="bs" & is.null(rand_ints) & is.null(num_rand_ints) ){
      num_rand_ints <- 1000
    }

    results <- find_cps( y, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_ints, rand_ints=rand_ints )

  }

  if ( method!="bs" ){
    rand_ints <- results$rand_ints
  }

  b <- results$results$b[ results$results$cp==1 ]
  d <- results$results$d[ results$results$cp==1 ]
  s <- results$results$s[ results$results$cp==1 ]
  e <- results$results$e[ results$results$cp==1 ]

  if (is.null(nu)){

    if (!is.null(h)){

      if (b[1]-h < 0){
        nu <- c(rep(1/b[1],b[1]), rep(-1/h,h), rep(0,n-b[1]-h))
        nu2 <- 1/b[1] + 1/h
      } else if (b[1]+h > n){
        nu <- c(rep(0,b[1]-h), rep(1/h,h), rep(-1/(n-b[1]),n-b[1]))
        nu2 <- 1/h + 1/(n-b[1])
      } else {
        nu <- c(rep(0,b[1]-h), rep(1/h,h), rep(-1/h,h), rep(0,n-b[1]-h))
        nu2 <- 2/h
      }

    } else {

      cps <- c( 0, sort( b, decreasing=FALSE ), n )
      j <- (1:length(cps))[cps==b[1]]
      nu <- rep( 0, n )
      nu[ (cps[j-1]+1):cps[j] ] <- 1/(cps[j] - cps[j-1])
      nu[ (cps[j]+1):cps[j+1] ] <- -1/(cps[j+1] - cps[j])
      nu2 <- 1/(cps[j] - cps[j-1]) + 1/(cps[j+1] - cps[j])

    }

  } else {

    nu2 <- sum(nu^2)

  }

  nuTy <- as.numeric( t(nu) %*% y )

  ## ?
  if (is.null(sigma2)){
    z <- 0
    b_sorted <- c(0, sort(b), n)
    for ( k in 1:(length(b)+1) ){
      mean_est <- mean(y[(b_sorted[k]+1):b_sorted[k+1]])
      z <- z + sum( (y[(b_sorted[k]+1):b_sorted[k+1]] - mean_est)^2 )
    }
    sigma2 <- z / n
  }





  if ( first_cp_only ){
    if ( method=="bs" ){
      cs <- cusum_phi_vec( y, nu, nu2, nuTy )
      cps2 <- NULL
    } else { #if ( method=="not" ){
      cps2 <- NULL
      cs <- NULL
    }
    interval <- calculate_interval( y, method=method, results=results, nu=nu, cs=cs, nu2=nu2, nuTy=nuTy,
                                    threshold=threshold, n.cp=1 )
    S <- data.frame( matrix( c( interval, as.matrix(c(b[1], d[1])) ), nrow=1 ) )
    colnames(S) <- c( "lower_lim", "upper_lim", "b1", "d1" )

    ncp_max <- 1
    eps <- eps0
    while ( max( S$upper_lim) < Inf ){
      phi <- max(S$upper_lim) + eps
      y2 <- y_phi(y, nu, phi, nu2=nu2, nuTy=nuTy)
      r2 <- find_cps( y2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_ints, rand_ints=rand_ints )
      b2 <- r2$results$b[ r2$results$cp==1 ]
      d2 <- r2$results$d[ r2$results$cp==1 ]
      s2 <- r2$results$s[ r2$results$cp==1 ]
      e2 <- r2$results$e[ r2$results$cp==1 ]

      if ( b[1] %in% b2 ){ ## then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2==b[1]]
      } else {
        n.cp <- maxiter
      }

      interval <- calculate_interval( y, method, nu=nu, results=r2, cs=cs, nu2=nu2, nuTy=nuTy, threshold=threshold, n.cp=n.cp )

      ### Check this interval is the next one
      ncps_found <- min( n.cp, length(b2) )
      if ( abs( interval[1] - max(S$upper_lim) ) < 10^(-10) ){
        if ( ncps_found == ncp_max ){
          S <- rbind( S, c( interval, b2[1:ncps_found], d2[1:ncps_found] ) )
        } else if ( ncps_found <= ncp_max ){
          S <- rbind( S, c( interval, b2[1:ncps_found], rep(NA, ncp_max - ncps_found), d2[1:ncps_found], rep(NA, ncp_max - ncps_found) ) )
        } else {
          S <- cbind( S[,1:(2+ncp_max)], matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max), S[,-(1:(2+ncp_max))],
                 matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max) )
          S <- rbind( S, c( interval, b2[1:ncps_found], d2[1:ncps_found] ) )
          ncp_max <- ncps_found
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }

    }

    eps <- eps0
    while ( min( S$lower_lim ) > -Inf ){
      phi <- min(S$lower_lim) - eps
      y2 <- y_phi(y, nu, phi, nu2=nu2, nuTy=nuTy)
      r2 <- find_cps( y2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_ints, rand_ints=rand_ints )
      b2 <- r2$results$b[ r2$results$cp==1 ]
      d2 <- r2$results$d[ r2$results$cp==1 ]
      s2 <- r2$results$s[ r2$results$cp==1 ]
      e2 <- r2$results$e[ r2$results$cp==1 ]

      if ( b[1] %in% b2 ){ ## then we only have to go as far as b[1] in calculating CPs
        n.cp <- (1:length(b2))[b2==b[1]]
      } else {
        n.cp <- maxiter
      }

      interval <- calculate_interval( y, method, nu=nu, results=r2, cs=cs, nu2=nu2, nuTy=nuTy, threshold=threshold, n.cp=n.cp )

      ### Check this interval is the next one
      ncps_found <- min( n.cp, length(b2) )
      if ( abs( interval[2] - min(S$lower_lim) ) < 10^(-10) ){
        if ( ncps_found == ncp_max ){
          S <- rbind( S, c( interval, b2[1:ncps_found], d2[1:ncps_found] ) )
        } else if ( ncps_found <= ncp_max ){
          S <- rbind( S, c( interval, b2[1:ncps_found], rep(NA, ncp_max - ncps_found), d2[1:ncps_found], rep(NA, ncp_max - ncps_found) ) )
        } else {
          S <- cbind( S[,1:(2+ncp_max)], matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max), S[,-(1:(2+ncp_max))],
                 matrix(NA, nrow=nrow(S), ncol=ncps_found - ncp_max) )
          S <- rbind( S, c( interval, b2[1:ncps_found], d2[1:ncps_found] ) )
          ncp_max <- ncps_found
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }

    }



  } else {



    ### Calculate interval s.t. the given values of b & d are obtained
    if ( method=="bs" ){
      cs <- cusum_phi_vec( y, nu, nu2, nuTy )
      cps2 <- NULL
    } else if ( method=="not" ){
      cps2 <- NULL
      cs <- NULL
    }
    interval <- calculate_interval( y, method, nu=nu, results=results, cs=cs, nu2=nu2, nuTy=nuTy, threshold=threshold, n.cp=maxiter )
    S <- data.frame( matrix( c( interval, as.matrix(c(b[1], d[1])) ), nrow=1 ) )
    colnames(S) <- c( "lower_lim", "upper_lim", "b1", "d1" )

    ncp_max <- length(b)

    ### Find other intervals
    eps <- eps0
    while ( max( S$upper_lim ) < Inf ){
      phi <- max(S$upper_lim) + eps
      y2 <- y_phi(y, nu, phi, nu2=nu2, nuTy=nuTy)
      r2 <- find_cps( y2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_ints, rand_ints=rand_ints )
      b2 <- r2$results$b[ r2$results$cp==1 ]
      d2 <- r2$results$d[ r2$results$cp==1 ]
      s2 <- r2$results$s[ r2$results$cp==1 ]
      e2 <- r2$results$e[ r2$results$cp==1 ]

      interval <- calculate_interval( y, method, nu=nu, results=r2, cs=cs, nu2=nu2, nuTy=nuTy, threshold=threshold, n.cp=maxiter )

      ### Check this interval is the next one
      if ( abs( interval[1] - max(S$upper_lim) ) < 10^(-10) ){
        if ( length(b2)==ncp_max ){
          S <- rbind( S, c( interval, b2, d2 ) )
        } else if ( length(b2)<=ncp_max ){
          S <- rbind( S, c( interval, b2, rep(NA, ncp_max - length(b2)), d2, rep(NA, ncp_max - length(b2)) ) )
        } else {
          S <- cbind( S[,1:(2+ncp_max)], matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max), S[,-(1:(2+ncp_max))],
                 matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max) )
          S <- rbind( S, c( interval, b2, d2 ) )
          ncp_max <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }

    }

    eps <- eps0
    while ( min( S$lower_lim ) > -Inf ){
      phi <- min(S$lower_lim) - eps
      y2 <- y_phi(y, nu, phi, nu2=nu2, nuTy=nuTy)
      r2 <- find_cps( y2, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_ints, rand_ints=rand_ints )
      b2 <- r2$results$b[ r2$results$cp==1 ]
      d2 <- r2$results$d[ r2$results$cp==1 ]
      s2 <- r2$results$s[ r2$results$cp==1 ]
      e2 <- r2$results$e[ r2$results$cp==1 ]

      interval <- calculate_interval( y, method, nu=nu, results=r2, cs=cs, nu2=nu2, nuTy=nuTy, threshold=threshold, n.cp=maxiter )

      ### Check this interval is the next one
      if ( abs( interval[2] - min(S$lower_lim) ) < 10^(-10) ){
        if ( length(b2)==ncp_max ){
          S <- rbind( S, c( interval, b2, d2 ) )
        } else if ( length(b2)<=ncp_max ){
          S <- rbind( S, c( interval, b2, rep(NA, ncp_max - length(b2)), d2, rep(NA, ncp_max - length(b2)) ) )
        } else {
          S <- cbind( S[,1:(2+ncp_max)], matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max), S[,-(1:(2+ncp_max))],
                 matrix(NA, nrow=nrow(S), ncol=length(b2) - ncp_max) )
          S <- rbind( S, c( interval, b2, d2 ) )
          ncp_max <- length(b2)
        }
        eps <- eps0
      } else {
        eps <- eps/2
      }

    }

  }

  ####### Calculate p-values

  ### Take intervals which are for the correct values of b
  if ( first_cp_only ){
    S2 <- S[ S[,3] == b[1], ]
    if ( ncp_max > 1 ){
      for ( j in 4:(2+ncp_max) ){
        S2 <- rbind( S2, S[ S[,j] == b[1], ] )
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
  P_phi_in_S <- sum( pnorm( S2[,2] / sqrt(nu2*sigma2) ) - pnorm( S2[,1] / sqrt(nu2*sigma2) ) )

  ### Calculate P(phi > |nuTy| & phi \in S)
  P_both <- 0
  for ( i in 1:nrow(S2) ){
    if ( abs(nuTy) >= S2[i,1] ){
      if ( abs(nuTy) <= S2[i,2] ){
        P_both <- P_both + pnorm( S2[i,2] / sqrt(nu2*sigma2) ) - pnorm( abs(nuTy) / sqrt(nu2*sigma2) )
      }
    }

    if ( (-1)*abs(nuTy) >= S2[i,1] ){
      if ( (-1)*abs(nuTy) <= S2[i,2] ){
        P_both <- P_both + pnorm( (-1)*abs(nuTy) / sqrt(nu2*sigma2) ) - pnorm( S2[i,1] / sqrt(nu2*sigma2) )
      }
    }
  }

  p_value <- P_both / P_phi_in_S

  ncp_max <- (ncol(S) - 2)/2
  colnames(S) <- colnames(S2) <- c( "lower_lim", "upper_lim", paste0("b", 1:ncp_max), paste0("d", 1:ncp_max) )

  return( list( b=b, d=d, nuTy=nuTy, p_value=p_value, S=S, S2=S2 ) )
}
