#' Calculate p-values for change in mean model
#'
#' @description Calculate p-values for detected changepoints, for the change in mean model. Changepoints can be detected using
#' either binary segmentation, wild binary segmentation, or narrowest over threshold.
#'
#' @param y Numeric vector of data.
#' @param method Method used to generate changepoints. Options: \code{"bs"} for binary segmentation,
#' \code{"wbs"} for wild or seeded binary segmentation, \code{"not"} for narrowest over threshold.
#' @param results Output of changepoint algorithm (either \code{binary_segmentation}, \code{wild_binary_segmentation}, or
#' \code{narrowest_over_threshold}).
#' @param N Positive integer. Number of random \eqn{\psi}'s to generate; defaults to 10.
#' @param threshold Minimum threshold for changepoint detection. Ignored if \code{results} is specified.
#' @param maxiter Positive integer. Maximum number of changepoints to find in changepoint algorithm. Ignored if
#' \code{results} is specified
#' @param h Positive integer (>=2) or vector of length 2 or \code{NULL}. Window size. If \code{NULL} then the changepoints either side of 
#' the changepoint of interest are used to define the window.
#' @param gamma Number in \code{(0, 1]}. If \code{h = NULL}, then \eqn{h} is taken to be \code{gamma} times the distance between the changepoint
#' of interest and the estimated changepoints on either side. Ignored if \code{h} specified; defaults to \code{1}.
#' @param cp_bound Logical. If \code{TRUE}, then if there is an estimated changepoint within the window that is fixed under the null hypothesis,
#' this changepoint will be used as the boundary of the window.
#' @param sigma2 Variance of \code{y}. Defaults to \code{1} if not specified.
#' @param eps0 Hyperparameter for calculating S.
#' @param include_original Logical. Whether to include the \eqn{\psi} value corresponding to the observed data in place of
#' one of the random samples. Defaults to \code{TRUE}.
#' @param num_pvals Integer or NA. Maximum number of p-values to calculate; if set to NA, it will default to \code{length(y) - 1}.
#' If the number of changepoints detected by the binary segmentation algorithm is less than or equal to \code{num_pvals}, then
#' p-values will be calculated for all changepoints. If \code{num_pvals} is less than the number of changepoints detected,
#' p-values will only be calculated for the first \code{num_pvals} changepoints, in order of detection.
#' @param random_samples Matrix containing random intervals. For wild binary segmentation and narrowest over
#' threshold only. Ignored if \code{results} specified or \code{method = "bs"}.
#' @param num_rand_samples Number of random intervals to use for wild binary segmentation or narrowest over threshold.
#' Ignored if \code{random_samples} or \code{results} is specified, \code{method = "bs"}, or \code{seeded = TRUE}.
#' @param seeded Logical. If \code{TRUE}, use seeded binary segmentation rather than wild binary segmentation. Only used if
#' \code{method = "wbs"}.
#' @param decay Decay parameter for seeded binary segmentation. Only used if \code{method = "wbs"} and \code{seeded = TRUE}.
#' @param return_probs Logical. If \code{TRUE}, the values of \eqn{Pr(\phi \in S \& |\phi| > |\phi_{obs}|)} and \eqn{Pr(\phi \in S)}
#' for each \eqn{\psi} will be included in the output.
#'
#' @details
#' Given a changepoint of interest \eqn{\tau_j}, there are two options for the null hypothesis:
#' \itemize{
#' \item There are no changepoints within a window size \eqn{h} of \eqn{\tau_j}. In this case \code{h} should be supplied. It is also possible
#' to use different values of \code{h} on either side of \eqn{\tau_j}, in which case \code{h} should be a vector of length 2.
#' If \code{cp_bound = TRUE}, then if there is another estimated changepoint within the \code{h} of \eqn{\tau_j} then the window size used
#' will be the minimum of \code{h} and the distance between \eqn{\tau_j} and the closest estimated changepoint to \eqn{\tau_j}.
#' \item There are no other changepoints between \eqn{\tau_{j-1}} and \eqn{\tau_{j+1}}. In this case \code{h} should be set to \code{NULL}
#' and \code{gamma} should be 1. (This is the default for the function.)
#' \item There are no other changepoints between the midpoint of \eqn{\tau_{j-1}} and \eqn{\tau_j} and the midpoint of \eqn{\tau_j} and 
#' \eqn{\tau_{j+1}}. In this case use \code{h = NULL} and \code{gamma = 0.5}.
#' }
#'
#' @return A list.
#' \itemize{
#' \item \code{p_value} A vector of p-values
#' \item \code{P_both} If \code{return_probs=TRUE}, a matrix containing values of \eqn{Pr(\phi \in S \& |\phi| > |\phi_{obs}|)} for each
#' \eqn{\psi}; otherwise \code{NA}.
#' \item \code{P_phi_in_S} If \code{return_probs=TRUE}, a matrix containing values of \eqn{Pr(\phi \in S)} for each \eqn{\psi};
#' otherwise \code{NA}.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(200) + c(rep(1,50), rep(-1,50), rep(1,50), rep(-1,50))
#' results <- binary_segmentation(y, threshold=4)
#' calculate_pvals(y, method="bs", results=results, h=10)
#'
#' y <- rnorm(200) + c(rep(1,50), rep(-1,50), rep(1,50), rep(-1,50))
#' results <- binary_segmentation(y, threshold=4)
#' calculate_pvals(y, method="bs", results=results, h=10)
#'
calculate_pvals <- function(y, method="bs", results=NULL, N=10, threshold=NULL, maxiter=NULL, h=NULL, gamma=1, cp_bound=TRUE, sigma2=1, eps0=0.01,
                            include_original=TRUE, num_pvals=NULL, random_samples=NULL, num_rand_samples=NULL, seeded=FALSE,
                            decay=NULL, return_probs=FALSE, autocor=FALSE){

  stopifnot( method %in% c("bs", "wbs", "not") )

  n <- length(y)

  # Implement changepoint algorithm
  if ( is.null(results) ){

    if ( method != "bs" & is.null(num_rand_samples) & is.null(random_samples) ){
      num_rand_samples <- 100
    } else if ( method != "bs" & is.null(num_rand_samples) ){
      num_rand_samples <- nrow(random_samples)
    }

    results <- find_cps(y, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_samples, rand_ints=random_samples,
                         seeded=seeded, decay=decay)

  }

  threshold <- results$threshold
  maxiter <- results$maxiter
  if ( method != "bs" ){
    random_samples <- results$rand_ints
  }
  b <- results$results$b[results$results$cp == 1]
  d <- results$results$d[results$results$cp == 1]

  if ( length(b) == 0 ){
    stop("No changepoints detected")
  } else if ( is.null(num_pvals) ){
    num_pvals <- length(b)
  } else if ( length(b) < num_pvals ){
    num_pvals <- length(b)
  }

  p_value <- rep(NA, num_pvals)
  P_both <- P_phi_in_S <- matrix(0, nrow=num_pvals, ncol=N)


  if ( autocor ){
    sig2 <- diag((1 - rho^seq(2, 2*n, by=2))/(1 - rho^2))
    for ( i in 1:(n - 1) ){
      sig2[(i + 1):n, i] <- sig2[i, (i + 1):n] <- rho^seq(1:(n - i)) * sig2[i, i]
    }
    D <- diag(svd(sig2)$d)
  }


  for ( jj in 1:num_pvals ){

    # Calculate nu
    if ( is.null(h) ){ # use changepoint on either side to determine h

      cps <- c(0, sort(b, decreasing=FALSE), n)
      j <- (1:length(cps))[cps==b[jj]]
      h1 <- ceiling( gamma*(cps[j] - cps[j-1]) )
      h2 <- ceiling( gamma*(cps[j+1] - cps[j]) )

    } else if ( length(h) >= 2 ) {

      if ( cp_bound ){
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
        if ( b[jj] - h[1] < 0 ){
          h1 <- b[jj]
          h2 <- h[2]
        } else if ( b[jj] + h[2] > n ){
          h1 <- h[1]
          h2 <- n - b[jj]
        } else {
          h1 <- h[1]
          h2 <- h[2]
        }
      }

    } else {

      stopifnot( h >= 2 )
      if ( cp_bound ){
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
      } else {
        if ( b[jj] - h < 0 ){
          h1 <- b[jj]
        } else {
          h1 <- h
        }
        if ( b[jj] + h > n ){
          h2 <- n - b[jj]
        } else {
          h2 <- h
        }
      }

    }


    nu <- c(rep(0, b[jj] - h1), rep(1/h1, h1), rep(-1/h2, h2), rep(0, n - b[jj] - h2))
    nu2 <- ifelse(autocor, c(t(nu) %*% D %*% nu), 1/h1 + 1/h2)
    nh <- h1 + h2    
    nuTy <- c(t(nu) %*% y)

    # Construct Y s.t. y_t = Y_t^T (1,phi,Psi)
    Z <- diag(nh) - matrix(1/nh, nrow=nh, ncol=nh) - 1/nu2 * crossprod(t(c(rep(1/h1, h1), rep(-1/h2, h2))))
    Z[abs(Z) < 10^(-10)] <- 0
    U <- svd(Z)$u[,1:(nh - 2)]
    U[abs(U) < 10^(-10)] <- 0

    Y <- cbind(y, matrix(0, nrow=n, ncol=nh - 1))
    Y[nu != 0, 1] <- mean(y[nu != 0]) # replace elements in (b-h, b+h) with mean
    Y[nu != 0, 2] <- 1/nu2 * nu[nu != 0] # constant of phi
    Y[nu != 0, 3:ncol(Y)] <- U # constants of Psi
    colnames(Y) <- c("y0", "phi", paste0("psi", 1:(nh-2)))
    
    # Generate psi values
    if ( include_original ){
      if ( N > 1 ){
        Psi <- rbind( matrix(rnorm((nh - 2) * (N - 1), sd=sqrt(sigma2)), nrow=(N - 1)), as.numeric(t(U) %*% y[(b[jj] - h1 + 1):(b[jj] + h2)]) )
      } else {
        Psi <- matrix(as.numeric(t(U) %*% y[(b[jj] - h1 + 1):(b[jj] + h2)]), nrow=1)
      }
    } else {
      Psi <- matrix(rnorm((nh - 2)*N, sd=sqrt(sigma2)), nrow=N)
    }


    for ( iter in 1:N ){

      # Calculate S
      y_new <- Y[,1] + nuTy / sum(nu^2) * nu # this is necessary b/c when calculating S it assumes we need to subtract it
      y_new[(b[jj] - h1 + 1):(b[jj] + h2)] <- y_new[(b[jj] - h1 + 1):(b[jj] + h2)] + U %*% t(Psi[iter,,drop=FALSE])

      r2 <- find_cps(y_new, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_samples, rand_ints=random_samples, seeded=seeded, decay=decay)
      b2 <- r2$results$b[r2$result$cp == 1]
      b2 <- b2[!is.na(b2)]
      if ( length(b2) >= 1 ){
        if ( !is.null(h) & (r2$results$b[1] == b[jj]) ){
          S <- calculate_S(y_new, results=r2, nu=nu, threshold=threshold, maxiter=maxiter, method=method, first_cp_only=TRUE, nuTy=nuTy, nu2=nu2, rand_ints=random_samples, seeded=seeded, decay=decay)
          # nuTy should be calculated using the original y!
        } else {
          S <- calculate_S(y_new, results=r2, nu=nu, nu2=nu2, threshold=threshold, maxiter=maxiter, method=method, nuTy=nuTy, rand_ints=random_samples, seeded=seeded, decay=decay)
        }
      } else {
        S <- calculate_S(y_new, results=r2, nu=nu, nu2=nu2, threshold=threshold, maxiter=maxiter, method=method, nuTy=nuTy, rand_ints=random_samples, seeded=seeded, decay=decay)
      }


      # Calculate probability

      if ( is.null(h) ){

        # Find which intervals contain the correct combination of changepoints
        ## First get rid of intervals which contain too many changepoints
        max_cps_found <- (ncol(S) - 2)/2
        if ( max_cps_found > length(results$changepoints) ){
          S2 <- S[ is.na(S[,(2+length(results$changepoints)+1)]), ]
        } else {
          S2 <- S
        }
        ## Get rid of intervals which contain too few changepoints
        S2 <- S2[ !is.na(S2[,(2+length(results$changepoints))]), ]
        ## Check we have the correct combination of changepoints
        z <- rep(0, nrow(S2))
        cps <- sort(results$changepoints)
        for ( j in 1:nrow(S2) ){
          if ( sum(abs(sort(S2[j, 3:(length(cps) + 2)]) - cps)) < 10^(-10) ){
            z[j] <- 1
          }
        }
        S2 <- S2[z==1,]

      } else {

        # Find intervals which contain b[jj]
        max_cps_found <- (ncol(S) - 2)/2
        S2 <- S[ S[,3]==b[jj], ]
        if ( max_cps_found > 1 ){
          for ( j in 2:max_cps_found ){
            S2 <- rbind(S2, S[ S[,j+2]==b[jj], ])
          }
        }
        S2 <- S2[ !is.na(S2$b1), ]

      }

      if ( autocor ){
        sigma_phi <- sqrt(c(t(nu) %*% sig2 %*% nu) * sigma2)
      } else {
        sigma_phi <- sqrt(nu2 * sigma2)
      }

      # Calculate P(phi \in S)
#      P_phi_in_S[jj, iter] <- sum( pnorm(S2[,2] / sqrt(nu2 * sigma2)) - pnorm(S2[,1] / sqrt(nu2 * sigma2)) )
      P_phi_in_S[jj, iter] <- sum(pnorm(S2[,2] / sigma_phi) - pnorm(S2[,1] / sigma_phi))

      # Calculate P(phi > |nuTy| & phi \in S)
      P_both[jj, iter] <- 0
      if ( nrow(S2) >= 1 ){
        for ( i in 1:nrow(S2) ){
          if ( S2$upper_lim[i] > abs(nuTy) ){
            P_both[jj, iter] <- P_both[jj, iter] +
              pnorm(S2$upper_lim[i] / sigma_phi) - pnorm(max( c(abs(nuTy), S2$lower_lim[i]) ) / sigma_phi)
          }
          if ( S2$lower_lim[i] < (-1)*abs(nuTy) ){
            P_both[jj, iter] <- P_both[jj, iter] +
              pnorm(min( c(-abs(nuTy), S2$upper_lim[i]) ) / sigma_phi) - pnorm(S2$lower_lim[i] / sigma_phi)
          }
        }
      }

    }

    p_value[jj] <- sum(P_both[jj,]) / sum(P_phi_in_S[jj,])
    
    print(paste0("Calculated ", jj, "th p-value"))

  }

  if ( !(return_probs) ){
    P_both <- NA
    P_phi_in_S <- NA
  }

  return( list(p_value=p_value, P_both=P_both, P_phi_in_S=P_phi_in_S) )

}
