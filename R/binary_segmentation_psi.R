#' Post-selection inference for binary segmentation
#'
#' @description Calculate p-values associated with detected changepoints for binary segmentation;
#' uses method of Jewell et al.
#'
#' @details
#' ...
#'
#' @param y A numeric vector of data.
#' @param results Output of binary_segmentation function; can be NULL.
#' @param nus List of nu's for each detected changepoint.
#' @param threshold Minimum threshold for CUSUM statistic for changepoint detection; ignored if results specified.
#' @param maxiter Max. number of changepoints to find; ignored if results specified; otherwise defaults to n-1.
#' @param eps0 Hyperparameter for calculating S.
#' @param sigma2 Variance of y. If unknown, it can be estimated within the function.
#' @param h Window size. If NULL, the window is taken to be ...
#' @param first_cp_only Logical. Should be TRUE if the condition is that the changepoint of interest is in the model;
#'  FALSE if the condition is that all changepoints are the same.
#' @param num_pvals Integer. Maximum number of p-values to calculate; defaults to n-1 if NULL.
#'
#' @return A list containing: the vector of changepoints and the vector of p-values.
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
  } else if ( is.null(num_pvals) ){
    num_pvals <- length(b)
  } else if ( length(b) < num_pvals ){
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

    nuTy <- as.numeric( t(nu) %*% y )

    if ( jj == 1 ){
      S <- calculate_S_all_methods( y, results=results, nu=nu, threshold=threshold, maxiter=maxiter, method="bs",
                                    nuTy=nuTy, first_cp_only=first_cp_only )
    } else {
      S <- calculate_S_all_methods( y, results=results, nu=nu, threshold=threshold, maxiter=maxiter, method="bs",
                                    nuTy=nuTy )
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
