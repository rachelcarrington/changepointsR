#' Calculate p-values for change in mean model
#'
#' @param y Numeric vector of data
#' @param method Method used to generate changepoints. Options: "bs" for binary segmentation,
#' "wbs" for wild or seeded binary segmentation, "not" for narrowest over threshold.
#' @param results Output of changepoint algorithm.
#' @param N Positive integer. Number of random psi's to generate.
#' @param nus List of nu's for each changepoint.
#' @param threshold Scalar. Minimum threshold for changepoint detection. Ignored if results is specified.
#' @param maxiter Positive integer. Maximum number of changepoints to find in changepoint algorithm. Ignored if
#' results is specified
#' @param h Positive integer (>=2) or NULL. Window size. If NULL then the changepoints either side of the
#' changepoint of interest are used to define the window.
#' @param sigma2 Scalar. Variance of y.
#' @param eps0 Scalar. Hyperparameter for calculating S.
#' @param include_original Logical. Whether to include the psi value corresponding to the observed data in place of
#' one of the ranom samples.
#' @param num_pvals Integer. Number of changepoints to calculate p-values for.
#' @param random_samples Matrix containing random intervals. For wild binary segmentation and narrowest over
#' threshold only.
#' @param num_rand_samples Number of random intervals to use for WBS/NOT. Ignored if random_samples is specified,
#' method = "bs", or seeded = TRUE.
#' @param seeded Logical. If TRUE, use seeded binary segmentation rather than wild binary segmentation. Only used if
#' method = "wbs".
#' @param decay Decay parameter for seeded binary segmentation. Only used if method = "wbs" and seeded = TRUE.
#'
#' @return A vector of p-values
#' @export
#'
#' @examples
#' # to do
#'
calculate_pvals <- function(y, method="bs", results=NULL, N=100, nus=NULL, threshold=NULL, maxiter=NULL, h=2,
                              sigma2=1, eps0=0.01, include_original=FALSE, num_pvals=NULL, random_samples=NULL,
                              num_rand_samples=NULL, seeded=FALSE, decay=NULL){

  ### Calculate p-values for all detected changepoints, for BS, W/SBS, or NOT

  #### if num_pvals > 1, nus should be a list

  n <- length(y)

  ## Implement changepoint algorithm
  if ( is.null(results) ){

    if ( method!="bs" & is.null(num_rand_samples) & is.null(random_samples) ){
      num_rand_samples <- 100
    } else if ( method!="bs" & is.null(num_rand_samples) ){
      num_rand_samples <- nrow(random_samples)
    }

    results <- find_cps( y, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_samples, rand_ints=random_samples,
                         seeded=seeded, decay=decay )

  }

  threshold <- results$threshold
  maxiter <- results$maxiter
  if ( method!="bs"){
    random_samples <- results$rand_ints
  }
  b <- results$results$b[ results$results$cp==1 ]
  d <- results$results$d[ results$results$cp==1 ]

  if ( length(b)==0 ){
    stop("No changepoints detected")
  } else if ( is.null(num_pvals) ){
    num_pvals <- length(b)
  } else if ( length(b) < num_pvals ){
    num_pvals <- length(b)
  }

  if ( is.null(nus) ){
    nus <- as.list( rep(NA, num_pvals) )
  }

  ## Estimate sigma^2 if not given (not sure if this is correct?)
  if (is.null(sigma2)){
    z <- 0
    b_sorted <- c(0, sort(b), n)
    for ( k in 1:(length(b)+1) ){
      mean_est <- mean(y[(b_sorted[k]+1):b_sorted[k+1]])
      z <- z + sum( (y[(b_sorted[k]+1):b_sorted[k+1]] - mean_est)^2 )
    }
    sigma2 <- z / n
  }


  Z <- diag(2*h) - matrix( 1/(2*h), nrow=2*h, ncol=2*h ) - h/2 * crossprod(t(c(rep(1/h,h), rep(-1/h,h))))
  Z[ abs(Z) < 10^(-10) ] <- 0
  U <- svd(Z)$u[,1:(2*h-2)]
  U[ abs(U) < 10^(-10) ] <- 0


  p_value <- rep(NA, num_pvals)
  P_both <- P_phi_in_S <- matrix(0, nrow=num_pvals, ncol=N)


  for ( jj in 1:num_pvals ){

    if ( b[jj] >= h & b[jj] <= n-h ){

      if (is.na(nus[[jj]])){
        if (b[jj]-h < 0){
          nu <- c(rep(1/b[jj],b[jj]), rep(-1/h,h), rep(0,n-b[jj]-h))
          nu2 <- 1/b[jj] + 1/h
        } else if (b[jj]+h > n){
          nu <- c(rep(0,b[jj]-h), rep(1/h,h), rep(-1/(n-b[jj]),n-b[jj]))
          nu2 <- 1/h + 1/(n-b[jj])
        } else {
          nu <- c(rep(0,b[jj]-h), rep(1/h,h), rep(-1/h,h), rep(0,n-b[jj]-h))
          nu2 <- 2/h
        }
      } else {
        nu <- nus[[jj]]
        nu2 <- sum(nu^2)
      }

      nuTy <- as.numeric( t(nu) %*% y )

      #### Construct Y s.t. y_t = Y_t^T (1,phi,Psi)
      Y <- cbind(y,matrix(0,nrow=n,ncol=2*h-1))
      Y[(b[jj]-h+1):(b[jj]+h),1] <- mean(y[(b[jj]-h+1):(b[jj]+h)]) ## replace elements in (b-h,b+h) with mean
      Y[(b[jj]-h+1):(b[jj]+h),2] <- 1/sum(nu^2)* nu[(b[jj]-h+1):(b[jj]+h)] ## constant of phi
      Y[(b[jj]-h+1):(b[jj]+h),3:ncol(Y)] <- U ## constants of Psi
      colnames(Y) <- c("y0","phi",paste0("psi",1:(2*h-2)))


      if ( include_original ){
        if ( N > 1 ){
          Psi <- rbind( matrix( rnorm((2*h-2)*(N-1), sd=sqrt(sigma2)), nrow=(N-1) ), as.numeric(t(U) %*% y[(b[jj]-h+1):(b[jj]+h)]) )
        } else {
          Psi <- matrix( as.numeric(t(U) %*% y[(b[jj]-h+1):(b[jj]+h)]), nrow=1 )
        }
      } else {
        Psi <- matrix( rnorm((2*h-2)*N, sd=sqrt(sigma2)), nrow=N )
      }


      ### Sampling
      for ( iter in 1:N ){

        y_new <- Y[,1] + nuTy / nu2 * nu ## this is necessary b/c when calculating S it assumes we need to subtract it
        y_new[(b[jj]-h+1):(b[jj]+h)] <- y_new[(b[jj]-h+1):(b[jj]+h)] + U %*% t(Psi[iter,,drop=FALSE])

        r2 <- find_cps( y_new, method=method, threshold=threshold, maxiter=maxiter, num_rand_ints=num_rand_samples, rand_ints=random_samples, seeded=seeded, decay=decay )
        b2 <- r2$results$b[ r2$result$cp==1 ]
        b2 <- b2[!is.na(b2)]
        if ( length( b2 ) >= 1 ){
          if ( r2$results$b[1] == b[jj] ){
            S <- calculate_S_all_methods( y_new, results=r2, nu=nu, threshold=threshold, maxiter=maxiter, method=method, first_cp_only=TRUE, nuTy=nuTy, rand_ints=random_samples, seeded=seeded, decay=decay )
            ## nuTy should be calculated using the original y!
          } else {
            S <- calculate_S_all_methods( y_new, results=r2, nu=nu, threshold=threshold, maxiter=maxiter, method=method, nuTy=nuTy, rand_ints=random_samples, seeded=seeded, decay=decay )
          }
        } else {
          S <- calculate_S_all_methods( y_new, results=r2, nu=nu, threshold=threshold, maxiter=maxiter, method=method, nuTy=nuTy, rand_ints=random_samples, seeded=seeded, decay=decay )
        }


        ### Calculate probability

        ### Take intervals which are for the correct values of b
        max_cps_found <- (ncol(S) - 2)/2
        S2 <- S[ S[,3]==b[jj], ]
        if ( max_cps_found > 1 ){
          for ( j in 2:max_cps_found ){
            S2 <- rbind( S2, S[ S[,j+2]==b[jj], ] )
          }
        }
        S2 <- S2[ !is.na(S2$b1), ]

        ### Calculate P(phi \in S)
        P_phi_in_S[jj,iter] <- sum( pnorm( S2[,2] / sqrt(nu2*sigma2) ) - pnorm( S2[,1] / sqrt(nu2*sigma2) ) )

        ### Calculate P(phi > |nuTy| & phi \in S)
        P_both[jj,iter] <- 0
        for ( i in 1:nrow(S2) ){
          if ( S2$upper_lim[i] > abs(nuTy) ){
            P_both[jj,iter] <- P_both[jj,iter] +
              pnorm( S2$upper_lim[i] / sqrt(nu2*sigma2) ) - pnorm( max( c( abs(nuTy), S2$lower_lim[i] ) ) / sqrt(nu2*sigma2) )
          }
          if ( S2$lower_lim[i] < (-1)*abs(nuTy) ){
            P_both[jj,iter] <- P_both[jj,iter] +
              pnorm( min( c( -abs(nuTy), S2$upper_lim[i] ) ) / sqrt(nu2*sigma2) ) - pnorm( S2$lower_lim[i] / sqrt(nu2*sigma2) )
          }

        }

      }

      p_value[jj] <- sum(P_both[jj,]) / sum(P_phi_in_S[jj,])

    }

  }

  return( p_value )

}


