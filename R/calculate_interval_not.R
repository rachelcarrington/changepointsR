#' Calculate interval for narrowest over threshold
#'
#' @description Find values of phi which satisfy the required inequalities so that applying narrowest over threshold to \eqn{y'(\phi)}
#' returns \code{(b, d)}.
#'
#' @details Used inside \code{calculate_S} if \code{method = "not"}.
#'
#' @param y Numeric vector of data.
#' @param nu Numeric vector.
#' @param results Output of \code{narrowest_over_threshold}.
#' @param nu2 Value of \eqn{||\nu||_2^2}.
#' @param nuTy Value of \eqn{\nu^T y}.
#'
#' @return A 2-dimensional vector.
#'
#' @export
#'
#' @examples
#' set.seed(100)
#' y <- rnorm(100) + c(rep(1,50), rep(-1,50))
#' results <- narrowest_over_threshold(y, lambda=4, N=50)
#' b <- results$results$b
#' h <- 10
#' nu <- c(rep(0, b[1] - h), rep(1/h, h), rep(-1/h, h), rep(0, length(y) - b[1] - h))
#' calculate_interval_not(y, nu, results=results)
#'
calculate_interval_not <- function(y, nu, results, nu2=NULL, nuTy=NULL){

  n <- length(y)

  b <- results$results$b
  d <- results$results$d
  s <- results$results$s
  e <- results$results$e

  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }
  if ( is.null(nuTy) ){
    nuTy <- as.numeric(t(nu) %*% y)
  }
  
  rand_ints <- results$rand_ints
  threshold <- results$threshold

  cps <- c(0, sort(b, decreasing=FALSE), n)
  j <- (1:length(cps))[cps == b[1]]
  cps2 <- cps[(j - 1):(j + 1)]

  if ( nrow(results$results) == 0 ){
    # If there are no CPs detected, we require |C_(s,e) (t)| < threshold for all (s,e) and all t
    max_lower_bound <- -Inf
    min_upper_bound <- Inf

    for ( j in 1:nrow(rand_ints) ){  
      cs <- cusum_phi_vec(y, nu, nu2=nu2, nuTy=nuTy, s=rand_ints[j,1], e=rand_ints[j,2])
      cs <- cs[abs(cs[,2]) > 10^(-10),, drop=FALSE]
      inequalities <- c(threshold - cs[,1], -threshold - cs[,1]) / cs[,2]
      signs <- c(ifelse(cs[,2] > 0, -1, 1), ifelse(cs[,2] > 0, 1, -1))
      max_lower_bound <- max(c(max_lower_bound, inequalities[signs == 1]))
      min_upper_bound <- min(c(min_upper_bound, inequalities[signs == -1]))
    }

  } else {

    for ( k in 1:length(b) ){
      # If k > 1, we can delete all intervals containing previously detected changepoints, which were discarded by the algorithm
      # We can also remove intervals that are narrower than previous changepoint-containing intervals, which have already been dealt with here
      
      if ( k > 1 ){
        rand_ints <- rand_ints[!(b[j] >= rand_ints[,1] & b[j] < rand_ints[,2]),,drop=FALSE]
        rand_ints <- rand_ints[rand_ints[,2] - rand_ints[,1] >= e[k - 1] - s[k - 1],]
      }
      
      # For remaining narrower intervals, we need |C(t)| < threshold for all t
      cp_interval_width <- e[k] - s[k]
      smaller_intervals <- rand_ints[rand_ints[,2] - rand_ints[,1] < cp_interval_width,,drop=FALSE]

      if ( nrow(smaller_intervals) >= 1 ){
        inequalities <- signs <- numeric(0)
        for ( j in 1:nrow(smaller_intervals) ){
          s1 <- smaller_intervals[j,1]
          e1 <- smaller_intervals[j,2]
          cs2 <- cusum_phi_vec(y, nu, nu2, nuTy, s=s1, e=e1)
          cs2 <- cs2[abs(cs2[,2]) > 10^(-10),,drop=FALSE] # if cs2[t,2] = 0 then C(t) is constant in phi
          inequalities <- c(inequalities, c(threshold - cs2[,1], -threshold - cs2[,1]) / cs2[,2])
          signs <- c(signs, c(ifelse(cs2[,2] > 0, -1, 1), ifelse(cs2[,2] > 0, 1, -1)))
        }
        max_lower_bound <- max(c(max_lower_bound, inequalities[signs == 1]))
        min_upper_bound <- min(c(min_upper_bound, inequalities[signs == -1]))
      }
  
      # For intervals of the same width, we need |C(t)| < C_(s,e) (b1) for all t in interval
      # For the interval containing the changepoint, we need |C(t)| < C(b1) for t \neq b1
      
      same_intervals <- rand_ints[rand_ints[,2] - rand_ints[,1] == cp_interval_width,,drop=FALSE]

      if ( nrow(same_intervals) >= 1 ){
        cs <- cusum_phi_vec(y, nu, nu2=nu2, nuTy=nuTy, s=s[k], e=e[k])
        for ( j in 1:nrow(same_intervals) ){
          interval <- same_intervals[j,]
   
          # If this is the CP-containing interval, we need |C(b)| > |C(t)| for t \neq b and |C(b)| >= threshold
          if ( sum(abs(interval - c(s[k], e[k]))) == 0 ){
            
            # C(b[k]) > (-d1) * threshold
            if ( abs(cs[b[k] - s[k] + 1, 2]) > 10^(-10) ){ # if denominator is 0, inequality is always satisfied so we can ignore it
              inequalities <- (ifelse(d[k] == 1, -threshold, threshold) - cs[b[k] - s[k] + 1, 1]) / cs[b[k] - s[k] + 1, 2]
              signs <- ifelse(cs[b[k] - s[k] + 1, 2] > 0, -d[k], d[k])
              max_lower_bound <- max(c(max_lower_bound, inequalities[signs == 1]))
              min_upper_bound <- min(c(min_upper_bound, inequalities[signs == -1]))
            }

            # |C(b[k])| > +/- C(t) for t \neq b[k]
            inequalities <- cbind((cs[-b[k], 1] - cs[b[k] - s[k] + 1, 1]) / (cs[b[k] - s[k] + 1, 2] - cs[-b[k], 2]), 
                                    (-1) * (cs[-b[k], 1] + cs[b[k] - s[k] + 1, 1]) / (cs[b[k] - s[k] + 1, 2] + cs[-b[k], 2]))
            signs <- cbind(ifelse(cs[b[k] - s[k] + 1, 2] - cs[-b[k], 2] > 0, -d[k], d[k]), ifelse(cs[b[k] - s[k] + 1, 2] + cs[-b[k], 2] > 0, -d[k], d[k]))

            # Remove entries where we have divided by 0
            signs[abs(cs[,2] - cs[b[k] - s[k] + 1, 2]) <= 10^(-10), 1] <- NA
            signs[abs(cs[,2] + cs[b[k] - s[k] + 1, 2]) <= 10^(-10), 2] <- NA
            
            # Update bounds
            max_lower_bound <- max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE)
            min_upper_bound <- min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE)

          # Otherwise (if this is not the CP-containing interval), we need |C(b)| > +/- C(t) for all t
          } else {
            cs2 <- cusum_phi_vec(y, nu, nu2, nuTy, s=interval[1], e=interval[2])
            inequalities <- cbind((cs2[,1] - cs[b[k] - s[k] + 1, 1]) / (cs[b[k] - s[k] + 1, 2] - cs2[,2]), 
                                    (-1) * (cs2[,1] + cs[b[k] - s[k] + 1, 1]) / (cs[b[k] - s[k] + 1, 2] + cs2[,2]))
            signs <- cbind(ifelse(cs[b[k] - s[k] + 1, 2] - cs2[,2] > 0, -d[k], d[k]), ifelse(cs[b[k] - s[k] + 1, 2] + cs2[,2] > 0, -d[k], d[k]))

            # Remove entries where we have divided by 0
            signs[abs(cs2[,2] - cs[b[k] - s[k] + 1, 2]) <= 10^(-10), 1] <- NA
            signs[abs(cs2[,2] + cs[b[k] - s[k] + 1, 2]) <= 10^(-10), 2] <- NA

            # Update bounds
            max_lower_bound <- max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE)
            min_upper_bound <- min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE)
          }
        }
      }
    }

    # If we have less than the maximum number of changepoints, make sure no more changepoints are detected
    if ( length(b) < results$maxiter ){
      k <- length(b) + 1
      
      # Remove intervals containing previous changepoint, or narrower than previous CP-containing interval
      rand_ints <- rand_ints[!(b[j] >= rand_ints[,1] & b[j] < rand_ints[,2]),,drop=FALSE]
      rand_ints <- rand_ints[rand_ints[,2] - rand_ints[,1] >= e[k - 1] - s[k - 1],]

      # If we have any remaining intervals, we need |C(t)| < threshold on them for all t
      if ( nrow(rand_ints) >= 1 ){
        inequalities <- signs <- numeric(0)
        for ( j in 1:nrow(remaining_intervals) ){  
          cs <- cusum_phi_vec(y, nu, nu2=nu2, nuTy=nuTy, s=rand_ints[j,1], e=rand_ints[j,2])
          cs <- cs[abs(cs[,2]) > 10^(-10),,drop=FALSE]
          inequalities <- c(inequalities, (threshold - cs[,1]) / cs[,2], (-threshold - cs[,1]) / cs[,2])
          signs <- c(signs, ifelse(cs[,2] > 0, -1, 1), ifelse(cs[,2] > 0, 1, -1))
        }
        max_lower_bound <- max(c(max_lower_bound, inequalities[signs == 1]), na.rm=TRUE)
        min_upper_bound <- min(c(min_upper_bound, inequalities[signs == -1]), na.rm=TRUE)
      }
    }
  }
  return( c(max_lower_bound, min_upper_bound) )
}