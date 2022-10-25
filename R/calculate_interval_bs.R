#' Calculate interval for binary segmentation
#'
#' @description Find values of phi which satisfy the required inequalities so that the binary segmentation algorithm returns (b, d);
#' Option to only consider part of (b, d), e.g. if we just want b1 & d1 to be specified, & don't care about
#' later values
#'
#' @details (possibly outdated)
#' At present this ignores the threshold and just matches the number of CPs found (currently updating)
#' If threshold is also supplied, this will also be used.
#' Otherwise, if threshold is supplied (but not n.cp), the interval returned will be that for which the complete
#' b and d match those given, when BS is given this threshold. If threshold is not supplied, it will be assumed
#' that BS is run with a fixed number of iterations. If both n.cp and threshold are supplied, it will be assumed that
#' n.cp is the maximum number of iterations, but the threshold will also be used as a minimum for C(t).
#'
#' @param y Numeric vector of data.
#' @param nu ...
#' @param b Vector of changepoints detected by binary segmentation algorithm.
#' @param d Directions of changepoints detected by binary segmentation algorithm. Vector whose entries are all equal to \code{1} or \code{-1}.
#' @param nu2 Value of the squared 2-norm of \code{nu}.
#' @param nuTy Value of \code{nu^T %*% y}.
#' @param threshold Minimum changepoint threshold used in binary segmentation algorithm.
#' @param n.cp Maximum number of changepoints to detect of binary segmentation algorithm.
#'
#' @return A 2-dimensional vector
#' @export
#'
#' @examples
#' # to do
#'
calculate_interval_bs <- function(y, nu, b, d, nu2=NULL, nuTy=NULL, threshold=NULL, n.cp=NULL ){

  if ( is.null(n.cp) & is.null(threshold) ){
    stop("At least one of n.cp and threshold must be supplied.")
  } else if ( is.null(n.cp) ){
    n.cp <- length(y)
  }

  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }

  if ( is.null(nuTy) ){
    nuTy <- as.numeric(t(nu)%*%y)
  }

  n <- length(y)

  cs <- cusum_phi_vec(y, nu, nu2, nuTy)

  ### If no CPs (threshold must be specified):
  if ( length(b)==0 ){

    ### We need all |Ct|'s to be < threshold
    ### We can ignore those for which cs[,2] = 0, since these are constant in phi
    x <- ( abs(cs[,2]) > 10^(-10) )
    inequalities <- c( (threshold - cs[x,1]) / cs[x,2], (-threshold - cs[x,1]) / cs[x,2] )
    signs <- c( ifelse( cs[x,2] > 0, -1, 1 ), ifelse( cs[x,2] > 0, 1, -1 ) )
    max_lower_bound <- max( inequalities[ signs==1] )
    min_upper_bound <- min( inequalities[ signs==-1] )

  } else if ( length(b)>=n.cp ){
    ## (in this case we ignore the threshold, except to check whether Ct at each CP is larger than it)

    ### For the first CP:
    #### (-d1) * C(b1) > threshold (0 if not given)
    if ( abs( cs[b[1],2] ) > 10^(-10) ){
      if ( is.null(threshold) ){
        inequalities <- (-1) * cs[b[1],1] / cs[b[1],2]
      } else {
        inequalities <- (-d[1]*threshold - cs[b[1],1]) / cs[b[1],2]
      }
      signs <- ifelse( cs[b[1],2] > 0, -d[1], d[1] )
      if ( signs == 1 ){
        max_lower_bound <- inequalities
        min_upper_bound <- Inf
      } else {
        max_lower_bound <- -Inf
        min_upper_bound <- inequalities
      }
    } else {
      max_lower_bound <- -Inf
      min_upper_bound <- Inf
    }

    #### | C( b1 ) | > +/- C( t ) for t \neq b1
    inequalities <- cbind( (cs[-b[1],1] - cs[b[1],1])/(-cs[-b[1],2] + cs[b[1],2] ),
                           (-1)*(cs[-b[1],1] + cs[b[1],1])/(cs[-b[1],2] + cs[b[1],2]) )
    signs <- cbind( ifelse( -cs[-b[1],2] + cs[b[1],2] > 0, -d[1], d[1] ),
                    ifelse( cs[-b[1],2] + cs[b[1],2] > 0, -d[1], d[1] ) )

    #### Remove entries where we have divided by 0
    signs[ abs(-cs[-b[1],2] + cs[b[1],2]) <= 10^(-10), 1 ] <- NA
    signs[ abs(cs[-b[1],2] + cs[b[1],2]) <= 10^(-10), 2 ] <- NA

    #### Update bounds
    max_lower_bound <- max( max_lower_bound, max( inequalities[ signs==1 ], na.rm=TRUE ) )
    min_upper_bound <- min( min_upper_bound, min( inequalities[ signs==-1 ], na.rm=TRUE ) )


    ### For later CPs, if any:
    if ( n.cp > 1 ){

      for ( k in 2:n.cp ){

        #### Split data at CPs
        previous_cps <- sort(b[1:(k-1)])
        s1 <- c(1, previous_cps + 1)
        e1 <- c(previous_cps, n)

        #### Calculate CUSUM statistic in terms of phi, for each interval
        cs2 <- matrix( NA, nrow=n-1, ncol=2 )
        for ( m in 1:length(s1) ){
          if ( e1[m] > s1[m] ){  ##(intervals may have length 0 if we have consecutive CPs)
            cs2[ s1[m]:(e1[m]-1), ] <- cusum_phi_vec( y, nu, nu2=nu2, nuTy=nuTy, s=s1[m], e=e1[m] )
          }
        }

        #### We need -d[k]*C(b[k]) to be greater than 0 or the threshold
        if ( abs( cs2[b[k],2] ) > 10^(-10) ){
          if ( is.null(threshold) ){
            inequalities <- (-1) * cs2[b[k],1] / cs2[b[k],2]
          } else {
            inequalities <- (-d[k]*threshold - cs2[b[k],1]) / cs2[b[k],2]
          }
          signs <- ifelse( cs2[b[k],2] > 0, -d[k], d[k] )
          if ( signs==1 ){
            max_lower_bound <- max( max_lower_bound, inequalities )
          } else {
            min_upper_bound <- min( min_upper_bound, inequalities )
          }
        }

        #### Also, need |C(b[k])| > |C(b[t])| for t \neq k
        inequalities <- cbind( (cs2[-b[k],1] - cs2[b[k],1])/(-cs2[-b[k],2] + cs2[b[k],2] ),
                               (-1)*(cs2[-b[k],1] + cs2[b[k],1])/(cs2[-b[k],2] + cs2[b[k],2]) )
        signs <- cbind( ifelse( -cs2[-b[k],2] + cs2[b[k],2] > 0, -d[k], d[k] ),
                        ifelse( cs2[-b[k],2] + cs2[b[k],2] > 0, -d[k], d[k] ) )

        #### Remove entries where we have divided by 0
        signs[ abs(-cs2[-b[k],2] + cs2[b[k],2]) <= 10^(-10), 1 ] <- NA
        signs[ abs(cs2[-b[k],2] + cs2[b[k],2]) <= 10^(-10), 2 ] <- NA

        #### Remove NAs & recalculate bounds
        max_lower_bound <- max( max_lower_bound, max(inequalities[ signs==1 ], na.rm=TRUE) )
        min_upper_bound <- min( min_upper_bound, min(inequalities[ signs==-1 ], na.rm=TRUE) )
      }
    }

  } else {
    ## length(b) < n.cp; we need to check that after first length(b) CPs, all |Ct|'s are smaller than lambda

    ### For the first CP:
    #### (-d1) * C(b1) > threshold (0 if not given)
    if ( abs(cs[b[1],2]) > 10^(-10) ){
      if ( is.null(threshold) ){
        inequalities <- (-1) * cs[b[1],1] / cs[b[1],2]
      } else {
        inequalities <- (-d[1]*threshold - cs[b[1],1]) / cs[b[1],2]
      }
      signs <- ifelse( cs[b[1],2] > 0, -d[1], d[1] )
      if ( signs == 1 ){
        max_lower_bound <- inequalities
        min_upper_bound <- Inf
      } else {
        max_lower_bound <- -Inf
        min_upper_bound <- inequalities
      }
    } else {
      max_lower_bound <- -Inf
      min_upper_bound <- Inf
    }

    #### | C( b1 ) | > +/- C( t ) for t \neq b1
    inequalities <- cbind( (cs[-b[1],1] - cs[b[1],1])/(-cs[-b[1],2] + cs[b[1],2] ),
                           (-1)*(cs[-b[1],1] + cs[b[1],1])/(cs[-b[1],2] + cs[b[1],2]) )
    signs <- cbind( ifelse( -cs[-b[1],2] + cs[b[1],2] > 0, -d[1], d[1] ),
                    ifelse( cs[-b[1],2] + cs[b[1],2] > 0, -d[1], d[1] ) )

    #### Remove entries where we have divided by 0
    signs[ abs(-cs[-b[1],2] + cs[b[1],2]) <= 10^(-10), 1 ] <- NA
    signs[ abs(cs[-b[1],2] + cs[b[1],2]) <= 10^(-10), 2 ] <- NA

    #### Update bounds
    max_lower_bound <- max( max_lower_bound, max( inequalities[ signs==1 ], na.rm=TRUE ) )
    min_upper_bound <- min( min_upper_bound, min( inequalities[ signs==-1 ], na.rm=TRUE ) )


    ### For later CPs, if any:
    if ( length(b) > 1 ){

      for ( k in 2:length(b) ){

        #### Split data at CPs
        previous_cps <- sort(b[1:(k - 1)])
        s1 <- c(1, previous_cps + 1)
        e1 <- c(previous_cps, n)

        #### Calculate CUSUM statistic in terms of phi, for each interval
        cs2 <- matrix(NA, nrow=n-1, ncol=2)
        for ( m in 1:length(s1) ){
          if ( e1[m] > s1[m] ){  ##(intervals may have length 0 if we have consecutive CPs)
            cs2[ s1[m]:(e1[m]-1), ] <- cusum_phi_vec(y, nu, nu2=nu2, nuTy=nuTy, s=s1[m], e=e1[m])
          }
        }

        #### We need -d[k]*C(b[k]) to be greater than 0 or the threshold
        if ( abs(cs2[b[k],2]) > 10^(-10) ){
          if ( is.null(threshold) ){
            inequalities <- (-1) * cs2[b[k], 1] / cs2[b[k], 2]
          } else {
            inequalities <- (-d[k]*threshold - cs2[b[k], 1]) / cs2[b[k], 2]
          }
          signs <- ifelse( cs2[b[k],2] > 0, -d[k], d[k] )
          if ( signs==1 ){
            max_lower_bound <- max(max_lower_bound, inequalities)
          } else {
            min_upper_bound <- min(min_upper_bound, inequalities)
          }
        }

        #### Also, need |C(b[k])| > |C(b[t])| for t \neq k
        inequalities <- cbind( (cs2[-b[k],1] - cs2[b[k],1])/(-cs2[-b[k],2] + cs2[b[k],2]),
                               (-1)*(cs2[-b[k],1] + cs2[b[k],1])/(cs2[-b[k],2] + cs2[b[k],2]) )
        signs <- cbind( ifelse( -cs2[-b[k],2] + cs2[b[k],2] > 0, -d[k], d[k] ),
                        ifelse( cs2[-b[k],2] + cs2[b[k],2] > 0, -d[k], d[k] ) )

        #### Remove entries where we have divided by 0
        signs[ abs(-cs2[-b[k],2] + cs2[b[k],2]) <= 10^(-10), 1 ] <- NA
        signs[ abs(cs2[-b[k],2] + cs2[b[k],2]) <= 10^(-10), 2 ] <- NA

        #### Remove NAs & recalculate bounds
        max_lower_bound <- max( max_lower_bound, max(inequalities[ signs==1 ], na.rm=TRUE) )
        min_upper_bound <- min( min_upper_bound, min(inequalities[ signs==-1 ], na.rm=TRUE) )
      }
    }

    ### Check that for the next k, all |Ct|'s are below the threshold
    k <- length(b) + 1

    #### Split data at CPs
    previous_cps <- sort(b[1:(k-1)])
    s1 <- c(1, previous_cps + 1)
    e1 <- c(previous_cps, n)

    #### Calculate CUSUM statistic in terms of phi, for each interval
    cs2 <- matrix( NA, nrow=n-1, ncol=2 )
    for ( m in 1:length(s1) ){
      if ( e1[m] > s1[m] ){  ##(intervals may have length 0 if we have consecutive CPs)
        cs2[ s1[m]:(e1[m]-1), ] <- cusum_phi_vec( y, nu, nu2=nu2, nuTy=nuTy, s=s1[m], e=e1[m] )
      }
    }

    x <- abs(cs2[,2]) > 10^(-10)
    inequalities <- c( (threshold - cs2[x,1]) / cs2[x,2], (-threshold - cs2[x,1]) / cs2[x,2] )
    signs <- c( ifelse( cs2[x,2] > 0, -1, 1 ), ifelse( cs2[x,2] > 0, 1, -1 ) )
    max_lower_bound <- max( c(max_lower_bound, inequalities[ signs==1 ]), na.rm=TRUE  )
    min_upper_bound <- min( c(min_upper_bound, inequalities[ signs==-1 ]), na.rm=TRUE  )

  }

  ### Inequalities are all satisfied by phi s.t. max_lower_bound < phi < min_upper_bound

  return( c( max_lower_bound, min_upper_bound ) )
}
