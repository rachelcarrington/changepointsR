#' Calculate interval for narrowest over threshold
#'
#' @description Find values of phi which satisfy the required inequalities so that the narrowest over threshold algorithm 
#' returns \code{(b, d)};
#' Option to only consider part of \code{(b, d)}, e.g. if we just want b1 & d1 to be specified, & don't care about
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
#' @param y vector of data (ignored if cs & nuTy are supplied)
#' @param nu ...
#' @param results Output of \code{wild_binary_segmentation}.
#' @param nu2 value of ||nu||_2^2
#' @param nuTy value of nu^T y
#'
#' @return A 2-dimensional vector
#' @export
#'
#' @examples
#' x <- 0
#'
calculate_interval_not <- function(y, nu, results, nu2=NULL, nuTy=NULL){
  
  ######## Find values of phi which satisfy the required inequalities so that the NOT (change in mean) algorithm returns given results

  n <- length(y)

  b <- results$results$b
  d <- results$results$d
  s <- results$results$s
  e <- results$results$e

  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }
  if ( is.null(nuTy) ){
    nuTy <- as.numeric( t(nu)%*%y )
  }
  
  rand_ints <- results$rand_ints
  lambda <- results$threshold

  cps <- c( 0, sort(b, decreasing=FALSE), n )
  j <- (1:length(cps))[ cps==b[1] ]
  cps2 <- cps[ (j-1):(j+1) ]

  inequalities_list <- inequalities_signs <- numeric(0)

  if ( nrow( results$results ) == 0 ){

  ### If there are no CPs detected, we require |C_(s,e) (t)| < lambda for all (s,e) and all t

    for ( j in 1:nrow(rand_ints) ){  
      cs <- cusum_phi_vec(y, nu, nu2=nu2, nuTy=nuTy, s=rand_ints[j,1], e=rand_ints[j,2])
      cs <- cs[ abs(cs[,2]) > 10^(-10),,drop=FALSE ]
      inequalities_list <- c( inequalities_list, (lambda - cs[,1])/cs[,2], ((-1)*lambda - cs[,1])/cs[,2] )
      inequalities_signs <- c( inequalities_signs, ifelse( cs[,2] > 0, -1, 1 ), ifelse( cs[,2] > 0, 1, -1 ) )
    }

  } else {
    ## If >=1 CP detected
  
    ### Inequalities:
    	
    #### For narrower intervals, need |C(t)| < lambda for all t in interval, discarding intervals containing previous CPs

    cp_widths <- e - s + 1
    widths <- rand_ints[,2] - rand_ints[,1] + 1
    for ( k in 1:length(b) ){  
      if ( k == 1 ){
        smaller_intervals <- rand_ints[ widths < cp_widths[k],,drop=FALSE ]
      } else { ## if k > 1
        smaller_intervals <- rand_ints[ widths < cp_widths[k] & widths >= cp_widths[k-1],,drop=FALSE ] 
          ## we don't need to consider intervals narrower than previous CP-containing intervals as these have already been checked
        for ( j in 1:k ){
          smaller_intervals <- smaller_intervals[ !( b[j] >= smaller_intervals[,1] & b[j] < smaller_intervals[,2] ),,drop=FALSE ] 
            ## delete intervals containing previous CPs
        }
      }

      if ( nrow(smaller_intervals) >= 1 ){
        for ( j in 1:nrow(smaller_intervals) ){
          s1 <- smaller_intervals[j,1]
          e1 <- smaller_intervals[j,2]

          cs2 <- cusum_phi_vec(y, nu, nu2, nuTy, s=s1, e=e1)
          cs2 <- cs2[ abs(cs2[,2]) > 10^(-10),,drop=FALSE ] ## if cs2[t,2] = 0 then C(t) is constant in phi
          inequalities_list <- c( inequalities_list, (lambda - cs2[,1])/cs2[,2], ((-1)*lambda - cs2[,1])/cs2[,2] )
          inequalities_signs <- c( inequalities_signs, (-1)*(2*(cs2[,2] > 0) - 1), 1*(2*(cs2[,2] > 0) - 1) )
        }
      }
  
    }
  
    #### For intervals of same width, need |C(t)| < C_(s,e) (b1) for all t in interval
    #### For interval containing changepoint, need |C(t)| < C(b1) for t \neq b1

    for ( k in 1:length(b) ){
      width <- e[k] - s[k] + 1
      widths <- rand_ints[,2] - rand_ints[,1] + 1
      same_intervals <- rand_ints[ widths==width,,drop=FALSE ]
  
      ### Delete intervals containing earlier CPs
      if ( k >= 2 & nrow(same_intervals) >= 1 ){
        x <- rep(1, nrow(same_intervals))
        for ( j in 1:nrow(same_intervals) ){
          for ( m in 1:(k-1) ){
            if ( b[m] >= same_intervals[j,1] & b[m] < same_intervals[j,2] ){
              x[j] <- 0
            }
          }
        }
        same_intervals <- same_intervals[ x==1,,drop=FALSE ]
      }

      if ( nrow(same_intervals) >= 1 ){
        cs <- cusum_phi_vec(y, nu, nu2=nu2, nuTy=nuTy, s=s[k], e=e[k])
        for ( j in 1:nrow(same_intervals) ){
          interval <- same_intervals[j,]
   
          ### If this is the CP-containing interval, we need |C(b)| > +/- C(t) for t \neq b
          ### and | C(t) | >= lambda
          if ( sum( abs(interval - c(s[k], e[k])) )==0 ){
            ### C(b[k]) > (-d1)*lambda
            if ( abs(cs[b[k] - s[k] + 1, 2]) > 10^(-10) ){ ## if denominator is 0, inequality is always satisfied so we can ignore it
              if ( d[k] == 1 ){
                inequalities_list <- c( inequalities_list, (-1)*(lambda + cs[b[k] - s[k] + 1, 1])/cs[b[k] - s[k] + 1, 2] )
              } else {
                inequalities_list <- c( inequalities_list, (lambda - cs[b[k] - s[k] + 1, 1])/cs[b[k] - s[k] + 1, 2] )
              }
              inequalities_signs <- c( inequalities_signs, ifelse( cs[b[k] - s[k] + 1, 2] > 0, -d[k], d[k] ) )
            }

            ### |C(b[k])| > +/- C(t) for t \neq b[k]
            inequalities <- cbind((cs[,1] - cs[b[k]-s[k]+1, 1])/(cs[b[k]-s[k]+1, 2] - cs[,2]), (-1)*(cs[,1] + cs[b[k]-s[k]+1, 1])/(cs[b[k]-s[k]+1, 2] + cs[,2]))
            signs <- cbind( ifelse( cs[b[k]-s[k]+1, 2] - cs[,2] > 0, -d[k], d[k] ), ifelse( cs[b[k]-s[k]+1, 2] + cs[,2] > 0, -d[k], d[k] ) )

            ### Remove b[k] since we don't need this (will probably return NaN or Inf otherwise)
            inequalities[ b[k]-s[k]+1, ] <- NA

            ### If the denominator of any of the inequalities is 0, replace calculated value (probably NaN) with NA
            nas1 <- (1:nrow(inequalities))[ abs(cs[,2] - cs[b[k]-s[k]+1, 2]) <= 10^(-10) ]
            nas2 <- (1:nrow(inequalities))[ abs(cs[,2] + cs[b[k]-s[k]+1, 2]) <= 10^(-10) ]
            inequalities[ nas1, 1 ] <- NA
            inequalities[ nas2, 2 ] <- NA
            inequalities_list <- c(inequalities_list, inequalities[ !is.na(inequalities) ])
            inequalities_signs <- c(inequalities_signs, signs[ !is.na(inequalities) ])

          ### Otherwise, we need |C(b)| > +/- C(t) for all t
          } else {
            s1 <- interval[1]
            e1 <- interval[2]
            cs2 <- cusum_phi_vec(y, nu, nu2, nuTy, s=s1, e=e1 )
            inequalities <- cbind((cs2[,1] - cs[b[k]-s[k]+1,1])/(cs[b[k]-s[k]+1,2] - cs2[,2]), (-1)*(cs2[,1] + cs[b[k]-s[k]+1,1])/(cs[b[k]-s[k]+1,2] + cs2[,2]))
            signs <- cbind( (cs[b[k]-s[k]+1,2] - cs2[,2]>0), (cs[b[k]-s[k]+1,2] + cs2[,2]>0) ) ## 1 if denominator > 0, 0 if < 0
            signs[ signs==0 ] <- -1 ## replace 0 with -1
            signs <- (-d[k]) * signs ## multiply by -d[1]

            ### If the denominator of any of the inequalities is 0, replace calculated value (probably NaN) with NA
            nas1 <- (1:nrow(inequalities))[ abs(cs2[,2] - cs[b[k] - s[k] + 1, 2]) <= 10^(-10) ]
            nas2 <- (1:nrow(inequalities))[ abs(cs2[,2] + cs[b[k] - s[k] + 1, 2]) <= 10^(-10) ]
            inequalities[nas1, 1] <- NA
            inequalities[nas2, 2] <- NA

            inequalities_list <- c(inequalities_list, inequalities[ !is.na(inequalities) ])
            inequalities_signs <- c(inequalities_signs, signs[ !is.na(inequalities) ])
          }
        }
      }
    }

  }

  ### Find which values of phi satisfy all inequalities
  max_lower_bound <- max(inequalities_list[ inequalities_signs==1 ], na.rm=TRUE)
  min_upper_bound <- min(inequalities_list[ inequalities_signs==-1 ], na.rm=TRUE)
  
  return( c(max_lower_bound, min_upper_bound) )
}