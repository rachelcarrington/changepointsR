calculate_interval_wbs <- function(y, nu, results=NULL, b=NULL, d=NULL, s=NULL, e=NULL, rand_ints=NULL, nu2=NULL, nuTy=NULL, 
                                     threshold=NULL, n.cp=NULL ){
  
  ##### Find values of phi which satisfy the required inequalities for results of WBS/SBS
  ##### Option to only consider part of (b, d), e.g. if we just want b1 & d1 to be specified, & don't care about
  ##### later values
  
  ### Note that if cs & nuTy are supplied y is ignored
  
  ### If stop.after is supplied, the function will return the interval in which the first stop.after changepoints 
  ### and signs found match those in b and d. If threshold is also supplied, this will also be used.
  ### Otherwise, if threshold is supplied (but not n.cp), the interval returned will be that for which the complete 
  ### b and d match those given, when BS is given this threshold. If threshold is not supplied, it will be assumed 
  ### that BS is run with a fixed number of iterations. If both n.cp and threshold are supplied, it will be assumed that
  ### n.cp is the maximum number of iterations, but the threshold will also be used as a minimum for C(t).

  if ( is.null(results) ){
    if ( is.null(b) || is.null(d) || is.null(s) || is.null(e) || is.null(rand_ints) ){
      stop("Not enough parameters given.")
    }
  } else {
    b <- results$results$b
    d <- results$results$d
    s <- results$results$s
    e <- results$results$e
    rand_ints <- results$rand_ints
    threshold <- results$threshold
  }

  
  if ( is.null( n.cp ) ){
    if ( is.null( threshold ) ){
      n.cp <- length(b) ## assume num. of iterations was specified since threshold is not given
    } else {
      n.cp <- length(y) ## just use threshold
    }
  }
  
  
  if ( is.null(nu2) ){
    nu2 <- sum(nu^2)
  }
  
  if ( is.null(nuTy) ){
    nuTy <- as.numeric(t(nu)%*%y)
  }
  
  n <- length(y)
  
  ### If no CPs (threshold must be specified):
  if ( sum(!is.na(b))==0 ){

    max_lower_bound <- -Inf
    min_upper_bound <- Inf

    # For all intervals (s,e), we need |C_{s,e}(t)| < \lambda for all t

    for ( ind in 1:nrow(rand_ints) ){

      cs <- cusum_phi_vec(y, nu, nu2, nuTy, s=rand_ints[ind,1], e=rand_ints[ind,2])
      x <- abs(cs[,2]) > 10^(-10) ## if cs[i,2] = 0, then this inequality is constant in phi
      inequalities <- c( (threshold - cs[x,1]) / cs[x,2], (-threshold - cs[x,1]) / cs[x,2] )
      signs <- c( ifelse( cs[x,2] > 0, -1, 1 ), ifelse( cs[x,2] > 0, 1, -1 ) )
      max_lower_bound <- max( c( max_lower_bound, inequalities[ signs==1] ) )
      min_upper_bound <- min( c( min_upper_bound, inequalities[ signs==-1] ) )
    
    }


  } else {
    
    ### For the first CP:
    cs <- cusum_phi_vec(y, nu, nu2, nuTy, s=s[1], e=e[1])

    #### (-d1) * C(b1) > threshold (0 if not given)
    tau <- b[1] - s[1] + 1
    if ( abs( cs[tau,2] ) > 10^(-10) ){
      if ( is.null(threshold) ){
        inequalities <- (-1) * cs[tau,1] / cs[tau,2]
      } else {
        inequalities <- (-d[1]*threshold - cs[tau,1]) / cs[tau,2]
      }
      signs <- ifelse( cs[tau,2] > 0, -d[1], d[1] )
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
    
    #### | C( b1 ) | > +/- C( t ) for t \neq b1, on interval (s[1], e[1])
    inequalities <- cbind( (cs[-tau,1] - cs[tau,1])/(-cs[-tau,2] + cs[tau,2] ), 
                           (-1)*(cs[-tau,1] + cs[tau,1])/(cs[-tau,2] + cs[tau,2]) )
    signs <- cbind( ifelse( -cs[-tau,2] + cs[tau,2] > 0, -d[1], d[1] ),
                    ifelse( cs[-tau,2] + cs[tau,2] > 0, -d[1], d[1] ) )
    
    #### Remove entries where we have divided by 0
    signs[ abs(-cs[-tau,2] + cs[tau,2]) <= 10^(-10), 1 ] <- NA
    signs[ abs(cs[-tau,2] + cs[tau,2]) <= 10^(-10), 2 ] <- NA
    
    #### Update bounds
    max_lower_bound <- max( c( max_lower_bound, inequalities[ signs==1 ] ), na.rm=TRUE )
    min_upper_bound <- min( c( min_upper_bound, inequalities[ signs==-1 ] ), na.rm=TRUE )


    #### Check other intervals: | C_{s1, e1}( b1 ) | > |C_{s,e} ( t )| for all s, e, t
    for ( ind in 1:nrow(rand_ints) ){

      s2 <- rand_ints[ind,1]
      e2 <- rand_ints[ind,2]
      if ( sum( abs( c(s[1],e[1]) - c(s2,e2) ) ) != 0 ){

        cs2 <- cusum_phi_vec(y, nu, nu2, nuTy, s=s2, e=e2)
        inequalities <- cbind( (cs2[,1] - cs[tau,1])/(-cs2[,2] + cs[tau,2] ), 
                               (-1)*(cs2[,1] + cs[tau,1])/(cs2[,2] + cs[tau,2]) )
        signs <- cbind( ifelse( -cs2[,2] + cs[tau,2] > 0, -d[1], d[1] ),
                        ifelse( cs2[,2] + cs[tau,2] > 0, -d[1], d[1] ) )
    
        #### Remove entries where we have divided by 0
        signs[ abs(-cs2[,2] + cs[tau,2]) <= 10^(-10), 1 ] <- NA
        signs[ abs(cs2[,2] + cs[tau,2]) <= 10^(-10), 2 ] <- NA
    
        #### Update bounds
        max_lower_bound <- max( c( max_lower_bound, inequalities[ signs==1 ] ), na.rm=TRUE )
        min_upper_bound <- min( c( min_upper_bound, inequalities[ signs==-1 ] ), na.rm=TRUE )

      }

    }
    

    ### For later CPs, if any:
    if ( n.cp > 1 & length(b) > 1 ){
      
      for ( k in 2:min( c( n.cp, length(b) ) ) ){

        #### We need -d[k]*C(b[k]) to be greater than 0 or the threshold
        cs <- cusum_phi_vec(y, nu, nu2, nuTy, s=s[k], e=e[k])
        tau <- b[k] - s[k] + 1
        if ( abs( cs[tau,2] ) > 10^(-10) ){
           if ( is.null(threshold) ){
            inequalities <- (-1) * cs[tau,1] / cs[tau,2]
          } else {
            inequalities <- (-d[k]*threshold - cs[tau,1]) / cs[tau,2]
          }
          signs <- ifelse( cs[tau,2] > 0, -d[k], d[k] )
          if ( signs==1 ){
            max_lower_bound <- max( max_lower_bound, inequalities )
          } else {
            min_upper_bound <- min( min_upper_bound, inequalities )
           }
        }
        
        #### Also, we need |C(b[k])| > |C(b[t])| for t \neq k
        inequalities <- cbind( (cs[-tau,1] - cs[tau,1])/(-cs[-tau,2] + cs[tau,2] ), 
                               (-1)*(cs[-tau,1] + cs[tau,1])/(cs[-tau,2] + cs[tau,2]) )
        signs <- cbind( ifelse( -cs[-tau,2] + cs[tau,2] > 0, -d[k], d[k] ),
                        ifelse( cs[-tau,2] + cs[tau,2] > 0, -d[k], d[k] ) )
        
        #### Remove entries where we have divided by 0
        signs[ abs(-cs[-tau,2] + cs[tau,2]) <= 10^(-10), 1 ] <- NA
        signs[ abs(cs[-tau,2] + cs[tau,2]) <= 10^(-10), 2 ] <- NA
        
        #### Remove NAs & recalculate bounds
        max_lower_bound <- max( max_lower_bound, max( inequalities[ signs==1 ] ), na.rm=TRUE )
        min_upper_bound <- min( min_upper_bound, min( inequalities[ signs==-1 ] ), na.rm=TRUE )

        
        #### Check other intervals: | C_{s1, e1}( b1 ) | > |C_{s,e} ( t )| for all s, e, t

        #### Remove intervals which contain the previous changepoint
        rand_ints <- rand_ints[ !( rand_ints[,1] <= b[k-1] & rand_ints[,2] > b[k-1] ),,drop=FALSE ]

        if ( nrow(rand_ints) >= 1 ){        
          for ( ind in 1:nrow(rand_ints) ){

            s2 <- rand_ints[ind,1]
            e2 <- rand_ints[ind,2]
            if ( sum( abs( c(s[1],e[1]) - c(s2,e2) ) ) != 0 ){

              cs2 <- cusum_phi_vec(y, nu, nu2, nuTy, s=s2, e=e2)
              inequalities <- cbind( (cs2[,1] - cs[tau,1])/(-cs2[,2] + cs[tau,2] ), 
                                     (-1)*(cs2[,1] + cs[tau,1])/(cs2[,2] + cs[tau,2]) )
              signs <- cbind( ifelse( -cs2[,2] + cs[tau,2] > 0, -d[k], d[k] ),
                              ifelse( cs2[,2] + cs[tau,2] > 0, -d[k], d[k] ) )
    
              #### Remove entries where we have divided by 0
              signs[ abs(-cs2[,2] + cs[tau,2]) <= 10^(-10), 1 ] <- NA
              signs[ abs(cs2[,2] + cs[tau,2]) <= 10^(-10), 2 ] <- NA
    
              #### Update bounds
              max_lower_bound <- max( c( max_lower_bound, inequalities[ signs==1 ] ), na.rm=TRUE )
              min_upper_bound <- min( c( min_upper_bound, inequalities[ signs==-1 ] ), na.rm=TRUE )

            }

          }
        }        


      }

    }
    
    if ( n.cp > length(b) ){

      ### Check that for the next k, all |Ct|'s are below the threshold
      k <- length(b) + 1

      #### Remove intervals which contain the previous changepoint
      rand_ints <- rand_ints[ !( rand_ints[,1] <= b[k-1] & rand_ints[,2] > b[k-1] ),,drop=FALSE ]
 
      #### Calculate bounds for which all C(t)'s in each RI will be below the threshold
      if ( nrow(rand_ints) >= 1 ){
        for ( ind in 1:nrow(rand_ints) ){

          cs <- cusum_phi_vec(y, nu, nu2, nuTy, s=rand_ints[ind,1], e=rand_ints[ind,2])
          x <- abs(cs[,2]) > 10^(-10) ## if cs[i,2] = 0, then this inequality is constant in phi
          inequalities <- c( (threshold - cs[x,1]) / cs[x,2], (-threshold - cs[x,1]) / cs[x,2] )
          signs <- c( ifelse( cs[x,2] > 0, -1, 1 ), ifelse( cs[x,2] > 0, 1, -1 ) )
          max_lower_bound <- max( c( max_lower_bound, inequalities[ signs==1] ) )
          min_upper_bound <- min( c( min_upper_bound, inequalities[ signs==-1] ) )

        }
      }

    }

  }
  
    
  ### Inequalities are all satisfied by phi s.t. max_lower_bound < phi < min_upper_bound
  
  return( c( max_lower_bound, min_upper_bound ) )
}