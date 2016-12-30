# functions for simulation and model validation

sim_birds <- function(N=10,N_solve=10,s,gamma,conf,pay,d1=2.5,d2=1,A0=c(0,4)) {
    solve <- matrix(0,nrow=N_solve,ncol=N)
    A <- matrix(0,nrow=N,ncol=2) # attractions
    # check lengths of inputs
    if ( length(s)<N ) s <- rep_len(s,N)
    if ( length(pay)<N ) pay <- rep_len(pay,N)
    if ( length(conf)<N ) conf <- rep_len(conf,N)
    if ( length(gamma)<N ) gamma <- rep_len(gamma,N)
    if ( length(d1)<N_solve ) d1 <- rep_len(d1,N_solve)
    if ( length(d2)<N_solve ) d2 <- rep_len(d2,N_solve)
    # init attractions
    for ( i in 1:N ) A[i,] <- A0
    # sim solves
    for ( j in 1:N_solve ) {
        for ( i in 1:N ) {
            PI <- softmax( A[i,1] , A[i,2] )[1]
            if ( j==1 ) {
                # first turn, so no social info
                PrA <- PI
            } else {    
                # count number hi solves from previous turn
                n_hi <- sum( solve[j-1,-i] )
                n_lo <- N - n_hi
                CONFI <- n_hi^conf[i]/(n_hi^conf[i] + n_lo^conf[i])
                if ( d1[j-1] > d2[j-1] ) 
                    do_payoff_bias <- ifelse( n_hi > 0 , 1 , 0.5 )
                if ( d1[j-1] < d2[j-1] ) 
                    do_payoff_bias <- ifelse( n_lo > 0 , 0 , 0.5 )
                SI <- (1-pay[i])*CONFI + pay[i]*do_payoff_bias
                PrA <- (1-s[i])*PI + s[i]*SI
            }
            solve[j,i] <- rbinom(1,size=1,prob=PrA)
            # update attractions
            payoffs <- c(0,0)
            if ( solve[j,i]==1 ) {
                payoffs[1] <- d1[j]
            } else {
                payoffs[2] <- d2[j]
            }
            A[i,] <- (1-gamma[i])*A[i,] + gamma[i]*payoffs
        }#i
    }#j
    # result has individuals in columns and solves in rows
    attr(solve,"d1") <- d1
    attr(solve,"d2") <- d2
    return(solve)
}

# avgpayoff per bird
avgpayoff <- function(simdat) {
    d1 <- attr(simdat,"d1")
    d2 <- attr(simdat,"d2")
    tmax <- nrow(simdat)
    n <- ncol(simdat)
    p <- 0
    for ( i in 1:tmax ) {
        p <- p + (simdat[i,]*d1[i])
        p <- p + ((1-simdat[i,])*d2[i])
    }
    return( p / tmax )
}

# compute fitness difference between individuals
fitdiff <- function(simdat,simdat2,mutant=1) {
    tmax <- nrow(simdat)
    n <- ncol(simdat)
    #non_mutant <- (1:n)[-mutant]
    non_mutant <- 1:n # common-types in simdat2
    d1 <- attr(simdat,"d1")
    d2 <- attr(simdat,"d2")
    p <- 0
    p2 <- 0
    for ( i in 1:tmax ) {
        p <- p + (simdat[i,]*d1[i])
        p <- p + ((1-simdat[i,])*d2[i])
        p2 <- p2 + (simdat2[i,]*d1[i])
        p2 <- p2 + ((1-simdat2[i,])*d2[i])
    }
    Wmutant <- mean(p[mutant])
    Wnmutant <- mean(p2[non_mutant])
    return( Wmutant - Wnmutant )
}
