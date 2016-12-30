\name{fitdiff}
\alias{fitdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute invasion fitness difference of learning strategy}
\description{
  Used with \code{\link{sim_birds}} to calculate selection gradients on social learning strategies.
}
\usage{
fitdiff(simdat,simdat2,mutant=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simdat}{Result of \code{\link{sim_birds}} that includes a mutant individual}
  \item{simdat2}{Result of \code{\link{sim_birds}} that includes only common-type strategies}
  \item{mutant}{Position of mutant in \code{simdat}}
}
\details{
  This function simulates groups of birds learning under the assumed statistical model. It produces data suitable to feed directly into the statistical model, for purpose of validating the analysis code. It can also be used to explore a variety of learning dynamics.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{}
\examples{

# calculate gradient over s and conf
s_seq <- c( 0 , 0.05 , 0.1 , 0.15 , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 )
l_seq <- c( 1 , 1.25 , 1.5 , 1.75 , 2.0 , 2.5 , 3.0 )
seq <- expand.grid( s=s_seq , l=l_seq , dl=NA , dlse=NA , ds=NA , dsse=NA )

N_bird <- 10
tmax <- 100
dhi_set <- rep(c(2.5,1),each=tmax/2)
dlo_set <- rep(c(1,2.5),each=tmax/2)

# function to compute common-type fitness
fsx0 <- function(ks,kl) sim_birds( N_bird , tmax , 
        s=rep( ks[1] , N_bird ) , 
        gamma=rep( 0.6 ,N_bird ) , 
        conf=rep( kl[1] , N_bird ) , 
        pay=rep(0,N_bird) , 
        d1=dhi_set , d2=dlo_set )
# function to compute mutant fitness in lambda* (kl)
fsx1 <- function(ks,kl) sim_birds( N_bird , tmax , 
        s=c( ks[1] , rep( ks[1] ,N_bird-1) ) , 
        gamma=c( 0.6 , rep( 0.6 ,N_bird-1) ) , 
        conf=c( kl[2] , rep( kl[1] ,N_bird-1) ) , 
        pay=c( 0 , rep(0,N_bird-1) ) , 
        d1=dhi_set , d2=dlo_set )
# function to compute mutant fitness in s* (ks)
fsx2 <- function(ks,kl) sim_birds( N_bird , tmax , 
        s=c( ks[2] , rep( ks[1] ,N_bird-1) ) , 
        gamma=c( 0.6 , rep( 0.6 ,N_bird-1) ) , 
        conf=c( kl[1] , rep( kl[1] ,N_bird-1) ) , 
        pay=c( 0 , rep(0,N_bird-1) ) , 
        d1=dhi_set , d2=dlo_set )
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

# my multi-core replicate wrapper
mcrepsilent <- function (n, expr, mc.cores = 2) 
{
    result <- simplify2array(mclapply(1:n, eval.parent(substitute(function(i, 
        ...) {
        expr
    })), mc.cores = mc.cores))
    result
}

# this is the top-level function to compute fitness differential
frep <- function(ks,kl,i) {
    n_rep <- 1e4
    # lambda gradient
    xl <- mcrepsilent( n_rep , fitdiff( fsx1(ks,kl) , fsx0(ks,kl) ) , mc.cores=63 )
    # s gradient
    xs <- mcrepsilent( n_rep , fitdiff( fsx2(ks,kl) , fsx0(ks,kl) ) , mc.cores=63 )
    print( concat("[",i,"] ",ks[1]," ",kl[1]) )
    return( c( mean(xl) , sd(xl)/sqrt(n_rep) , mean(xs) , sd(xs)/sqrt(n_rep) ) )
}

# DO NOT RUN THIS, UNLESS YOU HAVE CLUSTER OR A LOT OF TIME
# we used a 64 core cluster and it took 8 or 9 hours
# instead use data(wythamgradient) to load the resulting object into namespace
if ( FALSE ) {
time_start <- Sys.time()
for ( i in 1:nrow(seq) ) {
    seq[i,3:6] <- frep( c(seq[i,1],seq[i,1]+0.1) , c(seq[i,2],seq[i,2]+0.5) , i )
}
time_end <- Sys.time()
( time_end - time_start ) # report elapsed time

}

# write.csv( seq , file="gradientx.csv" , row.names=FALSE )
# can also use data(wythamgradient) to load this result into namespace
data(wythamgradient)
seq <- wythamgradient

#####################
# FIGURE 6
# make contour plot
blank(w=3)
par(mfrow=c(1,3))

n <- length(unique(seq$s))
m <- matrix( seq$dl , nrow = n )
contour( unique(seq$s) , unique(seq$l) , m , xlab="social learning weight" , ylab="conformity strength" )
mtext( "dl" )
lev0 <- contourLines(unique(seq$s) , unique(seq$l) , m , levels=0 )
isol <- lev0
for( i in 1:length(lev0) ) lines( lev0[[i]] , col="red" , lwd=2 )

m <- matrix( seq$ds , nrow = n )
contour( unique(seq$s) , unique(seq$l) , m , xlab="social learning weight" , ylab="conformity strength" )
mtext("ds")
lev0 <- contourLines(unique(seq$s) , unique(seq$l) , m , levels=0 )
isos <- lev0
for( i in 1:length(lev0) ) lines( lev0[[i]] , col="blue" , lwd=2 )

# make vector field plot
plot(NULL, type = "n", xlim=range(seq$s) , ylim=range(seq$l) , xlab="social learning" , ylab="conformity strength" )
fz <- function(x) sign(x)*abs(x)^(0.5)/20
arrows(seq[,1], seq[,2], seq[,1] + fz(seq[,5]), seq[,2] + fz(seq[,3]) , length=0.05 )

# add isoclines
for( i in 1:length(isol) ) lines( isol[[i]] , col="red" , lty=2 )
for( i in 1:length(isos) ) lines( isos[[i]] , col="blue" , lty=2 )


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

