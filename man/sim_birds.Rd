\name{sim_birds}
\alias{sim_birds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulated social learning data}
\description{
  Flexible simulation of social and individual learning in foraging groups.
}
\usage{
sim_birds(N=10,N_solve=10,s,gamma,conf,pay,d1=2.5,d2=1,A0=c(0,4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Number of individuals}
  \item{N_solve}{Number of time steps}
  \item{s}{Vector of social learning weight values}
  \item{gamma}{Vector of updating rate values}
  \item{conf}{Vector of conformity strength values}
  \item{pay}{Vector of payoff bias weight values}
  \item{d1}{Payoff to first behavioral option, or vector of payoffs}
  \item{d2}{Payoff to second option, or vector of payoffs}
  \item{A0}{Initial attraction values}
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

#####################################
# examples of model behavior
# this code reproduces Figure 5.

dhi_set <- rep(c(2.5,1),each=30)
dlo_set <- rep(c(1,2.5),each=30)

# if conformity too strong, population cannot track high payoff
set.seed(100)
sim_dat <- sim_birds( 10 , 60 , s=0.5 , gamma=0.6 , conf=10 , pay=0 , 
    d1=dhi_set , d2=dlo_set )
plot( NULL , xlim=c(1,nrow(sim_dat)) , ylim=c(0,ncol(sim_dat)+1) , xlab="turn" , ylab="bird" , yaxt="n" )
axis( 2 , at=1:ncol(sim_dat) , labels=1:ncol(sim_dat) )
for ( i in 1:ncol(sim_dat) ) 
    points( 1:nrow(sim_dat) , rep(i,nrow(sim_dat)) , pch=ifelse(sim_dat[,i]==1,16,1) )
abline(v=30.5,lty=2)
mtext( "Too much conformity" )

# moderate conformity, population can track high payoff
set.seed(101)
sim_dat <- sim_birds( 10 , 60 , s=0.4 , gamma=0.6 , conf=5 , pay=0 , 
    dhi=rep(c(2.5,1),each=30) , dlo=rep(c(1,2.5),each=30) )
plot( NULL , xlim=c(1,nrow(sim_dat)) , ylim=c(0,ncol(sim_dat)+1) , xlab="turn" , ylab="bird" , yaxt="n" )
axis( 2 , at=1:ncol(sim_dat) , labels=1:ncol(sim_dat) )
for ( i in 1:ncol(sim_dat) ) 
    points( 1:nrow(sim_dat) , rep(i,nrow(sim_dat)) , pch=ifelse(sim_dat[,i]==1,16,1) )
abline(v=30.5,lty=2)
mtext( "Enough conformity" )

# Too little conformity, slower tracking
set.seed(101)
sim_dat <- sim_birds( 10 , 60 , s=0.4 , gamma=0.6 , conf=1 , pay=0 , 
    dhi=rep(c(2.5,1),each=30) , dlo=rep(c(1,2.5),each=30) )
plot( NULL , xlim=c(1,nrow(sim_dat)) , ylim=c(0,ncol(sim_dat)+1) , xlab="turn" , ylab="bird" , yaxt="n" )
axis( 2 , at=1:ncol(sim_dat) , labels=1:ncol(sim_dat) )
for ( i in 1:ncol(sim_dat) ) 
    points( 1:nrow(sim_dat) , rep(i,nrow(sim_dat)) , pch=ifelse(sim_dat[,i]==1,16,1) )
abline(v=30.5,lty=2)
mtext( "Too little conformity" )

# now sample from actual posterior and simulate
# need posterior samples in object "post"
# see ?ewa_fit to produce it
data(WythamUnequal)
dat <- WythamUnequal
N_birds <- dat$N_birds
ks <- apply(post$ks,2,mean)
kg <- apply(post$kg,2,mean)
kl <- apply(post$kl,2,mean)
kp <- apply(post$kp,2,mean)

N <- 10
idx <- sample( 1:N_birds , size=N )
N_solve <- 60

sim_dat <- sim_birds( N , N_solve , s=ks[idx] , gamma=kg[idx] , conf=kl[idx] , pay=kp[idx] , 
    d1=rep(c(2.5,1),each=N_solve/2) , d2=rep(c(1,2.5),each=N_solve/2) )

plot( NULL , xlim=c(1,nrow(sim_dat)) , ylim=c(0,ncol(sim_dat)+1) , xlab="turn" , ylab="bird" , yaxt="n" )
axis( 2 , at=1:ncol(sim_dat) , labels=1:ncol(sim_dat) )
for ( i in 1:ncol(sim_dat) ) 
    points( 1:nrow(sim_dat) , rep(i,nrow(sim_dat)) , pch=ifelse(sim_dat[,i]==1,16,1) )
abline(v=N_solve/2 + 0.5 , lty=2 , col="red" )
mtext("Individuals sampled from posterior distribution")


###########################
# fake data for validating model fitting code

N <- 20
N_solve <- 400
s <- logistic(rnorm(N))
pay <- logistic(rnorm(N,-0.5,0.1))
conf <- 1.5 + runif(N,0,0.1)
sim_dat <- sim_birds(N=N,N_solve=N_solve,s=s,pay=pay,conf=conf,gamma=s/2)
x <- sim_dat

# put data into format that Stan code assumes
nrows <- N*N_solve
bird_id <- rep( 1:N , each=N_solve )
solve_hi <- as.vector(x)
turn <- rep( 1:N_solve , times=N )
s_prop <- sapply( 1:nrows , function(i) {
        if (turn[i]>1) {
            return( (sum(x[turn[i],]) - x[turn[i],bird_id[i]])/(N-1) )
        } else {
            return(-1)
        }
    } )
A_init <- rep( 4 , times=N )

# fit model
# prep data to pass to Stan
datx <- list(
    N = nrows,
    N_birds = N,
    bird = bird_id,
    solve_hi = solve_hi,
    s_prop = s_prop,
    age = rep(0,nrows),
    A_init = A_init
)
# define full model
lms <- list(
    "mu[1] + a_bird[bird[i],1] + b_age[1]*age[i]",
    "mu[2] + a_bird[bird[i],2] + b_age[2]*age[i]",
    "mu[3] + a_bird[bird[i],3] + b_age[3]*age[i]",
    "mu[4] + a_bird[bird[i],4] + b_age[4]*age[i]"
    )
links <- c("logit", "logit", "log", "logit", "")
prior <- "
    mu ~ normal(0,1);
    diff_hi ~ cauchy(0,1);
    b_age ~ normal(0,1);
    to_vector(z_bird) ~ normal(0,1);
    L_Rho_bird ~ lkj_corr_cholesky(3);
    sigma_bird ~ exponential(2);"
mod1 <- ewa_def( model=lms , prior=prior , link=links , data=datx )
# run the model
m <- ewa_fit( mod1 , warmup=500 , iter=1000 , chains=3 , cores=3 ,
    control=list( adapt_delta=0.99 , max_treedepth=12 ) )

# check convergence
precis( m , pars=c("mu","b_age","sigma_bird") , depth=2 )

# process and plot
post <- extract.samples(m)
ks <- apply(post$ks,2,mean)
kg <- apply(post$kg,2,mean)
kl <- apply(post$kl,2,mean)
kp <- apply(post$kp,2,mean)
ks_sd <- apply(post$ks,2,sd)
kg_sd <- apply(post$kg,2,sd)
kl_sd <- apply(post$kl,2,sd)
kp_sd <- apply(post$kp,2,sd)

# validation plots
the_col <- rangi2
blank(ex=1.66)
par(mfrow=c(2,2))

plot( s , ks , xlab="social learning weight (true)" , ylab="social learning weight (estimated)" , col=rangi2 , xlim=c(0,1) , ylim=c(0,1) )
for ( i in 1:length(ks) ) lines( c(s[i],s[i]) , ks[i]+c(1,-1)*ks_sd[i] , col=rangi2 )
abline(a=0,b=1,lty=2)

plot( s/2 , kg , xlab="updating rate (true)" , ylab="updating rate (estimated)" , col=rangi2 )
for ( i in 1:length(kg) ) lines( c(s[i]/2,s[i]/2) , kg[i]+c(1,-1)*kg_sd[i] , col=rangi2 )
abline(a=0,b=1,lty=2)

if ( length(conf)==1 ) conf <- rep(conf,length(kl)) + runif(length(kl),-0.1,0.1)
plot( conf , kl , xlab="conformity strength (true)" , ylab="conformity strength (estimated)" , col=rangi2 , xlim=c(1,3) , ylim=c(1,5) )
for ( i in 1:length(kl) ) lines( c(conf[i],conf[i]) , kl[i]+c(1,-1)*kl_sd[i] , col=rangi2 )
abline(a=0,b=1,lty=2)

if ( length(pay)==1 ) pay <- rep(pay,length(kp)) + runif(length(kp),-0.1,0.1)
plot( pay , kp , xlab="payoff bias weight (true)" , ylab="payoff bias weight (estimated)" , col=rangi2 )
for ( i in 1:length(kp) ) lines( c(pay[i],pay[i]) , kp[i]+c(1,-1)*kp_sd[i] , col=rangi2 )
abline(a=0,b=1,lty=2)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

