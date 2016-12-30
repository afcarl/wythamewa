\name{ewa_fit}
\alias{ewa_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a custom EWA dynamic learning model}
\description{
  This function takes a defined EWA dynamic learning model and passes it to Stan.
}
\usage{
ewa_fit(model,chains=3,cores=3,seed=1,
    control=list( adapt_delta=0.8 , max_treedepth=12 ), 
    iter=2200,warmup=200,refresh=100,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A list of objects as defined by \code{\link{ewa_def}}}
  \item{chains}{Number of chains to run}
  \item{cores}{Number of cores to run chains on}
  \item{seed}{Unused. Use \code{set.seed} before call to set seed.}
  \item{control}{Control parameters for Stan samplers}
  \item{iter}{Number of total iterations for each chain}
  \item{warmup}{Number of warmup samples for each chain}
  \item{refresh}{Interval to refresh sampling progress display}
}
\details{
  This function uses a list defined by \code{\link{ewa_def}} to sample from an EWA dynamic learning model. See examples below for code to replicate model fit from paper, as well as diagnostics, Table 1, and Figures 3 and 4.
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{ewa_fit}}}
\examples{

#######################################################
# This code defines the model with:
# (1) individual random effects and age effects for social learning strength (s_i)
# (2) individual random effects and age effects for updating rate (g_i)
# (3) individual random effects and age effects for conformist exponent (lambda_i)
# (4) individual random effects and age effects for payoff bias reliance (y_i, as defined in paper)

model_name <- "Sva_Gva_Lva_Pva"
# linear models:
# these are merely sums of mean of each parameter, zero-centered random effect for each bird, and then age effect (age is standardized, so also zero-centered)
# order is: (1) s, (2) g, (3) lambda, (4) y (kp internally)
lms <- list(
    "mu[1] + a_bird[bird[i],1] + b_age[1]*age[i]",
    "mu[2] + a_bird[bird[i],2] + b_age[2]*age[i]",
    "mu[3] + a_bird[bird[i],3] + b_age[3]*age[i]",
    "mu[4] + a_bird[bird[i],4] + b_age[4]*age[i]"
    )
links <- c("logit", "logit", "log", "logit", "")
# priors are weakly informative and regularizing
# note the non-centered parameterization of random effects, using the z_bird vector of standardized Gaussian priors
# for explanation, see pages 408-409 of the Statistical Rethinking textbook.
prior <- "
    mu ~ normal(0,1);
    diff_hi ~ cauchy(0,1);
    b_age ~ normal(0,1);
    to_vector(z_bird) ~ normal(0,1);
    L_Rho_bird ~ lkj_corr_cholesky(3);
    sigma_bird ~ exponential(2);"

# define and run the model
mod1 <- ewa_def( model=lms , prior=prior , link=links )
set.seed(1)
m <- ewa_fit( mod1 , warmup=500 , iter=1000 , chains=3 , cores=3 ,
    control=list( adapt_delta=0.99 , max_treedepth=12 ) )

# diagnostics
# this code extracts and sorts Rhat values for parameters
# all should be close to 1, or less than 1.1 in any event
# note that parameter y from paper is called "kp" in code

pars <- mod1$pars
xpars <- pars[ -which(pars=="log_lik") ]
xpars <- xpars[ -which(xpars=="ks") ]
xpars <- xpars[ -which(xpars=="kg") ]
xpars <- xpars[ -which(xpars=="kl") ]
xpars <- xpars[ -which(xpars=="kp") ]
xpars <- xpars[ -which(xpars=="kpx") ]
xpars <- xpars[ -which(pars=="log_lik") ]
xpars <- xpars[ -which(xpars=="a_bird") ]
x <- summary(m,pars=xpars)$summary
m_diag <- x[, c(1:3,9:10) ]
sort( m_diag[,5] )

# traceplots
tracerplot( m , pars=c("mu","b_age","sigma_bird") )

# summarize mean learning effects
( pt <- precis( m , pars=c("mu","b_age","sigma_bird") , depth=2 ) )
# convert to LaTeX table, to replicate Table 1 in paper
library(xtable); xtable(pt@output)

# process estimates
post <- extract.samples(m)

# prep data for plotting
ns <- dim(post$mu)[1]
nj <- dim(post$a_bird)[2]
data(WythamUnequal)
dat <- WythamUnequal
agej <- sapply( 1:dat$N_birds , function(id) unique(dat$age[dat$bird==id]) )

# extract posterior means of individual learning parameters
ks <- apply(post$ks,2,mean)
kg <- apply(post$kg,2,mean)
kl <- apply(post$kl,2,mean)
ky <- apply(post$kp,2,mean)

##################
# this code replicates Figure 3
blank(ex=1.66)
par(mfrow=c(2,2))
bcol <- rangi2 # point color

plot(ks,kg,col=bcol , xlab="social learning weight (s)" , ylab="updating rate (g)" )
text( 0.03 , 0.87 , "(a)" , cex=1.3 )

plot(ks,kl,col=bcol, xlab="social learning weight (s)" , ylab="conformity exponent (lambda)" ); 
abline(h=1,lty=2)
text( 0.67 , 4.15 , "(b)" , cex=1.3 )

plot(ks,ky,col=bcol, xlab="social learning weight (s)" , ylab="payoff bias weight (y)" )
text( 0.67 , 0.08 , "(c)" , cex=1.3 )

# sigmoid learning curves
flearn <- function(x,i) {
    # x: proportion observed
    # i: individual id
    (1-ky[i])*x^kl[i]/(x^kl[i]+(1-x)^kl[i]) + ky[i]*ifelse(x>0,1,0)
}
plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="observed freq high option" , ylab="probability of high option" )
abline(a=0,b=1,lty=2)
fseq <- seq(from=0,to=1,length.out=50)
lcol <- col.alpha(bcol,0.15)
for ( i in 1:length(ks) ) {
    f <- sapply( fseq , function(x) flearn(x,i) )
    lines( fseq , f , col=lcol )
}
text( 0.05 , 0.95 , "(d)" , cex=1.3 )

##########
# sigmoid curves from a single individual, showing uncertainty
# not in paper
id <- 20
ns <- 100 # number or curves to plot
plot( NULL , xlim=c(0,1) , ylim=c(0,1) , xlab="observed freq high option" , ylab="probability of high option" )
abline(a=0,b=1,lty=2)
fseq <- seq(from=0,to=1,length.out=50)
lcol <- col.alpha(bcol,0.15)
for ( i in 1:ns ) {
    f <- sapply( fseq , function(x) (1-post$kp[i,id])*x^post$kl[i,id]/(x^post$kl[i,id]+(1-x)^post$kl[i,id]) + post$kp[i,id]*ifelse(x>0,1,0) )
    lines( fseq , f , col=lcol )
}
f <- sapply( fseq , function(x) flearn(x,id) )
lines( fseq , f , col="black" )
mtext( concat("Individual ",id) )

###################
# plot age affects - replicates Figure 4

# extract actual bird ages
data(WythamUnequal)
dat <- WythamUnequal
age_id <- sapply( 1:dat$N_birds , function(id) unique(dat$age[dat$bird==id]) )
age_idj <- jitter(age_id) # jitter for easy vizualization
age_r <- range(dat$age)
age_seq <- seq(from=age_r[1]-0.5,to=age_r[2]+0.5,length.out=30)
Ns <- dim(post$mu)[1]

blank(w=2)
par(mfrow=c(1,2))
a <- rangi2 # point color
kp <- ky # renaming so cosistent with old code convention

##### s
plot( age_idj , ks , xlab="age (standardized)" , ylab="social learning weight" , col=a )
socV <- sapply( age_seq , function(age) 
        logistic(post$mu[,1] + post$b_age[,1]*age + 0*rnorm(Ns,0,post$sigma_bird[,1]) ) )
soc_med <- apply( socV ,2,median)
soc_PI <- apply( socV ,2,PI)
lines( age_seq , soc_med )
shade( soc_PI , age_seq )

##### g
plot( age_idj , kg , xlab="age (standardized)" , ylab="updating rate" , col=a )
lV <- sapply( age_seq , function(age) 
        logistic(post$mu[,2] + post$b_age[,2]*age + 0*rnorm(Ns,0,post$sigma_bird[,2]) ) )
l_med <- apply( lV ,2,median)
l_PI <- apply( lV ,2,PI)
lines( age_seq , l_med )
shade( l_PI , age_seq )

##### lambda - not in paper

plot( age_idj , kl , xlab="age (std)" , ylab="conformity strength" , col=a )
abline( h=1 , lty=2 )
lV <- sapply( age_seq , function(age) 
        exp(post$mu[,3] + post$b_age[,3]*age + 0*rnorm(Ns,0,post$sigma_bird[,3]) ) )
l_med <- apply( lV ,2,median)
l_PI <- apply( lV ,2,PI)
lines( age_seq , l_med )
# shade( l_PI , age_seq )
lines( age_seq , l_PI[1,] , col=col.alpha("black",0.5) )
lines( age_seq , l_PI[2,] , col=col.alpha("black",0.5) )

##### y - not in paper
plot( age_idj , kp , xlab="age (standardized)" , ylab="payoff bias weight" , col=a ,
    ylim=c(0,0.1) )
lV <- sapply( age_seq , function(age) 
        logistic(post$mu[,4] + post$b_age[,4]*age + 0*rnorm(Ns,0,post$sigma_bird[,4]) ) )
l_med <- apply( lV ,2,median)
l_PI <- apply( lV ,2,PI)
lines( age_seq , l_med )
# shade( l_PI , age_seq )
lines( age_seq , l_PI[1,] , col=col.alpha("black",0.5) )
lines( age_seq , l_PI[2,] , col=col.alpha("black",0.5) )


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

