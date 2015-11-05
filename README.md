wythamewa
==========

This package contains models and data to replicate the payoff and frequency dependent social learning analysis of the Wytham Woods Parus major data with unequal payoff options.

The full model can be replicated with the following code:
```
lms <- list(
    "mu[1] + a_bird[bird[i],1] + b_age[1]*age[i]",
    "mu[2] + a_bird[bird[i],2] + b_age[2]*age[i]",
    "mu[3] + a_bird[bird[i],3] + b_age[3]*age[i]",
    "mu[4] + a_bird[bird[i],4] + b_age[4]*age[i]",
    "mu[5] + a_bird[bird[i],5] + b_age[5]*age[i]"
    )
links <- c("logit", "logit", "log", "logit", "logit")
prior <- "
    mu ~ normal(0,1);
    diff_hi ~ cauchy(0,1);
    b_age ~ normal(0,1);
    to_vector(z_bird) ~ normal(0,1);
    L_Rho_bird ~ lkj_corr_cholesky(3);
    sigma_bird ~ exponential(2);"
library(wythamewa)
mod1 <- ewa_def( model=lms , prior=prior , link=links )
m <- ewa_fit( mod1 , warmup=1000 , iter=2000 , chains=4 , cores=4 ,
    control=list( adapt_delta=0.99 , max_treedepth=12 ) )
```
