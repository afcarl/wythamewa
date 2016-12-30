###
# functions for fitting the GRETI models

# define a model and its related structures
ewa_def <- function(model,prior,data,file,pars,start,link=c("logit","logit","log","logit","logit"),N_effects,N_veffects,...) {
    if ( missing(file) ) {
        # use default template
        utils::data(template_raw)
        model_template <- template_raw
    } else {
        template_file <- scan(file,character(),sep="\n")
        model_template <- paste(template_file,collapse="\n")
    }
    
    #########################
    # define linear models
    # model should be a list with
    # (1) social learning linear model
    # (2) gamma updating linear model
    # (3) conformity exponent linear model
    # (4) payoff bias convex weight linear model
    # (5) payoff bias Xweight linear model

    if (missing(model)) stop("Need a model.")
    if (missing(prior)) stop("Need to specify prior.")
    if ( class(model)=="character" ) model <- as.list(model)

    # process inverse link functions
    invlink <- list(logit="inv_logit",log="exp",identity="")

    lm <- list("0","0","0","0","0")
    if ( missing(N_effects) ) N_effects <- length(model)
    if ( missing(N_veffects) ) N_veffects <- N_effects
    for ( i in 1:length(model) ) {
        inv_link <- invlink[[link[i]]]
        if ( model[[i]]=="" ) inv_link <- "" # null linear model for this component
        lm[[i]] <- concat(inv_link,"(",model[[i]],")")
    }

    ###################
    # prep data and inits
    if ( missing(data) ) {
        utils::data(WythamUnequal)
        dat <- WythamUnequal
        dat$N_effects <- N_effects
        dat$N_veffects <- N_veffects
        data <- dat
    } else {
        dat <- data
        dat$N_effects <- N_effects
        dat$N_veffects <- N_veffects
        data <- dat
    }

    if (missing(start)) {
        start <- list(
            mu = rep(0,dat$N_effects),
            sigma_bird = rep(0.5,dat$N_veffects),
            z_bird = matrix(0,ncol=dat$N_birds,nrow=dat$N_veffects),
            L_Rho_bird = diag(dat$N_veffects),
            diff_hi = 1,
            b_age = rep(0,dat$N_effects)
        )
    }

    # build per-individual versions of linear models, for generated quantities
    # need to replace bird[i] with individual index values j
    # need to replace age[i] with unique age for individual j in agej[j]
    # we construct agej here and include in the model
    lmj <- lm
    for ( i in 1:length(lmj) ) {
        lmj[[i]] <- gsub("bird[i]","j",lmj[[i]],fixed=TRUE)
        lmj[[i]] <- gsub("age[i]","agej[j]",lmj[[i]],fixed=TRUE)
    }
    # build agej variable
    n_birds <- data$N_birds
    agej <- rep(0,n_birds)
    for ( j in 1:n_birds ) {
        ages_found <- unique(data$age[data$bird==j])
        if ( length(ages_found)>1 ) message(concat("Warning: More than one unique age value found for individual ",j,". Using only first for computing generated quantities."))
        agej[j] <- ages_found[1]
    }
    data$agej <- agej

    # substitutions
    # linear models and priors
    x <- gsub("!s_i!",lm[[1]],model_template,fixed=TRUE)
    x <- gsub("!g_i!",lm[[2]],x,fixed=TRUE)
    x <- gsub("!lambda_i!",lm[[3]],x,fixed=TRUE)
    x <- gsub("!pay_i!",lm[[4]],x,fixed=TRUE)
    x <- gsub("!payX_i!",lm[[5]],x,fixed=TRUE)
    x <- gsub("!priors!",prior,x,fixed=TRUE)
    # generated quantities linear models
    x <- gsub("!s_j!",lmj[[1]],x,fixed=TRUE)
    x <- gsub("!g_j!",lmj[[2]],x,fixed=TRUE)
    x <- gsub("!lambda_j!",lmj[[3]],x,fixed=TRUE)
    x <- gsub("!pay_j!",lmj[[4]],x,fixed=TRUE)
    x <- gsub("!payX_j!",lmj[[5]],x,fixed=TRUE)
    # sew it up
    model_text <- x

    # parameters to sample
    if (missing(pars))
        pars <- c( "mu", "b_age", "diff_hi", "sigma_bird", "Rho_bird", "a_bird", "log_lik", "ks", "kg", "kl", "kp" , "kpx" )

    # build result suitable to pass to ewa_fit
    result <- list(
            stan_code = model_text,
            data = data,
            start = start,
            pars = pars
        )

    return( result )
}

ewa_fit <- function(model,chains=3,cores=3,seed=1,
    control=list( adapt_delta=0.8 , max_treedepth=12 ), 
    iter=2200,warmup=200,refresh=100,...) {

    message("Precompiling Stan model...")
    # precompile model
    m_prep <- stan( 
        model_code=model$stan_code , 
        data=model$data , 
        init=list(model$start) , 
        chains=1 , warmup=1 , iter=2 )

    message("Starting up sampling chains...")
    # hand off to mclapply to simultanesouly sample from 3 chains on 3 cores
    n_chains <- chains
    n_cores <- min( parallel::detectCores() , cores )
    seed <- seed
    require(parallel)
    ctrl <- control
    init_list <- list()
    for ( i in 1:n_chains ) init_list[[i]] <- model$start

    # fit
    return( 
        stan( fit=m_prep , data=model$data , init=init_list , iter=iter , warmup=warmup , chains=n_chains , cores=n_cores , pars=model$pars , control=ctrl , refresh=refresh )
        )
}

# test code
if ( FALSE ) {

# for reading in model tempate and writing out as rda, in order to include in package
# load model template
# collapse into single string
template_raw <- scan("model_EWA_NC_GRETI_templateW.stan",character(),sep="\n")
template_raw <- paste(template_raw,collapse="\n")
save(template_raw,file="template_raw.rda")

library(wythamewa)

# Sva_Gva_Lva_Pva_PXva
lms <- list(
    "mu[1] + a_bird[bird[i],1] + b_age[1]*age[i]",
    "mu[2] + a_bird[bird[i],2] + b_age[2]*age[i]",
    "mu[3] + a_bird[bird[i],3] + b_age[3]*age[i]",
    "mu[4] + a_bird[bird[i],4] + b_age[4]*age[i]",
    "mu[5] + a_bird[bird[i],5] + b_age[5]*age[i]"
    )
prior <- "
    mu ~ normal(0,1);
    diff_hi ~ cauchy(0,1);
    b_age ~ normal(0,1);
    to_vector(z_bird) ~ normal(0,1);
    L_Rho_bird ~ lkj_corr_cholesky(3);
    for ( i in 1:3 ) sigma_bird[i] ~ exponential(1);
    sigma_bird[4] ~ exponential(2);
    sigma_bird[5] ~ exponential(1);"
mod1 <- ewa_def( model=lms , prior=prior )

m <- ewa_fit( mod1 , warmup=50 , iter=100 , chains=2 )

WAIC(m)

post <- extract.samples(m)

ks <- apply(post$ks,2,median)
kg <- apply(post$kg,2,median)
kl <- apply(post$kl,2,median)
kp <- apply(post$kp,2,median)
kpx <- apply(post$kpx,2,median)

ks_sd <- apply(post$ks,2,sd)
kg_sd <- apply(post$kg,2,sd)
kl_sd <- apply(post$kl,2,sd)
kp_sd <- apply(post$kp,2,sd)
kpx_sd <- apply(post$kpx,2,sd)

cav <- Vectorize( col.alpha , "alpha" )

blank(ex=2)
par(mfrow=c(2,2))
a <- cav("orange",1-kg_sd/max(kg_sd))
plot(ks,kg,col=a)
a <- cav("orange",1-kl_sd/max(kl_sd))
plot(ks,kl,col=a); abline(h=1,lty=2)
a <- cav("orange",1-kp_sd/max(kp_sd))
plot(ks,kp,col=a)
a <- cav("orange",1-kp_sd/max(kp_sd))
plot(kp,kpx,col=a)

}