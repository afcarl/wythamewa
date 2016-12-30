\name{ewa_def}
\alias{ewa_def}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Define a custom EWA dynamic learning model}
\description{
  This function take linear models and priors as inputs and defines a custom EWA dynamic learning model.
}
\usage{
ewa_def(model,prior,data,file,pars,start,
  link=c("logit","logit","log","logit","logit"),
  N_effects,N_veffects,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A \code{list} of linear models for each learning parameter. See examples.}
  \item{prior}{A quoted string of priors to insert into the model definition. See examples.}
  \item{data}{Optional. A \code{list} of data. If missing, uses \code{data(WythamUnequal)}.}
  \item{file}{Optional model template. If missing, uses \code{data(template_raw)}}
  \item{pars}{Optional character vector of parameter names to sample}
  \item{start}{Optional named \code{list} of starting values for parameters}
  \item{link}{Character vector of link function names for each linear model in \code{model}}
  \item{N_effects}{Optional number of linear models. Defaults to \code{length(model)}.}
  \item{N_veffects}{Optional number of varying effects. Defaults to \code{N_effects}.}
}
\details{
  This function uses a model template, by default \code{data(template_raw)}, to define a custom dynamic learning model to fit in Stan. The template contains fields for linear models for each learning parameter. The first four parameters are (1) s, (2) g, (3) lambda, and (4) y, as defined in the paper. The fifth is an alternative parameterization of the payoff bias, called px in the code.

  The expressions in \code{model} define these paramters as combinations of means and varying effects and any other effects desired. The distributions in \code{prior} define priors for each parameter. 
}
\value{
}
\references{}
\author{Richard McElreath}
\seealso{\code{\link{ewa_fit}}}
\examples{

# define full model from paper
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

# display Stan model code
cat(mod1$stan_code)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }

