
library(testthat)
library(rstan)
library(Deriv)
library(truncnorm)
library(testit)
################################ INTERMEDIATE TERMS ################################

.d_key <- c(m= ".mean", s=".sd" )


#' Return the conditional density of a normal random variable, given that its value is between .a and .b, evaluated at Za.
#'
#' @param .mean Mean of normal distribution.
#' @param .sd Standard deviation of normal distribution.
#' @param .a Unstandardized lower truncation limit.
#' @param .b Unstandardized upper truncation limit.
#' @importFrom stats dnorm pnorm

alpha_a <- function(.mean, .sd, .a, .b) {

  Z_score <- function(mean, sd, x){ (x-mean)/sd }

  Za = Z_score(.mean, .sd, .a)
  Zb = Z_score(.mean, .sd, .b)

  dnorm(Za) / (pnorm(Zb) - pnorm(Za))
}

#' Return the conditional density of a normal random variable, given that its value is between .a and .b, evaluated at Zb.
#'
#' @param .mean Mean of normal distribution.
#' @param .sd Standard deviation of normal distribution.
#' @param .a Unstandardized lower truncation limit.
#' @param .b Unstandardized upper truncation limit.
#' @importFrom stats dnorm pnorm


alpha_b <- function(.mean, .sd, .a, .b) {

  Z_score <- function(mean, sd, x){ (x-mean)/sd }

  Za = Z_score(.mean, .sd, .a)
  Zb = Z_score(.mean, .sd, .b)

  dnorm(Za) / (pnorm(Zb) - pnorm(Za))
}



#' Calculate the negative log posterior under the Jeffreys prior
#'
#' @param mean Mean of normal distribution.
#' @param sd Standard deviation of normal distribution.
#' @param .x Data to use log likelihood calculation.
#' @param .a Left truncation limit.
#' @param .b Right truncation limit.
#'
#' @example
#' nlpost_jeffreys(mean = 1, sd = 1, .x = 2, .a = -5, .b = 5)
#'
#' @importFrom mvtnorm dmvnorm pmvnorm

nlpost_jeffreys = function(mean, sd, .x, .a, .b) {
  assert("Right truncation point must be larger than left truncation point" , .a < .b )

  #MM: Leon, please change .mu and .sigma below to mean and sd instead of this hack I introduced
  .mu = mean
  .sigma = sd

  if ( .sigma < 0 ) return(.Machine$integer.max)

  # as in nll()
  term1 = dmvnorm(x = as.matrix(.x, nrow = 1),
                  mean = as.matrix(.mu, nrow = 1),
                  # sigma here is covariance matrix,
                  sigma = as.matrix(.sigma^2, nrow=1),
                  log = TRUE)


  term2 = length(.x) * log( pmvnorm(lower = .a,
                                    upper = .b,
                                    mean = .mu,
                                    # remember sigma here is covariance matrix, not the SD
                                    sigma = .sigma^2 ) )

  term3 = log( sqrt( det( E_fisher(.mean = .mu, .sd = .sigma, .n = length(.x), .a = .a, .b = .b) ) ) )

  nlp.value = -( sum(term1) - term2 + term3 )

    if ( is.infinite(nlp.value) | is.na(nlp.value) ) {
      return(.Machine$integer.max)
    }

  nlp.value

}

#' Finds the MAP estimates for mu and sigma of the full normal distribution, given data from a truncated normal.
#'
#' @param x Vector of observations from truncated normal
#' @param mean.start Initial value for mu.
#' @param sd.start Initial value for sigma.
#' @param ci.level Number between 0.5 and 1. Gives a 100(ci.level)% confidence interval.
#' @param a Left truncation limit.
#' @param b Right truncation limit.
#' @param ... Parameters to pass to sampling()
#' @example
#'
#' a <- 1
#' b <- 5
#'
#' mean.start <- 1
#' sd.start <- 1
#' iter <- 5000
#' max_treedepth <- 1
#'
#'
#' # Notice that everything following b is there as part of the ellipsis.
#' # These are optional, and will be passed directly into sampling(). See references
#' # for additional information regarding sampling().
#' estimate_jeffreys_mcmc(c(1.2,2.2,3.2), mean.start, sd.start, ci.level=0.05, a, b,
#'                        iter = iter)
#'
#' @importFrom rstan stan_model sampling
#' @importFrom stats median quantile coef
#' @importFrom stats4 mle
#'
#' @references
#' https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html

estimate_jeffreys_mcmc <- function(x,
                                  mean.start = 0,
                                  sd.start = 1,
                                  ci.level = 0.95,
                                  a,
                                  b,
                                  ...) {
  assert("Feasible standard deviation starting point ", sd.start > 0)
  model.text <- "
functions{
	real jeffreys_prior(real mu, real sigma, real a, real b, int n){
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;

		mustarL = (a - mu) / sigma;
		mustarU = (b - mu) / sigma;
		// note that normal_lpdf, etc., parameterize in terms of SD, not var
		//  the (0,1) below are *not* start values for MCMC
		alphaL = exp( normal_lpdf(mustarL | 0, 1) -
	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
	                normal_lcdf(mustarL | 0, 1) ) );

		alphaU = exp( normal_lpdf(mustarU | 0, 1) -
 	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
 	                normal_lcdf(mustarL | 0, 1) ) );

		// second derivatives for Fisher info
		kmm = -n/sigma^2 + n/sigma^2 * ((alphaU-alphaL)^2 + alphaU*mustarU- alphaL*mustarL);
		kms = -2*n/sigma^2 * (alphaL - alphaU) +
	   		  n/sigma^2 * (alphaL - alphaU + (alphaU*mustarU^2 - alphaL*mustarL^2) +
			  				(alphaL-alphaU) * (alphaL*mustarL - alphaU*mustarU));
		kss = n/sigma^2 - 3*n/sigma^2 * (1 + mustarL*alphaL - mustarU*alphaU) +
	   			n/sigma^2 * (mustarU*alphaU*(mustarU^2 - 2) - mustarL*alphaL*(mustarL^2 - 2) +
								(alphaU*mustarU - alphaL*mustarL)^2);

		fishinfo[1,1] = -kmm;
		fishinfo[1,2] = -kms;
		fishinfo[2,1] = -kms;
		fishinfo[2,2] = -kss;

		return sqrt(determinant(fishinfo));
	}
}
data{
	int<lower=0> n;
    real a;
	real b;
	real<lower=a,upper=b> y[n];
}
parameters{
    real mu;
	real<lower=0> sigma;
}
model{
	target += log( jeffreys_prior(mu, sigma, a, b, n) );
	for(i in 1:n)
        y[i] ~ normal(mu, sigma)T[a,b];
}
generated quantities{
  real log_lik;
  real log_prior = log(jeffreys_prior(mu, sigma, a, b, n));
  real log_post;
  log_lik = normal_lpdf(y | mu, sigma);
  log_lik += -n * log_diff_exp( normal_lcdf(b | mu, sigma), normal_lcdf(a | mu, sigma) );
  log_post = log_lik + log_prior;
}
"

# prepare to capture warnings from Stan
stan.warned = 0
stan.warning = NA

# set start values for sampler
init.fcn <- function(o){ list(mean=mean.start, sd=sd.start) }

# like tryCatch, but captures warnings without stopping the function from
#  returning its results
withCallingHandlers({

  # "isystem" arg is just a placeholder to avoid Stan's not understanding special characters
  #  in getwd(), even though we don't actually use the dir at all
  stan.model <- stan_model(model_code = model.text,
                           isystem = "~/Desktop")

  post = sampling(stan.model,
                  cores = 1,
                  refresh = 0,
                  init = init.fcn, data=c(n=length(x), a=a, b=b, y=x), ...)


}, warning = function(condition){
  stan.warned <<- 1
  stan.warning <<- condition$message
} )


postSumm <- summary(post)$summary

nlpost_simple = function(.mean, .sd, x, a, b) {
  nlpost.value = nlpost_jeffreys(.pars = c(.mean, .sd),
                                 .x = x, .a = a, .b = b)
  return(nlpost.value)
}

#bm

# as confirmed in "2021-8-19 Investigate profile penalized LRT inference",
#  the "MLEs" from this match those from optim() above
# https://stat.ethz.ch/R-manual/R-patched/library/stats4/html/mle.html
res <- mle( minuslogl = nlpost_simple,
           start = list( .mu=mean.start, .sigma=sd.start) )


maps <- as.numeric(coef(res))

# posterior means, then medians
Mhat = c( postSumm["mu", "mean"], median( rstan::extract(post, "mu")[[1]] ) )
Shat = c( postSumm["sigma", "mean"], median( rstan::extract(post, "sigma")[[1]] ) )

# SEs
MhatSE = postSumm["mu", "se_mean"]
ShatSE = postSumm["sigma", "se_mean"]

# convert the numeric ci.level to strings of percentages.
l.lim.str <- paste0(toString(1-ci.level * 100), "%")
r.lim.str <- paste0(toString((ci.level) * 100), "%")

# CI limits
S.CI = c( postSumm["sigma", l.lim.str], postSumm["sigma", r.lim.str] )
M.CI = c( postSumm["mu", l.lim.str], postSumm["mu", r.lim.str] )


myMhatCI = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], 1-ci.level ),
                          quantile( rstan::extract(post, "mu")[[1]], ci.level ) ) )


# the point estimates are length 2 (post means, then medians),
#  but the inference is the same for each type of point estimate
return( list( post = post,
              mean.est = Mhat,
              sd.est = Shat,

              MhatSE = rep(MhatSE, 2),
              ShatSE = rep(ShatSE, 2),

              M.CI = M.CI,
              S.CI = S.CI,

              stan.warned = stan.warned,
              stan.warning = stan.warning,

              MhatRhat = postSumm["mu", "Rhat"],
              ShatRhat = postSumm["sigma", "Rhat"]
) )
}

#' Finds the Fisher information matrix contained in n samples from a truncated normal distribution.

#' @param .mean Mean of underlying normal distribution.
#' @param .sd Standard deviation of underlying normal distribution.
#' @param .n Number of observations.
#' @param .a Lower truncation limit.
#' @param .b Upper truncation limit.
#' @importFrom assert assert
#' @importFrom stats dnorm pnorm

E_fisher = function(.mean, .sd, .n, .a, .b) {
  #assert("Positive standard deviation: ", .sd > 0)
  #assert("Left truncation before right truncation: ", .a < .b)

  Za = (.a - .mean) / .sd
  Zb = (.b - .mean) / .sd

  alpha.a = dnorm(Za) / ( pnorm(Zb) - pnorm(Za) )
  alpha.b = dnorm(Zb) / ( pnorm(Zb) - pnorm(Za) )

  k11 = -(.n/.sd^2) + (.n/.sd^2)*( (alpha.b - alpha.a)^2 + (alpha.b*Zb - alpha.a*Za) )

  k12 = -( 2*.n*(alpha.a - alpha.b) / .sd^2 ) +
    (.n/.sd^2)*( alpha.a - alpha.b + alpha.b*Zb^2 - alpha.a*Za^2 +
                      (alpha.a - alpha.b)*(alpha.a*Za - alpha.b*Zb) )

  k22 = (.n/.sd^2) - (3*.n*(1 + alpha.a*Za - alpha.b*Zb) / .sd^2) +
    (.n/.sd^2)*( Zb*alpha.b*(Zb^2 - 2) - Za*alpha.a*(Za^2 - 2) +
                      (alpha.b*Zb - alpha.a*Za)^2 )

  return( matrix( c(-k11, -k12, -k12, -k22),
                  nrow = 2,
                  byrow = TRUE ) )
}

#' Returns the value of the Jeffreys prior.
#'
#' @param .pars Vector of parameters specifying mean and standard deviation.
#' @param .x Data to use log likelihood calculation.
#' @param .a Left truncation limit.
#' @param .b Right truncation limit.

prior = function(.pars, .x, .a, .b) {
  assert("Left truncation before right truncation: ", .a < .b)
  .mu = .pars[1]
  .sigma = .pars[2]


  return (log( sqrt( det( E_fisher(.mean = .mu, .sd = .sigma, .n = length(.x), .a = .a, .b = .b) ) ) ))
}

#' Returns the value of the negative log-posterior with the Jeffrey's prior.
#'
#' @param .pars Vector specifying mu and sigma.
#' @param .x Data to use log likelihood calculation.
#' @param .a Left truncation limit.
#' @param .b Right truncation limit.
#'
#' @example
#' mu <- 1
#' sigma <- 1
#' neg_log_post(c(mu, sigma), 1, -5, 5)
#'
#' @importFrom truncnorm dtruncnorm

neg_log_post = function(.pars, .x, .a, .b) {
  assert("Positive standard deviation: ", .pars[2] > 0)
  assert("Left truncation before right truncation: ", .a < .b)
  .mu = .pars[1]
  .sigma = .pars[2]

  # regular LL part
  -sum( log( dtruncnorm(x = .x, a = .a, b = .b, mean = .mu, sd = .sigma) ) ) -
    # Jeffreys part
    log( sqrt( det( E_fisher(.mean = .mu, .sd = .sigma, .n = length(.x), .a = .a, .b = .b) ) ) )
}
