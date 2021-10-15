

library(testthat)
library(rstan)
library(Deriv)
library(truncnorm)
library(testit)
################################ INTERMEDIATE TERMS ################################

#' Returns the standardized lower truncation limit for a normal distribution.
#'
#' @param m Mean of normal distribution.
#' @param s Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @example
#' #Returns the z-score of the value 1, for a normal distribution of mean 0, standard deviation 2.
#' usl(0, 2, 1)
#'

usl <- function(m,s,ul,uh) {
  (ul-m)/s
}

#' Returns the standardized upper truncation limit for a normal distribution.
#'
#' @param m Mean of normal distribution.
#' @param s Standard deviation of normal distribution.
#' @param uh Unstandardized upper truncation limit.
#' @example
#' # Returns the z-score of the value 1, for a normal distribution of mean 0, standard deviation 2.
#' ush(0, 2, 1)

ush <- function(m,s,ul,uh) {
  (uh-m)/s
}

#' Return the conditional density of a normal random variable, given that its value is between ul and uh, evaluated at usl.
#'
#' @param m Mean of normal distribution.
#' @param s Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.


alphal <- function(m, s, ul, uh) {
  dnorm(usl(m,s,ul,uh)) / (pnorm(ush(m,s,ul,uh)) - pnorm(usl(m,s,ul,uh)))
}

#' Return the conditional density of a normal random variable, given that its value is between ul and uh, evaluated at ush.
#'
#' @param m Mean of normal distribution.
#' @param s Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.


alphah <- function(m, s, ul, uh) {
  dnorm(ush(m,s,ul,uh)) / (pnorm(ush(m,s,ul,uh)) - pnorm(usl(m,s,ul,uh)))
}

#' Calculates factor of log-likelihood independent of the data.
#' We assume that the data are distributed according to a normal distribution of mean m, standard
#' deviation s, truncated at lower and upper limits, ul and uh, respectively.
#'
#' @param m Mean of normal distribution.
#' @param s Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#' @param n Number of observations.


ll <- function(m, s, ul, uh, n) {
  -n * log(s * sqrt(2 * pi)) - n * log(pnorm((uh-m)/s) - pnorm((ul-m)/s))
}

#' Takes samples from a normal distribution of mean m and standard deviation s.
#' Rejects samples greater than uh or less than ul.
#'
#' @param nsamp Number of observations to draw from normal distribution.
#' @param mu Mean of normal distribution
#' @param sigma Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.

samp_trunc <- function(nsamp, mu, sigma, ul, uh) {
  p_low <- pnorm(ul, mu, sigma)
  p_high <- pnorm(uh, mu, sigma)

  # draw quantiles uniformly between the limits and pass these
  # to the relevant quantile function.
  qnorm(runif(nsamp, p_low, p_high), mu, sigma)
}



#' Finds the derivative of the factor of log-likelihood independent of data.
#'
#' @param dstr Type of derivative to take. Consists of characters "m" and "s"
#' @example
#' # d_str should consist of the characters "m" and "s" indicating what types of partial derivatives to take.
#' # The length of d_str indicates the order of the derivative.
#' # For example, the following calculates the third-order derivative of log-likelihood with respect to mu, then sigma, then sigma.
#' d_str = "mss"
#' get_deriv(d_str)

get_deriv <- function(d_str) {
  # First order derivatives of the ll without the X term.
  if (nchar(d_str) == 1) {
    Deriv(ll, substr(d_str, 1,1))
  }
  # Second order derivatives of the ll without the X term.
  else if (nchar(d_str) == 2) {
    Deriv(Deriv(ll, substr(d_str, 1,1)), substr(d_str, 2,2))
  }
  # Third order derivatives of the ll without the X term.
  else if (nchar(d_str) == 3) {
    Deriv(Deriv(Deriv(ll, substr(d_str, 1,1)),substr(d_str, 2,2)), substr(d_str, 3,3))
  }
  else {
    NULL
  }
}



#' Finds cumulants of the full log-likelihood.
#'
#' @param d_str Type of derivative to take. Consists of characters "m" and "s"
#' @param m Mean of underlying normal distribution.
#' @param s Standard deviation of underlying normal distributin.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#' @example
#' # d_str should consist of the characters "m" and "s" indicating what types of partial derivatives to take.
#' # The length of d_str indicates the order of the derivative.
#' # For example, the following calculates the third-order derivative of log-likelihood with respect to mu, then sigma, then sigma.
#' d_str = "mss"
#' get_cumulants(d_str, 0, 1, 0, 1)

get_cumulants <- function(d_str, m, s, ul, uh, n) {
  # This is the NUMERIC VALUE of the derivative of the ll without the X term.
  deriv_no_x = get_deriv(d_str)(m, s, ul, uh, n)
  # Second order cumulants.
  if (nchar(d_str) == 2) {
    if (d_str == "mm"){
      deriv_no_x - (n/s**2)
    }
    else if (d_str %in% c("ms", "sm")){
      deriv_no_x - (2*n/s**3) * (etruncnorm(mean=m, sd=s, a=ul, b=uh) - m)
    }
    else if (d_str == "ss") {
      deriv_no_x - (3*n/s**4) * (s**2 * (1 + usl(m, s, ul, uh) * alphal(m, s, ul, uh) - ush(m, s, ul, uh) * alphah(m, s, ul, uh)))
    }
    else {
      NULL
    }
  }
  # Third order cumulants.
  else if (nchar(d_str) == 3) {
    if (d_str == "mmm") {
      deriv_no_x
    }
    else if (d_str %in% c("mms", "msm", "smm")) {
      deriv_no_x + (2*n/s**3)
    }
    else if (d_str %in% c("mss", "ssm", "sms")) {
      deriv_no_x + (6*n/s**4) * (etruncnorm(mean=m, sd=s, a=ul, b=uh) - m)
    }
    else if (d_str == "sss"){
      deriv_no_x + (12*n/s**5) * (s**2 * (1 + usl(m, s, ul, uh) * alphal(m, s, ul, uh) - ush(m, s, ul, uh) * alphah(m, s, ul, uh)) )
    }
    else {
      NULL
    }
  }
  else {
    NULL
  }
}

#' Get derivatives w.r.t var of the cumulants.
#' Evaluates derivatives of cumulants at specific value.
#'
#'
#' @param d_str Cumulant to take the derivative of
#' @param var Derivative of cumulant with respect to var
#'
#' @example
#' # In each case, we make a function f. f is the sum of
#' # (1) get_deriv(d_str), which gives an expression for the appropriate derivative of the LOG LIKELIHOOD WITHOUT THE X TERM.
#' # (2) expectation of the term involving x, in the appropriate derivative of the FULL LOG LIKELIHOOD.
#'
#' #For example:
#' # d^3\sigma/d\sigma^3 = (some terms not dependent on x) - (3/\sigma^4) * (\sum(x_i - \mu)^2).
#'
#' # Therefore, taking expectations will keep the terms not dependent on x the same, while changing the latter term dependent on x
#' # to -(3n/\sigma^4) * (\sigma(1 + usl * alphal - ush * alphah)). This sum is our cumulant, which we denote as f in the function.
#'
#' #We then take the derivative of f w.r.t var, since f is equal to the appropriate cumulant.

get_cumulants_deriv <- function(d_str, var, m, s, ul, uh, n) {


  cumulant_deriv <- function(d_str, var, m, s, ul, uh, n) {
    if (d_str == "mmm") {

      Deriv(get_deriv(d_str), var)
    }
    else if (d_str %in% c("mms", "msm", "smm")) {
      f <- function (m, s, ul, uh, n) {
        get_deriv(d_str)  + (2*n/s**3)
      }
      Deriv(f, var)
    }
    else if (d_str == "mm") {
      f <- function (m, s, ul, uh, n){ get_deriv(d_str)  - (n/s**2)}
      Deriv(f, var)
    }
    else if (d_str == "ss") {
      f <- function (m, s, ul, uh, n) {get_deriv(d_str)  - (3*n/s**4) * (s**2 * (1 + usl(m, s, ul, uh) * alphal(m, s, ul, uh) - ush(m, s, ul, uh) * alphah(m, s, ul, uh)) )}
      Deriv(f, var)
    }
    else if (d_str == "sss") {
      f <- function (m, s, ul, uh, n) {
        get_deriv(d_str)  + (12*n/s**5) * (s ** 2) * (s**2 * (1 + usl(m, s, ul, uh) * alphal(m, s, ul, uh) - ush(m, s, ul, uh) * alphah(m, s, ul, uh)))
      }
      Deriv(f, var)
    }
    else if (d_str %in% c("sm", "ms")) {
      f <- function (m, s, ul, uh, n) {
        get_deriv(d_str)- (2*n/s**3) * (s * (alphal(m, s, ul, uh) - alphah(m, s, ul, uh)) )
      }
      Deriv(f, var)
    }
    else if (d_str %in% c("sms", "mss", "ssm")) {
      f <- function (m,s,ul,uh, n) {
        get_deriv(d_str)  + (6*n/s**4) * (s * (alphal(m, s, ul, uh) - alphah(m, s, ul, uh)))
      }
      Deriv(f, var)
    }
    else {
      NULL
    }
  }
  cumulant_deriv(d_str, var, m, s, ul, uh, n)(m,s,ul,uh, n)
}

#' Finds the Fisher Information matrix on n observations from a truncated normal distribution.
#'
#' @param m Mean of underlying normal distribution.
#' @param s Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.

find_K <- function( m, s, ul, uh, n) {
  A <- -get_cumulants("mm",  m, s, ul, uh, n)
  B <- -get_cumulants("sm",  m, s, ul, uh, n)
  C <- -get_cumulants("ss",  m, s, ul, uh, n)
  K <- matrix(c(A, B ,B ,C), nrow=2, ncol=2, byrow=TRUE)
  K
}

#' Finds the A matrix from Godwin et al.
#'
#' @param m Mean of underlying normal distribution.
#' @param s Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.
#'
#' @references
#' Ryan T. Godwin (2016) Bias reduction for the maximum likelihood estimator of the doubly-truncated Poisson distribution,
#' Communications in Statistics - Theory and Methods, 45:7, 1887-1901, DOI: 10.1080/03610926.2013.867999

find_A <- function( m, s, ul, uh, n) {
  a_mm_s = get_cumulants_deriv("mm", "s",  m, s, ul, uh, n) - 0.5 * get_cumulants("mms",  m, s, ul, uh, n)
  a_ms_s = get_cumulants_deriv("ms", "s",  m, s, ul, uh, n) - 0.5 * get_cumulants("mss",  m, s, ul, uh, n)
  a_ss_s = get_cumulants_deriv("ss", "s",  m, s, ul, uh, n) - 0.5 * get_cumulants("sss",  m, s, ul, uh, n)

  a_mm_m = get_cumulants_deriv("mm", "m",  m, s, ul, uh, n) - 0.5 * get_cumulants("mmm",  m, s, ul, uh, n)
  a_ms_m = get_cumulants_deriv("ms", "m",  m, s, ul, uh, n) - 0.5 * get_cumulants("msm",  m, s, ul, uh, n)
  a_ss_m = get_cumulants_deriv("ss", "m",  m, s, ul, uh, n) - 0.5 * get_cumulants("ssm",  m, s, ul, uh, n)

  A_mu = matrix(c(a_mm_m, a_ms_m, a_ms_m, a_ss_m), nrow=2, ncol=2, byrow=TRUE)
  A_sigma = matrix(c(a_mm_s, a_ms_s, a_ms_s, a_ss_s), nrow=2, ncol=2, byrow=TRUE)
  cbind(A_mu, A_sigma)
}

#' Applies the bias correction, K^{-1}AK, of Godwin et al. and
#' calculates bias of the bias-corrected MLE.
#'
#' @param m Mean of underlying normal distribution.
#' @param s Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.
#' @param sim.reps Number of repetitions to estimate bias of the bias-corrected MLE.
#'
#' @references
#' Ryan T. Godwin (2016) Bias reduction for the maximum likelihood estimator of the doubly-truncated Poisson distribution,
#' Communications in Statistics - Theory and Methods, 45:7, 1887-1901, DOI: 10.1080/03610926.2013.867999

calculate_bias <- function( m, s, ul, uh, n, sim.reps) {

  K <- find_K( m, s, ul, uh, n)
  A <- find_A( m, s, ul, uh, n)
  bias <- inv(K) %*% A %*% c(inv(K))
  # vectors of standard MLEs
  means_b <- c()
  stddevs_b <- c()

  # vectors of Cordeira-corrected MLEs
  means_ub <- c()
  stddevs_ub <- c()

  # run simulation
  for (i in 1:sim.reps) {
    samps <- samp_trunc(n, m, s, ul, uh)
    mles = coef(mle.tmvnorm(as.matrix(samps, ncol = 1),
                            lower = ul,
                            upper = uh ))
    # MLE for stddev and mean of the truncated normal distribution
    mle_mean <- mles["mu_1"]
    mle_stddev <- sqrt(mles["sigma_1.1"])
    params_b <- c(mle_mean, mle_stddev)
    params_ub <- c(mle_mean, mle_stddev) - bias

    means_b <- c(means_b, params_b[1])
    stddevs_b <- c(stddevs_b, params_b[2])

    means_ub <- c(means_ub, params_ub[1])
    stddevs_ub <- c(stddevs_ub, params_ub[2])
  }

  if (type == "meanub") {
    mean(means_ub-m)

  }
  else if (type == "meanb") {
    mean(means_b-m)
  }
  else if (type == "stddevub") {
    mean(stddevs_ub-s)
  }
  else {
    mean(stddevs_b-s)
  }
}



#' Calculates the Cordiero bias-correction for the MLE of the mean and standard deviation of a
#' doubly-truncated normal.
#'
#' @param m Mean of underlying normal distribution.
#' @param s Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.
#'
cordeiro_bias = function(.m, .s, .ul, .uh, .n) {
  if (.s <= 0) {

  }
  K <- find_K( m=.m, s=.s, ul=.ul, uh=.uh, n=.n)
  A <- find_A( m=.m, s=.s, ul=.ul, uh=.uh, n=.n)
  bias <- inv(K) %*% A %*% c(inv(K))
  return( list(bias = bias, K = K, A = A) )
}

#' Auxillary function to calculate the negative log posterior with Jeffrey's prior.
#'
#' @param .pars Vector of parameters specifying mu and sigma.
#' @param par2is "sd" if the second entry in .pars is the standard deviation, and "var" otherwise.
#' @param .x Data to use log likelihood calculation.
#' @param .a Left truncation limit.
#' @param .b Right truncation limit.

nlpost_Jeffreys = function(.pars, par2is = "sd", .x, .a, .b) {
  assert("Left truncation is before right truncation" , .a < b )
  # variance parameterization
  if (par2is == "var") {

    .mu = .pars[1]
    .var = .pars[2]

    if ( .var < 0 ) return(.Machine$integer.max)

    # as in nll()
    term1 = dmvnorm(x = as.matrix(.x, nrow = 1),
                    mean = as.matrix(.mu, nrow = 1),
                    # sigma here is covariance matrix
                    sigma = as.matrix(.var, nrow=1),
                    log = TRUE)


    term2 = length(.x) * log( pmvnorm(lower = .a,
                                      upper = .b,
                                      mean = .mu,
                                      # remember sigma here is covariance matrix, not the SD
                                      sigma = .var ) )

    term3 = log( sqrt( det( E_fisher(.mu = .mu, .sigma = sqrt(.var), .n = length(.x), .a = .a, .b = .b) ) ) )

    nlp.value = -( sum(term1) - term2 + term3 )

    if ( is.infinite(nlp.value) | is.na(nlp.value) ) return(.Machine$integer.max)
  }

  # SD parameterization
  if (par2is == "sd") {

    .mu = .pars[1]
    .sigma = .pars[2]

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

    term3 = log( sqrt( det( E_fisher(.mu = .mu, .sigma = .sigma, .n = length(.x), .a = .a, .b = .b) ) ) )

    nlp.value = -( sum(term1) - term2 + term3 )

    if ( is.infinite(nlp.value) | is.na(nlp.value) ) return(.Machine$integer.max)
  }

  nlp.value

}

#' Finds the MAP estimates for mu and sigma of the full normal distribution, given data from a truncated normal.
#'
#' @param x Data from truncated normal
#' @param p Vector of parameters containing truncation points and number of observations.
#' @param mu.start Initial value for mu.
#' @param sigma.start Initial value for sigma.
#' @param ci.left String formatted as "X%". Left end of a confidence interval for each parameter estimate.
#' @param ci.right String formatted as "X%". Right end of a confidence interval for each parameter estimate.
#'
#' @examples
#'
#' n <- 3
#' a <- 1
#' b <- 2
#' iter <- 5000
#' max_treedepth <- 10
#'
#' # Notice that everything following ci.right is there as part of the ellipsis.
#' # These are optional, and will be passed directly into sampling(). See references
#' # for additional information regarding sampling().
#' estimate_jeffreys_mcmc(x=c(1.2,2.2,3.2), mu.start, sigma.start, "5%", ci.right="95%",
#'                        data = list( n = n, LL = a, UU = b),
#'                        iter = 10)
#'
#' @references
#' https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html

estimate_jeffreys_mcmc <- function(x,
                                  #p,
                                  mu.start,
                                  sigma.start,
                                  ci.left,
                                  ci.right,
                                  ...) {
  assert("Feasible standard deviation starting point ", sigma.start > 0)
  # LL and UU: cutpoints on RAW scale, not Z-scores
  # sigma: SD, not variance
  model.text <- "
functions{
	real jeffreys_prior(real mu, real sigma, real LL, real UU, int n){
		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;

		mustarL = (LL - mu) / sigma;
		mustarU = (UU - mu) / sigma;
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
    real LL;
	real UU;
	real<lower=LL,upper=UU> y[n];
}
parameters{
    real mu;
	real<lower=0> sigma;
}
model{
	target += log( jeffreys_prior(mu, sigma, LL, UU, n) );
	for(i in 1:n)
        y[i] ~ normal(mu, sigma)T[LL,UU];
}
generated quantities{
  real log_lik;
  real log_prior = log(jeffreys_prior(mu, sigma, LL, UU, n));
  real log_post;
  log_lik = normal_lpdf(y | mu, sigma);
  log_lik += -n * log_diff_exp( normal_lcdf(UU | mu, sigma), normal_lcdf(LL | mu, sigma) );
  log_post = log_lik + log_prior;
}
"

# prepare to capture warnings from Stan
stan.warned = 0
stan.warning = NA

# set start values for sampler
init.fcn <- function(o){ list(mu=mu.start, sigma=sigma.start) }

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
                  #data = list( n = p$n, LL = p$a, UU = p$b, y = x ),

                  #iter = p$stan.iter,
                  #control = list(max_treedepth = p$stan.maxtreedepth,
                                 #adapt_delta = p$stan.adapt_delta),
                  init = init.fcn, ...)


}, warning = function(condition){
  stan.warned <<- 1
  stan.warning <<- condition$message
} )


postSumm = summary(post)$summary


nlpost_simple = function(.mu, .sigma, par2is, x, a, b) {
  nlpost.value = nlpost_Jeffreys(.pars = c(.mu, .sigma),
                                 par2is = par2is,
                                 .x = x, .a = a, .b = b)
  return(nlpost.value)
}

#bm

# as confirmed in "2021-8-19 Investigate profile penalized LRT inference",
#  the "MLEs" from this match those from optim() above
# https://stat.ethz.ch/R-manual/R-patched/library/stats4/html/mle.html
res = mle( minuslogl = nlpost_simple,
           start = list( .mu=mu.start, .sigma=sigma.start) )

# not actually MLEs, of course, but rather MAPs
maps = as.numeric(coef(res))

# posterior means, then medians
Mhat = c( postSumm["mu", "mean"], median( rstan::extract(post, "mu")[[1]] ) )
Shat = c( postSumm["sigma", "mean"], median( rstan::extract(post, "sigma")[[1]] ) )
Vhat = Shat^2
# sanity check


# SEs
MhatSE = postSumm["mu", "se_mean"]
ShatSE = postSumm["sigma", "se_mean"]
# because VhatSE uses delta method, VhatSE will be length 2 because Shat is length 2
VhatSE = ShatSE * 2 * Shat
# how Stan estimates the SE: https://discourse.mc-stan.org/t/se-mean-in-print-stanfit/2869


# CI limits
S.CI = c( postSumm["sigma", "5%"], postSumm["sigma", ci.right] )
V.CI = S.CI^2
M.CI = c( postSumm["mu", ci.left], postSumm["mu", ci.right] )
# sanity check:
l.lim <- as.numeric(substr(ci.left, 1, nchar(ci.left) - 1)) / 100
r.lim <- as.numeric(substr(ci.right, 1, ncchar(ci.right) - 1)) / 100

assert("Left endpoint is less than right endpoint", l.lim < r.lim)

myMhatCI = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], l.lim ),
                          quantile( rstan::extract(post, "mu")[[1]], r.lim ) ) )


# the point estimates are length 2 (post means, then medians),
#  but the inference is the same for each type of point estimate
return( list( post = post,
              Mhat = Mhat,
              Vhat = Vhat,
              Shat = Shat,

              map = maps,

              MhatSE = rep(MhatSE, 2),
              VhatSE = VhatSE,  # already length 2 (see above)
              ShatSE = rep(ShatSE, 2),

              M.CI = M.CI,
              V.CI = V.CI,
              S.CI = S.CI,

              stan.warned = stan.warned,
              stan.warning = stan.warning,

              MhatRhat = postSumm["mu", "Rhat"],
              ShatRhat = postSumm["sigma", "Rhat"]
) )
}

#' Finds the Fisher information matrix contained in n samples from a truncated normal distribution.
#'
#' @param .mu Mean of underlying normal distribution.
#' @param .sigma Standard deviation of underlying normal distribution.
#' @param .n Number of observations.
#' @param .a Lower truncation limit.
#' @param .b Upper truncation limit.

E_fisher = function(.mu, .sigma, .n, .a, .b) {
  assert("Positive standard deviation: ",.sigma > 0)
  assert("Left truncation before right truncation: ", .a < .b)

  Za = (.a - .mu) / .sigma
  Zb = (.b - .mu) / .sigma

  alpha.a = dnorm(Za) / ( pnorm(Zb) - pnorm(Za) )
  alpha.b = dnorm(Zb) / ( pnorm(Zb) - pnorm(Za) )

  k11 = -(.n/.sigma^2) + (.n/.sigma^2)*( (alpha.b - alpha.a)^2 + (alpha.b*Zb - alpha.a*Za) )

  k12 = -( 2*.n*(alpha.a - alpha.b) / .sigma^2 ) +
    (.n/.sigma^2)*( alpha.a - alpha.b + alpha.b*Zb^2 - alpha.a*Za^2 +
                      (alpha.a - alpha.b)*(alpha.a*Za - alpha.b*Zb) )

  k22 = (.n/.sigma^2) - (3*.n*(1 + alpha.a*Za - alpha.b*Zb) / .sigma^2) +
    (.n/.sigma^2)*( Zb*alpha.b*(Zb^2 - 2) - Za*alpha.a*(Za^2 - 2) +
                      (alpha.b*Zb - alpha.a*Za)^2 )

  return( matrix( c(-k11, -k12, -k12, -k22),
                  nrow = 2,
                  byrow = TRUE ) )
}

#' Returns the value of the Jeffrey's prior.
#'
#' @param .pars Vector of parameters specifying mu and sigma.
#' @param .x Data to use log likelihood calculation.
#' @param .a Left truncation limit.
#' @param .b Right truncation limit.

prior = function(.pars, .x, .a, .b) {
  assert("Left truncation before right truncation: ", .a < .b)
  .mu = .pars[1]
  .sigma = .pars[2]


  return (log( sqrt( det( E_fisher(.mu = .mu, .sigma = .sigma, .n = length(.x), .a = .a, .b = .b) ) ) ))
}

#' Returns the value of the negative log-posterior with the Jeffrey's prior.
#'
#' @param .pars Vector of parameters specifying mu and sigma.
#' @param .x Data to use log likelihood calculation.
#' @param .a Left truncation limit.
#' @param .b Right truncation limit.

neg_log_post = function(.pars, .x, .a, .b) {
  assert("Positive standard deviation: ", .pars[2] > 0)
  assert("Left truncation before right truncation: ", .a < .b)
  .mu = .pars[1]
  .sigma = .pars[2]

  # regular LL part
  -sum( log( dtruncnorm(x = .x, a = .a, b = .b, mean = .mu, sd = .sigma) ) ) -
    # Jeffreys part
    log( sqrt( det( E_fisher(.mu = .mu, .sigma = .sigma, .n = length(.x), .a = .a, .b = .b) ) ) )
}
