
library(testthat)
library(rstan)
library(Deriv)
library(truncnorm)
library(testit)
################################ INTERMEDIATE TERMS ################################

.d_key <- c(m= ".mean", s=".sd" )


#' Returns the standardized lower truncation limit for a normal distribution.
#'
#' @param .mean Mean of normal distribution.
#' @param .sd Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#' @example
#' # Returns the z-score of the value 1, for a normal distribution of mean 0, standard deviation 2.
#' usl(0, 2, 1)
#'

usl <- function(.mean,.sd,ul,uh) {
  (ul-.mean)/.sd
}

#' Returns the standardized upper truncation limit for a normal distribution.
#'
#' @param .mean Mean of normal distribution.
#' @param .sd Standard deviation of normal distribution.
#' @param uh Unstandardized upper truncation limit.
#' @param ul Unstandardized lower truncation limit.
#' @example
#' # Returns the z-score of the value 1, for a normal distribution of mean 0, standard deviation 2.
#' ush(0, 2, 1)

ush <- function(.mean, .sd,ul,uh) {
  (uh-.mean)/.sd
}

#' Return the conditional density of a normal random variable, given that its value is between ul and uh, evaluated at usl.
#'
#' @param .mean Mean of normal distribution.
#' @param .sd Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#' @importFrom stats dnorm pnorm


alphal <- function(.mean, .sd, ul, uh) {
  dnorm(usl(.mean,.sd,ul,uh)) / (pnorm(ush(.mean,.sd,ul,uh)) - pnorm(usl(.mean,.sd,ul,uh)))
}

#' Return the conditional density of a normal random variable, given that its value is between ul and uh, evaluated at ush.
#'
#' @param .mean Mean of normal distribution.
#' @param .sd Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#' @importFrom stats dnorm pnorm


alphah <- function(.mean, .sd, ul, uh) {
  dnorm(ush(.mean,.sd,ul,uh)) / (pnorm(ush(.mean,.sd,ul,uh)) - pnorm(usl(.mean,.sd,ul,uh)))
}

#' Calculates factor of log-likelihood independent of the data.
#' We assume that the data are distributed according to a normal distribution of mean m, standard
#' deviation s, truncated at lower and upper limits, ul and uh, respectively.
#'
#' @param .mean Mean of normal distribution.
#' @param .sd Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#' @param n Number of observations.


ll <- function(.mean, .sd, ul, uh, n) {
  -n * log(.sd * sqrt(2 * pi)) - n * log(pnorm((uh-.mean)/.sd) - pnorm((ul-.mean)/.sd))
}



#' Finds the derivative of the factor of log-likelihood independent of data.
#'
#' @param d_str Type of derivative to take. Consists of characters "m" and "s"
#' @example
#' # d_str should consist of the characters "m" and "s" indicating what types of partial derivatives to take.
#' # The length of d_str indicates the order of the derivative.
#' # For example, the following calculates the third-order derivative of log-likelihood with respect to mu, then sigma, then sigma.
#' d_str = "mss"
#' get_deriv(d_str)
#' @importFrom Deriv Deriv

get_deriv <- function(d_str) {
  # First order derivatives of the ll without the X term.
  if (nchar(d_str) == 1) {
    Deriv(ll, .d_key[substr(d_str, 1,1)][[1]])
  }
  # Second order derivatives of the ll without the X term.
  else if (nchar(d_str) == 2) {
    Deriv(Deriv(ll, .d_key[substr(d_str, 1,1)][[1]]), .d_key[substr(d_str, 2,2)][[1]])
  }
  # Third order derivatives of the ll without the X term.
  else if (nchar(d_str) == 3) {
    Deriv(Deriv(Deriv(ll, .d_key[substr(d_str, 1,1)][[1]]), .d_key[substr(d_str, 2,2)][[1]]), .d_key[substr(d_str, 3,3)][[1]])
  }
  else {
    NULL
  }
}



#' Finds cumulants of the full log-likelihood.
#'
#' @param d_str Type of derivative to take. Consists of characters "m" and "s"
#' @param .mean Mean of underlying normal distribution.
#' @param .sd Standard deviation of underlying normal distributin.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#' @param n Number of observations.
#' @example
#' # d_str should consist of the characters "m" and "s" indicating what types of partial derivatives to take.
#' # The length of d_str indicates the order of the derivative.
#' # For example, the following calculates the third-order derivative of log-likelihood with respect to mu, then sigma, then sigma.
#' d_str = "mss"
#' get_cumulants(d_str, 0, 1, 0, 1)
#'
#' @importFrom truncnorm etruncnorm

get_cumulants <- function(d_str, .mean, .sd, ul, uh, n) {
  # This is the NUMERIC VALUE of the derivative of the ll without the X term.
  deriv_no_x = get_deriv(d_str)(.mean, .sd, ul, uh, n)
  # Second order cumulants.
  if (nchar(d_str) == 2) {
    if (d_str == "mm"){
      deriv_no_x - (n/(.sd)**2)
    }
    else if (d_str %in% c("ms", "sm")){
      deriv_no_x - (2*n/(.sd)**3) * (etruncnorm(mean=.mean, sd=.sd, a=ul, b=uh) - .mean)
    }
    else if (d_str == "ss") {
      deriv_no_x - (3*n/(.sd)**4) * ((.sd)**2 * (1 + usl(.mean, .sd, ul, uh) * alphal(.mean, .sd, ul, uh) - ush(.mean, .sd, ul, uh) * alphah(.mean, .sd, ul, uh)))
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
      deriv_no_x + (2*n/(.sd)**3)
    }
    else if (d_str %in% c("mss", "ssm", "sms")) {
      deriv_no_x + (6*n/(.sd)**4) * (etruncnorm(mean=.mean, sd=.sd, a=ul, b=uh) - .mean)
    }
    else if (d_str == "sss"){
      deriv_no_x + (12*n/(.sd)**5) * ((.sd)**2 * (1 + usl(.mean, .sd, ul, uh) * alphal(.mean, .sd, ul, uh) - ush(.mean, .sd, ul, uh) * alphah(.mean, .sd, ul, uh)) )
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
#' @param .mean Mean of untruncated normal distribution
#' @param .sd Standard deviation of untruncated normal distribution.
#' @param ul Lower truncation limit
#' @param uh Upper truncation limit
#' @param n Number of observations
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
#' # We then take the derivative of f w.r.t var, since f is equal to the appropriate cumulant.
#'
#' @importFrom Deriv Deriv

get_cumulants_deriv <- function(d_str, var, .mean, .sd, ul, uh, n) {


  cumulant_deriv <- function(d_str, var, .mean, .sd, ul, uh, n) {
    if (d_str == "mmm") {

      Deriv(get_deriv(d_str), .d_key[var][[1]])
    }
    else if (d_str %in% c("mms", "msm", "smm")) {
      f <- function (.mean, .sd, ul, uh, n) {
        get_deriv(d_str)  + (2*n/(.sd)**3)
      }
      Deriv(f, .d_key[var][[1]])
    }
    else if (d_str == "mm") {
      f <- function (.mean, .sd, ul, uh, n){ get_deriv(d_str)  - (n/(.sd)**2)}
      Deriv(f, .d_key[var][[1]])
    }
    else if (d_str == "ss") {
      f <- function (.mean, .sd, ul, uh, n) {get_deriv(d_str)  - (3*n/(.sd)**4) * ((.sd)**2 * (1 + usl(.mean, .sd, ul, uh) * alphal(.mean, .sd, ul, uh) - ush(.mean, .sd, ul, uh) * alphah(.mean, .sd, ul, uh)) )}
      Deriv(f, .d_key[var][[1]])
    }
    else if (d_str == "sss") {
      f <- function (.mean, .sd, ul, uh, n) {
        get_deriv(d_str)  + (12*n/(.sd)**5) * ((.sd) ** 2) * ((.sd)**2 * (1 + usl(.mean, .sd, ul, uh) * alphal(.mean, .sd, ul, uh) - ush(.mean, .sd, ul, uh) * alphah(.mean, .sd, ul, uh)))
      }
      Deriv(f, .d_key[var][[1]])
    }
    else if (d_str %in% c("sm", "ms")) {
      f <- function (.mean, .sd, ul, uh, n) {
        get_deriv(d_str)- (2*n/(.sd)**3) * ((.sd) * (alphal(.mean, .sd, ul, uh) - alphah(.mean, .sd, ul, uh)) )
      }
      Deriv(f, .d_key[var][[1]])
    }
    else if (d_str %in% c("sms", "mss", "ssm")) {
      f <- function (.mean, .sd, ul, uh, n) {
        get_deriv(d_str)  + (6*n/(.sd)**4) * ((.sd) * (alphal(.mean, .sd, ul, uh) - alphah(.mean, .sd, ul, uh)))
      }
      Deriv(f, .d_key[var][[1]])
    }
    else {
      NULL
    }
  }
  cumulant_deriv(d_str, var, .mean, .sd, ul, uh, n)(.mean, .sd, ul, uh, n)
}

#' Finds the Fisher Information matrix on n observations from a truncated normal distribution.
#'
#' @param .mean Mean of underlying normal distribution.
#' @param .sd Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.

find_K <- function(.mean, .sd, ul, uh, n) {
  A <- -get_cumulants("mm",  .mean, .sd, ul, uh, n)
  B <- -get_cumulants("sm",  .mean, .sd, ul, uh, n)
  C <- -get_cumulants("ss",  .mean, .sd, ul, uh, n)
  K <- matrix(c(A, B ,B ,C), nrow=2, ncol=2, byrow=TRUE)
  K
}

#' Finds the A matrix from Godwin et al.
#'
#' @param .mean Mean of underlying normal distribution.
#' @param .sd Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.
#'
#' @references
#' Ryan T. Godwin (2016) Bias reduction for the maximum likelihood estimator of the doubly-truncated Poisson distribution,
#' Communications in Statistics - Theory and Methods, 45:7, 1887-1901, DOI: 10.1080/03610926.2013.867999

find_A <- function( .mean, .sd, ul, uh, n) {
  a_mm_s = get_cumulants_deriv("mm", "s",  .mean, .sd, ul, uh, n) - 0.5 * get_cumulants("mms",  .mean, .sd, ul, uh, n)
  a_ms_s = get_cumulants_deriv("ms", "s",  .mean, .sd, ul, uh, n) - 0.5 * get_cumulants("mss",  .mean, .sd, ul, uh, n)
  a_ss_s = get_cumulants_deriv("ss", "s",  .mean, .sd, ul, uh, n) - 0.5 * get_cumulants("sss",  .mean, .sd, ul, uh, n)

  a_mm_m = get_cumulants_deriv("mm", "m",  .mean, .sd, ul, uh, n) - 0.5 * get_cumulants("mmm",  .mean, .sd, ul, uh, n)
  a_ms_m = get_cumulants_deriv("ms", "m",  .mean, .sd, ul, uh, n) - 0.5 * get_cumulants("msm",  .mean, .sd, ul, uh, n)
  a_ss_m = get_cumulants_deriv("ss", "m",  .mean, .sd, ul, uh, n) - 0.5 * get_cumulants("ssm",  .mean, .sd, ul, uh, n)

  A_mu = matrix(c(a_mm_m, a_ms_m, a_ms_m, a_ss_m), nrow=2, ncol=2, byrow=TRUE)
  A_sigma = matrix(c(a_mm_s, a_ms_s, a_ms_s, a_ss_s), nrow=2, ncol=2, byrow=TRUE)
  cbind(A_mu, A_sigma)
}

#' Takes samples from a normal distribution of mean m and standard deviation s.
#' Rejects samples greater than uh or less than ul.
#'
#' @param nsamp Number of observations to draw from normal distribution.
#' @param mean Mean of normal distribution
#' @param sd Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @param uh Unstandardized upper truncation limit.
#'
#' @importFrom stats pnorm qnorm runif
#'

samp_trunc <- function(nsamp, mean, sd, ul, uh) {
  p_low <- pnorm(ul, mean, sd)
  p_high <- pnorm(uh, mean, sd)

  # draw quantiles uniformly between the limits and pass these
  # to the relevant quantile function.
  qnorm(runif(nsamp, p_low, p_high), mean, sd)
}

#' Applies the bias correction, K^{-1}AK, of Godwin et al. and
#' calculates bias of the bias-corrected MLE.
#'
#' @param .mean Mean of underlying normal distribution.
#' @param .sd Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.
#' @param sim.reps Number of repetitions to estimate bias of the bias-corrected MLE.
#' @param type "meanub", "stddevub" for bias of Cordiero estimator for mean or standard deviation respectively.
#'             "meanb", "stddevb" for bias of MLE of mean and standard deviation.
#' @importFrom pracma inv
#' @importFrom stats coef
#' @importFrom tmvtnorm mle.tmvnorm

#' @references
#' Ryan T. Godwin (2016) Bias reduction for the maximum likelihood estimator of the doubly-truncated Poisson distribution,
#' Communications in Statistics - Theory and Methods, 45:7, 1887-1901, DOI: 10.1080/03610926.2013.867999

calculate_bias <- function(.mean, .sd, ul, uh, n, sim.reps, type) {

  K <- find_K(.mean, .sd, ul, uh, n)
  A <- find_A(.mean, .sd, ul, uh, n)
  bias <- inv(K) %*% A %*% c(inv(K))
  # vectors of standard MLEs
  means_b <- c()
  stddevs_b <- c()

  # vectors of Cordeiro-corrected MLEs
  means_ub <- c()
  stddevs_ub <- c()

  # run simulation
  for (i in 1:sim.reps) {
    samps <- samp_trunc(n, .mean, .sd, ul, uh)
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
    mean(means_ub-.mean)

  }
  else if (type == "meanb") {
    mean(means_b-.mean)
  }
  else if (type == "stddevub") {
    mean(stddevs_ub-.sd)
  }
  else {
    mean(stddevs_b-.sd)
  }
}



#' Calculates the Cordiero bias-correction for the MLE of the mean and standard deviation of a
#' doubly-truncated normal.
#'
#' @param .mean Mean of underlying normal distribution.
#' @param .sd Standard deviation of underlying normal distribution.
#' @param .ul Lower truncation limit.
#' @param .uh Upper truncation limit.
#' @param .n Number of observations.

cordeiro_bias = function(.mean, .sd, .ul, .uh, .n) {
  K <- find_K( .mean=.mean, .sd=.sd, ul=.ul, uh=.uh, n=.n)
  A <- find_A( .mean=.mean, .sd=.sd, ul=.ul, uh=.uh, n=.n)
  bias <- inv(K) %*% A %*% c(inv(K))
  return( list(bias = bias, K = K, A = A) )
}

#' Auxillary function to calculate the negative log posterior with Jeffrey's prior.
#'
#' @param .pars Vector of parameters specifying mean and standard deviation.
#' @param .x Data to use log likelihood calculation.
#' @param .a Left truncation limit.
#' @param .b Right truncation limit.
#'
#' @example
#' mu <- 1
#' sigma <- 1
#' x <- 2
#' a <- -5
#' b <- 5
#' nlpost_Jeffreys(pars=c(mu, sigma), x, a, b)
#'
#' @importFrom mvtnorm dmvnorm pmvnorm

nlpost_Jeffreys = function(.pars, .x, .a, .b) {
  assert("Left truncation is before right truncation" , .a < .b )
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
#' MM: Let's format ci.left and ci.right as a single numerical arg, ci.level (default 0.95), and assume they want symmetric
#' @param ci.left Number between 0 and 1 . Left end of a confidence interval for each parameter estimate.
#' @param ci.right Number between 0 and 1. Right end of a confidence interval for each parameter estimate. By default 1-ci.left.
#' @param ... Parameters to pass to sampling()
#' @example
#'
#' a <- 1
#' b <- 5
#' mean.start <- 1
#' sd.start <- 1
#' iter <- 5000
#' max_treedepth <- 1
#'
#'
#' # Notice that everything following ci.right is there as part of the ellipsis.
#' # These are optional, and will be passed directly into sampling(). See references
#' # for additional information regarding sampling().
#' estimate_jeffreys_mcmc(c(1.2,2.2,3.2), mean.start, sd.start, "5%", ci.right="95%",
#'                        data = list( n = 3, LL = a, UU = b),
#'                        iter = iter)
#'
#' @importFrom rstan stan_model sampling
#' @importFrom stats median quantile
#' @importFrom stats4 mle
#' @references
#' https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html

estimate_jeffreys_mcmc <- function(x,
                                  mean.start = 0,
                                  sd.start = 1,
                                  ci.left,
                                  ci.right,
                                  ...) {
  assert("Feasible standard deviation starting point ", sd.start > 0)
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
init.fcn <- function(o){ list(mu=mean.start, sigma=sd.start) }

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
                  init = init.fcn, ...)


}, warning = function(condition){
  stan.warned <<- 1
  stan.warning <<- condition$message
} )


postSumm <- summary(post)$summary
print(postSumm)

nlpost_simple = function(.mean, .sd, x, a, b) {
  nlpost.value = nlpost_Jeffreys(.pars = c(.mean, .sd),
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

# convert the numeric ci.left and ci.right to strings of percentages.
l.lim.str <- paste0(toString(ci.left * 100), "%")
r.lim.str <- paste0(toString(ci.right * 100), "%")

# CI limits
S.CI = c( postSumm["sigma", l.lim.str], postSumm["sigma", r.lim.str] )
M.CI = c( postSumm["mu", l.lim.str], postSumm["mu", r.lim.str] )
# sanity check:


assert("Left endpoint is less than right endpoint", ci.right < ci.left)

myMhatCI = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], ci.left ),
                          quantile( rstan::extract(post, "mu")[[1]], ci.right ) ) )


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
  assert("Positive standard deviation: ",.sd > 0)
  assert("Left truncation before right truncation: ", .a < .b)

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
