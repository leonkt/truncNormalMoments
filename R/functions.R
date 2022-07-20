#' @keywords internal
z_score <- function(mean, sd, x) (x - mean) / sd

#' Return the conditional density of a normal random variable, given that its
#' value is between a and b, evaluated at z_a.
#'
#' @param mean Mean of normal distribution.
#' @param sd Standard deviation of normal distribution.
#' @param a Unstandardized lower truncation limit.
#' @param b Unstandardized upper truncation limit.
#'
#' @importFrom stats dnorm pnorm
#' @export
#'
#' @examples
#' alpha_a(mean = 1, sd = 1, a = -5, b = 5)
alpha_a <- function(mean, sd, a, b) {
  stopifnot(a < b)
  stopifnot(sd > 0)

  z_a <- z_score(mean, sd, a)
  z_b <- z_score(mean, sd, b)

  dnorm(z_a) / (pnorm(z_b) - pnorm(z_a))
}

#' Return the conditional density of a normal random variable, given that its
#' value is between a and b, evaluated at z_b.
#'
#' @param mean Mean of normal distribution.
#' @param sd Standard deviation of normal distribution.
#' @param a Unstandardized lower truncation limit.
#' @param b Unstandardized upper truncation limit.
#'
#' @importFrom stats dnorm pnorm
#' @export
#'
#' @examples
#' alpha_b(mean = 1, sd = 1, a = -5, b = 5)
alpha_b <- function(mean, sd, a, b) {
  stopifnot(a < b)
  stopifnot(sd > 0)

  z_a <- z_score(mean, sd, a)
  z_b <- z_score(mean, sd, b)

  # dnorm(z_a) / (pnorm(z_b) - pnorm(z_a))
  dnorm(z_b) / (pnorm(z_b) - pnorm(z_a)) # TODO: is this right?
}

#' Calculate the negative log posterior under the Jeffreys prior
#'
#' @param mean Mean of normal distribution.
#' @param sd Standard deviation of normal distribution.
#' @param x Data to use log likelihood calculation.
#' @param a Left truncation limit.
#' @param b Right truncation limit.
#'
#' @export
#' @examples
#' nlpost_jeffreys(mean = 1, sd = 1, x = 2, a = -5, b = 5)
nlpost_jeffreys <- function(mean, sd, x, a, b) {
  stopifnot(a < b)

  if (sd < 0) return(.Machine$integer.max)

  # as in nll()
  term1 <- mvtnorm::dmvnorm(x = as.matrix(x, nrow = 1),
                            mean = as.matrix(mean, nrow = 1),
                            # sigma here is covariance matrix
                            sigma = as.matrix(sd ^ 2, nrow = 1),
                            log = TRUE)

  # remember sigma here is covariance matrix, not the SD
  term2 <- length(x) * log(mvtnorm::pmvnorm(lower = a,
                                            upper = b,
                                            mean = mean,
                                            sigma = sd ^ 2))

  ef <- e_fisher(mean = mean, sd = sd, n = length(x), a = a, b = b)
  term3 <- log(sqrt(det(ef)))

  nlp_value <- -(sum(term1) - term2 + term3)

  if (is.infinite(nlp_value) || is.na(nlp_value)) return(.Machine$integer.max)
  return(as.numeric(nlp_value))
}

#' Estimates the posterior modes for the parameters of the underlying normal
#' distribution, given truncated data
#'
#' @param x Vector of observations from truncated normal
#' @param mean_start Initial value for mu.
#' @param sd_start Initial value for sigma.
#' @param ci_level Number between 0.5 and 1. Gives a 100(ci_level)% confidence
#'   interval.
#' @param a Left truncation limit.
#' @param b Right truncation limit.
#' @param ... Parameters to pass to sampling()
#'
#' @importFrom stats median quantile coef
#' @export
#'
#' @examples
#' x <- truncnorm::rtruncnorm(100, a = 0, b = 2, mean = 0.5, sd = 0.5)
#' trunc_est(x, a = 0, b = 2)
#'
#' @references
#' https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html
trunc_est <- function(x,
                      mean_start = 0,
                      sd_start = 1,
                      ci_level = 0.95,
                      a,
                      b,
                      ...) {
  stopifnot(a < b)
  stopifnot(sd_start > 0)
  # TODO: assert that x values within truncation limits

  model_text <- "
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
stan_warned <- 0
stan_warning <- NA

# set start values for sampler
init_fcn <- function() list(mean = mean_start, sd = sd_start)

# like tryCatch, but captures warnings without stopping the function from
#  returning its results
withCallingHandlers({

  # "isystem" arg is just a placeholder to avoid Stan's not understanding
  # special characters in getwd(), even though we don't actually use the dir at
  # all
  stan_model <- rstan::stan_model(model_code = model_text,
                                  isystem = "~/Desktop")

  stan_fit <- rstan::sampling(stan_model,
                              cores = 1,
                              init = init_fcn,
                              data = list(n = length(x), a = a, b = b, y = x),
                              ...)


}, warning = function(condition) {
  stan_warned <<- 1
  stan_warning <<- condition$message
})

stan_extract <- rstan::extract(stan_fit)
stan_summary <- as.data.frame(rstan::summary(stan_fit)$summary[c("mu", "sigma"),])
means <- stan_summary$mean
ses <- stan_summary$se_mean
rhats <- stan_summary$Rhat

medians <- c(median(stan_extract$mu), median(stan_extract$sigma))

index_maxlp <- which.max(stan_extract$log_post)
maxlps <- c(stan_extract$mu[index_maxlp], stan_extract$sigma[index_maxlp])

stan_ci <- function(param, q) as.numeric(quantile(stan_extract[[param]], q))
cil <- c(stan_ci("mu", 1 - ci_level), stan_ci("sigma", 1 - ci_level))
ciu <- c(stan_ci("mu", ci_level), stan_ci("sigma", ci_level))

stan_stats <- data.frame(param = c("mean", "sd"), mean = means,
                         median = medians, maxlp = maxlps, se = ses,
                         ci_lower = cil, ci_upper = ciu, rhat = rhats)

return(list(stats = stan_stats, fit = stan_fit))
}

#' Finds the Fisher information matrix contained in n samples from a truncated
#' normal distribution.

#' @param mean Mean of underlying normal distribution.
#' @param sd Standard deviation of underlying normal distribution.
#' @param n Number of observations.
#' @param a Lower truncation limit.
#' @param b Upper truncation limit.
#'
#' @export
#' @examples
#' e_fisher(mean = 1, sd = 1, n = 10, a = -5, b = 5)
e_fisher <- function(mean, sd, n, a, b) {
  stopifnot(sd > 0)
  stopifnot(a < b)

  z_a <- z_score(mean, sd, a)
  z_b <- z_score(mean, sd, b)

  alp_a <- alpha_a(mean, sd, a, b)
  alp_b <- alpha_b(mean, sd, a, b)

  k11 <- -(n / sd ^ 2) + (n / sd ^ 2) *
    ((alp_b - alp_a) ^ 2 + (alp_b * z_b - alp_a * z_a))

  k12 <- -(2 * n * (alp_a - alp_b) / sd ^ 2) +
    (n / sd ^ 2) * (alp_a - alp_b + alp_b * z_b ^ 2 - alp_a * z_a ^ 2 +
                      (alp_a - alp_b) * (alp_a * z_a - alp_b * z_b))

  k22 <- (n / sd ^ 2) -
    (3 * n * (1 + alp_a * z_a - alp_b * z_b) / sd ^ 2) +
    (n / sd ^ 2) *
    (z_b * alp_b * (z_b ^ 2 - 2) - z_a * alp_a * (z_a ^ 2 - 2) +
       (alp_b * z_b - alp_a * z_a) ^ 2)

  return(matrix(c(-k11, -k12, -k12, -k22), nrow = 2, byrow = TRUE))
}

#' Returns the value of the Jeffreys prior.
#'
#' @param mean Mean of untruncated normal distribution
#' @param sd Standard deviation of untruncated normal distribution.
#' @param x Data to use log likelihood calculation.
#' @param a Left truncation limit.
#' @param b Right truncation limit.
#'
#' @export
#' @examples
#' prior_jeffreys(mean = 1, sd = 1, x = 2, a = -5, b = 5)
prior_jeffreys <- function(mean, sd, x, a, b) {
  stopifnot(a < b)
  stopifnot(sd > 0)

  ef <- e_fisher(mean = mean, sd = sd, n = length(x), a = a, b = b)
  return(log(sqrt(det(ef))))
}
