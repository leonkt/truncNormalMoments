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

  z_score <- function(mean, sd, x) (x - mean) / sd

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

  z_score <- function(mean, sd, x) (x - mean) / sd

  z_a <- z_score(mean, sd, a)
  z_b <- z_score(mean, sd, b)

  dnorm(z_a) / (pnorm(z_b) - pnorm(z_a))
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
#'
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
#' trunc_est(x = c(-1,-2,1,2), mean_start = 1, sd_start = 0.5, ci_level = 0.975,
#'           a = 1, b = 5)
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
init_fcn <- function(o) list(mean = mean_start, sd = sd_start)

# like tryCatch, but captures warnings without stopping the function from
#  returning its results
withCallingHandlers({

  # "isystem" arg is just a placeholder to avoid Stan's not understanding
  # special characters in getwd(), even though we don't actually use the dir at
  # all
  stan_model <- rstan::stan_model(model_code = model_text,
                                  isystem = "~/Desktop")


  post <- rstan::sampling(stan_model,
                          cores = 1,
                          #refresh = 0,
                          init = init_fcn,
                          data = list(n = length(x), a = a, b = b, y = x), ...)


}, warning = function(condition) {
  stan_warned <<- 1
  stan_warning <<- condition$message
})

ext <- rstan::extract(post)
best_ind <- which.max(ext$log_post)

post_summ <- summary(post)$summary

# nlpost_simple <- function(mean, sd) nlpost_jeffreys(mean, sd, x, a, b)

# res <- stats4::mle(minuslogl = nlpost_simple,
#                    start = list(mean = ext$mu[best_ind],
#                                 sd = ext$sigma[best_ind]),
#                    method = "Nelder-Mead")

# maps <- as.numeric(stats4::coef(res))

# posterior means, then medians
mean_est <- median(rstan::extract(post, "mu")[[1]])
sd_est <- median(rstan::extract(post, "sigma")[[1]])

mean_maxlp <- ext$mu[best_ind]
sd_maxlp <- ext$sigma[best_ind]
# SEs
mean_se <- post_summ["mu", "se_mean"]
sd_se <- post_summ["sigma", "se_mean"]

# convert the numeric ci_level to strings of percentages.
# l.lim.str <- paste0(toString((1 - ci_level) * 100), "%")
# r.lim.str <- paste0(toString((ci_level) * 100), "%")

# CI limits
# sd_ci_lims <- c(post_summ["sigma", l.lim.str], post_summ["sigma", r.lim.str])
# mean_ci_lims <- c(post_summ["mu", l.lim.str], post_summ["mu", r.lim.str])

mean_ci <- as.numeric(
  c(quantile(rstan::extract(post, "mu")[[1]], 1 - ci_level),
    quantile(rstan::extract(post, "mu")[[1]], ci_level))
)

sd_ci <- as.numeric(
  c(quantile(rstan::extract(post, "sigma")[[1]], 1 - ci_level),
    quantile(rstan::extract(post, "sigma")[[1]], ci_level))
)


# the point estimates are length 2 (post means, then medians), but the inference
# is the same for each type of point estimate

return(list(post = post,
            mean_est = mean_est,
            sd_est = sd_est,

            mean_maxlp = mean_maxlp,
            sd_maxlp = sd_maxlp,
            mean_se = rep(mean_se, 2),
            sd_se = rep(sd_se, 2),

            mean_ci = mean_ci,
            sd_ci = sd_ci,

            stan_warned = stan_warned,
            stan_warning = stan_warning,

            mean_rhat = post_summ["mu", "Rhat"],
            sd_rhat = post_summ["sigma", "Rhat"]))
}

#' Finds the Fisher information matrix contained in n samples from a truncated
#' normal distribution.

#' @param mean Mean of underlying normal distribution.
#' @param sd Standard deviation of underlying normal distribution.
#' @param n Number of observations.
#' @param a Lower truncation limit.
#' @param b Upper truncation limit.
#'
#' @importFrom stats dnorm pnorm
#' @export
#'
#' @examples
#' e_fisher(mean = 1, sd = 1, n = 10, a = -5, b = 5)
e_fisher <- function(mean, sd, n, a, b) {
  stopifnot(sd > 0)
  stopifnot(a < b)

  z_a <- (a - mean) / sd
  z_b <- (b - mean) / sd

  alpha_a <- alpha_a(mean, sd, a, b)
  alpha_b <- alpha_b(mean, sd, a, b)

  k11 <- -(n / sd ^ 2) + (n / sd ^ 2) * ((alpha_b - alpha_a) ^ 2 +
                                           (alpha_b * z_b - alpha_a * z_a))

  k12 <- -(2 * n * (alpha_a - alpha_b) / sd ^ 2) +
    (n / sd ^ 2) * (alpha_a - alpha_b + alpha_b * z_b ^ 2 - alpha_a * z_a ^ 2 +
                      (alpha_a - alpha_b) * (alpha_a * z_a - alpha_b * z_b))

  k22 <- (n / sd ^ 2) - (3 * n * (1 + alpha_a * z_a - alpha_b * z_b) / sd ^ 2) +
    (n / sd ^ 2) * (z_b * alpha_b * (z_b ^ 2 - 2) - z_a * alpha_a *
                      (z_a ^ 2 - 2) +
                      (alpha_b * z_b - alpha_a * z_a) ^ 2)

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
