################################ INTERMEDIATE TERMS ################################

#' Return the conditional density of a normal random variable, given that its value is between a and b, evaluated at Za.
#'
#' @param mean Mean of normal distribution.
#' @param sd Standard deviation of normal distribution.
#' @param a Unstandardized lower truncation limit.
#' @param b Unstandardized upper truncation limit.
#'
#' @importFrom stats dnorm pnorm
#' @export
alpha_a <- function(mean, sd, a, b) {
  stopifnot(a < b)
  stopifnot(sd > 0)

  Z_score <- function(mean, sd, x){ (x-mean)/sd }

  Za = Z_score(mean, sd, a)
  Zb = Z_score(mean, sd, b)

  dnorm(Za) / (pnorm(Zb) - pnorm(Za))
}

#' Return the conditional density of a normal random variable, given that its value is between a and b, evaluated at Zb.
#'
#' @param mean Mean of normal distribution.
#' @param sd Standard deviation of normal distribution.
#' @param a Unstandardized lower truncation limit.
#' @param b Unstandardized upper truncation limit.
#'
#' @importFrom stats dnorm pnorm
#' @export
alpha_b <- function(mean, sd, a, b) {
  stopifnot(a < b)
  stopifnot(sd > 0)

  Z_score <- function(mean, sd, x){ (x-mean)/sd }

  Za = Z_score(mean, sd, a)
  Zb = Z_score(mean, sd, b)

  dnorm(Za) / (pnorm(Zb) - pnorm(Za))
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
#' @example
#' nlpost_jeffreys(mean = 1, sd = 1, x = c(2), a = -5, b = 5)

nlpost_jeffreys <- function(mean, sd, x, a, b) {
  stopifnot(a < b)

  if ( sd < 0 ) return(.Machine$integer.max)

  # as in nll()
  term1 <- mvtnorm::dmvnorm(x = as.matrix(x, nrow = 1),
                            mean = as.matrix(mean, nrow = 1),
                            # sigma here is covariance matrix,
                            sigma = as.matrix(sd^2, nrow=1),
                            log = TRUE)


  term2 <- length(x) * log( mvtnorm::pmvnorm(lower = a,
                                             upper = b,
                                             mean = mean,
                                             # remember sigma here is covariance matrix, not the SD
                                             sigma = sd^2 ) )

  term3 <- log( sqrt( det( E_fisher(mean = mean, sd = sd, n = length(x), a = a, b = b) ) ) )

  nlp.value <- -( sum(term1) - term2 + term3 )

  if ( is.infinite(nlp.value) | is.na(nlp.value) ) {
    return(.Machine$integer.max)
  }

  nlp.value

}

#' Estimates the posterior modes for the parameters of the underlying normal distribution, given truncated data
#'
#' @param x Vector of observations from truncated normal
#' @param mean.start Initial value for mu.
#' @param sd.start Initial value for sigma.
#' @param ci.level Number between 0.5 and 1. Gives a 100(ci.level)% confidence interval.
#' @param a Left truncation limit.
#' @param b Right truncation limit.
#' @param ... Parameters to pass to sampling()
#'
#' @importFrom stats median quantile coef
#' @export
#'
#' @example
#' trunc_est(x=c(-1,-2,1,2), mean.start=1, sd.start=0.5, ci.level=0.975, a=1, b=5)
#'
#' @references
#' https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html
trunc_est <- function(x,
                      mean.start = 0,
                      sd.start = 1,
                      ci.level = 0.95,
                      a,
                      b,
                      ...) {
  stopifnot(a < b)
  stopifnot(sd.start > 0)

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
  stan.model <- rstan::stan_model(model_code = model.text,
                                  isystem = "~/Desktop")


  post <- rstan::sampling(stan.model,
                          cores = 1,
                          #refresh = 0,
                          init = init.fcn, data=list(n=length(x), a=a, b=b, y=x), ...)


}, warning = function(condition){
  stan.warned <<- 1
  stan.warning <<- condition$message
} )

ext <- rstan::extract(post)
best.ind = which.max(ext$log_post)


postSumm <- summary(post)$summary

print(postSumm)
nlpost_simple = function(mean, sd) {
  nlpost.value = nlpost_jeffreys(mean, sd, x, a, b)
  return(nlpost.value)
}


res <- stats4::mle( minuslogl = nlpost_simple,
                    start = list(mean=ext$mu[best.ind], sd= ext$sigma[best.ind]), method="Nelder-Mead" )


maps <- as.numeric(stats4::coef(res))

# posterior means, then medians
mean.est <- median(rstan::extract(post, "mu")[[1]])
sd.est <- median(rstan::extract(post, "sigma")[[1]])

mean.maxlp <- ext$mu[best.ind]
sd.maxlp <- ext$sigma[best.ind]
# SEs
mean.se = postSumm["mu", "se_mean"]
sd.se = postSumm["sigma", "se_mean"]

# convert the numeric ci.level to strings of percentages.
l.lim.str <- paste0(toString((1-ci.level) * 100), "%")
r.lim.str <- paste0(toString((ci.level) * 100), "%")

# CI limits
sd.ci.lims = c( postSumm["sigma", l.lim.str], postSumm["sigma", r.lim.str] )
mean.ci.lims = c( postSumm["mu", l.lim.str], postSumm["mu", r.lim.str] )


mean.ci = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], 1-ci.level ),
                         quantile( rstan::extract(post, "mu")[[1]], ci.level ) ) )

sd.ci = as.numeric( c( quantile( rstan::extract(post, "sigma")[[1]], 1-ci.level ),
                       quantile( rstan::extract(post, "sigma")[[1]], ci.level ) ) )
# the point estimates are length 2 (post means, then medians),
#  but the inference is the same for each type of point estimate
return( list( post = post,
              mean.est = mean.est,
              sd.est = sd.est,

              mean.maxlp = mean.maxlp,
              sd.maxlp = sd.maxlp,
              mean.se = rep(mean.se, 2),
              sd.se= rep(sd.se, 2),

              mean.ci = mean.ci,
              sd.ci = sd.ci,

              stan.warned = stan.warned,
              stan.warning = stan.warning,

              mean.rhat = postSumm["mu", "Rhat"],
              sd.rhat = postSumm["sigma", "Rhat"]
) )
}

#' Finds the Fisher information matrix contained in n samples from a truncated normal distribution.

#' @param mean Mean of underlying normal distribution.
#' @param sd Standard deviation of underlying normal distribution.
#' @param n Number of observations.
#' @param a Lower truncation limit.
#' @param b Upper truncation limit.
#'
#' @importFrom stats dnorm pnorm
#' @export
E_fisher <- function(mean, sd, n, a, b) {
  stopifnot(sd > 0)
  stopifnot(a < b)

  Za = (a - mean) / sd
  Zb = (b - mean) / sd

  alpha.a = alpha_a(mean, sd, a, b)
  alpha.b = alpha_b(mean, sd, a, b)

  k11 = -(n/sd^2) + (n/sd^2)*( (alpha.b - alpha.a)^2 + (alpha.b*Zb - alpha.a*Za) )

  k12 = -( 2*n*(alpha.a - alpha.b) / sd^2 ) +
    (n/sd^2)*( alpha.a - alpha.b + alpha.b*Zb^2 - alpha.a*Za^2 +
                 (alpha.a - alpha.b)*(alpha.a*Za - alpha.b*Zb) )

  k22 = (n/sd^2) - (3*n*(1 + alpha.a*Za - alpha.b*Zb) / sd^2) +
    (n/sd^2)*( Zb*alpha.b*(Zb^2 - 2) - Za*alpha.a*(Za^2 - 2) +
                 (alpha.b*Zb - alpha.a*Za)^2 )

  return( matrix( c(-k11, -k12, -k12, -k22),
                  nrow = 2,
                  byrow = TRUE ) )
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
prior <- function(mean, sd, x, a, b) {
  stopifnot(a < b)
  stopifnot(sd > 0)

  return (log( sqrt( det( E_fisher(mean = mean, sd = sd, n = length(x), a = a, b = b) ) ) ))
}


