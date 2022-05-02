gc()

# test issue:
# - cannot use library(rstan) here
# - NAMESPACE needs to use import(rstan), not just importFrom(rstan, ...)

# https://stackoverflow.com/questions/69919465/r-package-development-tests-pass-in-console-but-fail-via-devtoolstest#69922844


######## DEFINING VARIABLES #########

# compare fitting model manually vs. with trunc_est

n <- 1000
m <-2
s <- 0.5
a.alt <- -1
b.alt <- 5
x.alt <- rtruncnorm(500, a=a.alt, b=b.alt, mean=m, sd=s)
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



# need to have isystem arg to avoid "syntax error"
# but either way, it says parsing has been successful
stan.model <- stan_model(model_code = model.text, verbose = TRUE)

post.alt <- sampling(stan.model, data = list( n = length(x.alt),
                                              a = a.alt, b = b.alt,
                                              y = x.alt, iters=100))

postSumm.alt <- summary(post.alt)$summary
meanSE.alt <- postSumm.alt["mu", "se_mean"]
meanCI.75 <- c( postSumm.alt["mu", "25%"], postSumm.alt["mu", "75%"] )

res.alt <- trunc_est(
  x = x.alt,
  mean.start = m,
  sd.start = s,
  a = a.alt,
  b = b.alt,
  ci.level=0.975)



##### TESTS #########


test_that("trunc_est approx. matches STAN MAPs.", {
  expect_lt(abs(meanSE.alt[1]- res.alt$mean.se[1]),  2)
  expect_lt(abs(meanSE.alt[1]- res.alt$mean.se[1]),  2)
})


gc()
