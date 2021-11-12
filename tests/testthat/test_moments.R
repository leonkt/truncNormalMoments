

# test issue:
# - cannot use library(rstan) here
# - NAMESPACE needs to use import(rstan), not just importFrom(rstan, ...)

# https://stackoverflow.com/questions/69919465/r-package-development-tests-pass-in-console-but-fail-via-devtoolstest#69922844


######## DEFINING VARIABLES #########

# compare fitting model manually vs. with trunc_est

n <- 20
m <-0
s <- 1
a <- -3
b <- 3
x <- c(0.385258397941442, 1.68066267127739, 0.729227742032434, 0.479120432291688,
       0.897279068695914, 0.0575356881970433, 0.165652783807015, 0.875647820475464,
       0.380168104717168, 0.825551957468494, 0.589842597791253, 0.402854395794205,
       1.7668857263465, 0.649703054651576, 0.83074621395, 0.407612065235468,
       1.49884534475478, 0.420770708293162, 0.699931456061883, 1.06188921169992)

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
l.lim <- 0.025
r.lim <- 0.975


# need to have isystem arg to avoid "syntax error"
# but either way, it says parsing has been successful
stan.model <- stan_model(model_code = model.text, verbose = TRUE)
post <- sampling(stan.model, data = list( n = length(x), a = a, b = b, y = x, iters=100))


postSumm <- summary(post)$summary
M.CI <- c( postSumm["mu", "2.5%"], postSumm["mu", "97.5%"] )
# Mhat <- c( postSumm["mu", "mean"], median( rstan::extract(post, "mu")[[1]] ) )
MhatSE <- postSumm["mu", "se_mean"]

# now with trunc_est: to be filled in
#MM: the tests need to compare the manual fit to the package, not the manual fit to itself
# remember also to test the other arguments, like ci.level, a, b, by trying differently values, as we discussed earlier
# remember also that we need to test that the functions generate the appropriate warnings when users pass bad input, etc.
res = trunc_est(
  x = x,
  mean.start = m,
  sd.start = s,
  a = a,
  b = b )



######## TESTS #########





