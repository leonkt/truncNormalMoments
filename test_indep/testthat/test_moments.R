library(testthat)
library(rstan)

######## DEFINING VARIABLES #########
n <- 20
m <- 1.2
s <- 8
ul <- 0
uh <- 3
x <- 2
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
l.lim <- 0.025
r.lim <- 0.975

stan.model <- stan_model(model_code = model.text,
                         isystem = "~/Desktop")
post <- sampling(stan.model, data = list( n = n, LL = ul, UU = uh, y =c(0.385258397941442, 1.68066267127739, 0.729227742032434, 0.479120432291688,
                                                   0.897279068695914, 0.0575356881970433, 0.165652783807015, 0.875647820475464,
                                                   0.380168104717168, 0.825551957468494, 0.589842597791253, 0.402854395794205,
                                                   1.7668857263465, 0.649703054651576, 0.83074621395, 0.407612065235468,
                                                   1.49884534475478, 0.420770708293162, 0.699931456061883, 1.06188921169992)),  iter=100)
print(post)
print(summary(post))
postSumm <- summary(post)$summary

myMhatCI <- as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], l.lim ),
                          quantile( rstan::extract(post, "mu")[[1]], r.lim ) ) )
M.CI <- c( postSumm["mu", "2.5%"], postSumm["mu", "97.5%"] )
Mhat <- c( postSumm["mu", "mean"], median( rstan::extract(post, "mu")[[1]] ) )
MhatSE <- postSumm["mu", "se_mean"]

######## AUXILLIARY FUNCTIONS #########

# Auxilliary function. Defined within function, so testthat cannot call the function.
# Calculates E(x-mu)
Ex_min_m <- function( m, s, ul, uh) {
  (s * (alphal( m, s, ul, uh) - alphah( m, s, ul, uh)) )
}

# Auxilliary function. Defined within function, so testthat cannot call the function.
# Calculates E(x-mu)^2
Ex_min_m_sq <- function(m, s, ul, uh) {
  (s**2 * (1 + usl(m, s, ul, uh) * alphal(m, s, ul, uh) - ush(m, s, ul, uh) * alphah(m, s, ul, uh)) )
}

######## TESTS #########

test_that("Cumulants evaluated at point matches analytical expression.", {
  expect_equal(get_cumulants("mm", m, s, ul, uh, n), -(n/s**2) + (n/s**2) * ((alphah(m, s, ul, uh) - alphal(m, s, ul, uh)) ** 2 + alphah(m, s, ul, uh) * ush(m, s, ul, uh) - alphal(m, s, ul, uh) * usl(m, s, ul, uh)))
  expect_equal(get_cumulants("ms", m, s, ul, uh, n), (-2 * n/s**3) * Ex_min_m(m,s,ul,uh) + (n/s**2) * (alphal(m, s, ul, uh) - alphah(m, s, ul, uh) + alphah(m, s, ul, uh) * ush(m, s, ul, uh) **2 - alphal(m, s, ul, uh) * usl(m, s, ul, uh)**2 + (alphal(m, s, ul, uh) - alphah(m, s, ul, uh)) * (alphal(m, s, ul, uh)*usl(m, s, ul, uh) - alphah(m, s, ul, uh) * ush(m, s, ul, uh)) ))
  expect_equal(get_cumulants("ss", m, s, ul, uh, n), (n/s**2) - (3 * n/s**4) * Ex_min_m_sq(m, s, ul, uh) + n * ((ush(m, s, ul, uh) * alphah(m, s, ul, uh)/s**2 ) * (ush(m, s, ul, uh) **2 - 2) - (usl(m, s, ul, uh)  * alphal(m, s, ul, uh) / s**2) * (usl(m, s, ul, uh) **2 - 2) + ((alphah(m, s, ul, uh) * ush(m, s, ul, uh) - alphal(m, s, ul, uh) * usl(m, s, ul, uh) )**2)/s**2))
})

test_that("Fisher Information Matrix matches K.", {
  expect_equal(find_K(m, s, ul, uh, n), E_fisher(m, s, n, ul, uh))
})

test_that("Confidence Interval for Mhat matches", {
  expect_equal(M.CI, myMhatCI)
})

test_that("Parameter and SE estimates from the posterior and from extract() match.", {
  expect_equal(Mhat[1], mean( rstan::extract(post, "mu")[[1]] ) )
  expect_equal( postSumm["mu", "sd"], sd( rstan::extract(post, "mu")[[1]] ) )
  expect_equal( MhatSE,postSumm["mu", "sd"] / sqrt( postSumm["mu", "n_eff"] ) )
})



