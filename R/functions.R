library(testthat)

################################ INTERMEDIATE TERMS ################################

#' Returns the standardized lower truncation limit for a normal distribution.
#'
#' @param m Mean of normal distribution.
#' @param s Standard deviation of normal distribution.
#' @param ul Unstandardized lower truncation limit.
#' @examples
#' Returns the z-score of the value 1, for a normal distribution of mean 0, standard deviation 2.
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
#' @examples
#' Returns the z-score of the value 1, for a normal distribution of mean 0, standard deviation 2.
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
#' d_str should consist of the characters "m" and "s" indicating what types of partial derivatives to take.
#' The length of d_str indicates the order of the derivative.
#' For example, the following calculates the third-order derivative of log-likelihood with respect to mu, then sigma, then sigma.
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
#' d_str should consist of the characters "m" and "s" indicating what types of partial derivatives to take.
#' The length of d_str indicates the order of the derivative.
#' For example, the following calculates the third-order derivative of log-likelihood with respect to mu, then sigma, then sigma.
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
#' In each case, we make a function f. f is the sum of 
#' (1) get_deriv(d_str), which gives an expression for the appropriate derivative of the LOG LIKELIHOOD WITHOUT THE X TERM.
#' (2) expectation of the term involving x, in the appropriate derivative of the FULL LOG LIKELIHOOD.
#' 
#' For example:
#' d^3\sigma/d\sigma^3 = (some terms not dependent on x) - (3/\sigma^4) * (\sum(x_i - \mu)^2).
#' 
#' Therefore, taking expectations will keep the terms not dependent on x the same, while changing the latter term dependent on x
#' to -(3n/\sigma^4) * (\sigma(1 + usl * alphal - ush * alphah)). This sum is our cumulant, which we denote as f in the function. 
#'
#' We then take the derivative of f w.r.t var, since f is equal to the appropriate cumulant.

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

#' Ensures that the cumulants calculated in the Deriv package match those calculated by hand.
#' 
#' @param m Mean of underlying normal distribution.
#' @param s Standard deviation of underlying normal distribution.
#' @param ul Lower truncation limit.
#' @param uh Upper truncation limit.
#' @param n Number of observations.

cumulants_check <- function(m, s, ul, uh, n) {
  epsilon = 0.000001
  # Calculates E(x-mu)
  Ex_min_m <- function( m, s, ul, uh) {
    (s * (alphal( m, s, ul, uh) - alphah( m, s, ul, uh)) )
  }
  # Calculates E(x-mu)^2
  Ex_min_m_sq <- function(m, s, ul, uh) {
    (s**2 * (1 + usl(m, s, ul, uh) * alphal(m, s, ul, uh) - ush(m, s, ul, uh) * alphah(m, s, ul, uh)) )
  }
  # Theoretical derivation of cumulants
  d_mm = -(n/s**2) + (n/s**2) * ((alphah(m, s, ul, uh) - alphal(m, s, ul, uh)) ** 2 + alphah(m, s, ul, uh) * ush(m, s, ul, uh) - alphal(m, s, ul, uh) * usl(m, s, ul, uh))
  
  d_sm =  (-2 * n/s**3) * Ex_min_m(m,s,ul,uh) + (n/s**2) * (alphal(m, s, ul, uh) - alphah(m, s, ul, uh) + alphah(m, s, ul, uh) * ush(m, s, ul, uh) **2 - alphal(m, s, ul, uh) * usl(m, s, ul, uh)**2 + (alphal(m, s, ul, uh) - alphah(m, s, ul, uh)) * (alphal(m, s, ul, uh)*usl(m, s, ul, uh) - alphah(m, s, ul, uh) * ush(m, s, ul, uh)) )
  
  d_ss = (n/s**2) - (3 * n/s**4) * Ex_min_m_sq(m, s, ul, uh) + n * ((ush(m, s, ul, uh) * alphah(m, s, ul, uh)/s**2 ) * (ush(m, s, ul, uh) **2 - 2) - (usl(m, s, ul, uh)  * alphal(m, s, ul, uh) / s**2) * (usl(m, s, ul, uh) **2 - 2) + ((alphah(m, s, ul, uh) * ush(m, s, ul, uh) - alphal(m, s, ul, uh) * usl(m, s, ul, uh) )**2)/s**2)
  
  # Ensures that the value of the cumulants we get from get_cumulants match the derived value of the cumulants.
  assert("mm check:", get_cumulants("mm", m, s, ul, uh, n) -  epsilon <= d_mm && d_mm <= get_cumulants("mm", m, s, ul, uh, n) + epsilon)
  assert("ss check:", get_cumulants("ss", m, s, ul, uh, n) -  epsilon <= d_ss && d_ss <= get_cumulants("ss", m, s, ul, uh, n) + epsilon)
  assert("sm check:", get_cumulants("sm", m, s, ul, uh, n) -  epsilon <= d_sm && d_sm <= get_cumulants("sm", m, s, ul, uh, n) + epsilon)
  
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
  K <- find_K( m=.m, s=.s, ul=.ul, uh=.uh, n=.n)
  A <- find_A( m=.m, s=.s, ul=.ul, uh=.uh, n=.n)
  bias <- inv(K) %*% A %*% c(inv(K))
  return( list(bias = bias, K = K, A = A) )
}

