#' Estimate truncated normal distribution
#'
#' Estimates the posterior modes for the mean (\code{mu}) and standard deviation (\code{sigma}) of the underlying normal
#' distribution, given truncated data with known truncation point(s).
#'
#' @param x Vector of observations from truncated normal
#' @param mu_start Initial value for mu.
#' @param sigma_start Initial value for sigma.
#' @param ci_level Number between 0.5 and 1. Gives a 100*2*(ci_level - 0.5)% symmetric HPD interval.
#' @param a Left truncation limit.
#' @param b Right truncation limit.
#' @param ... Parameters to pass to sampling()
#'
#' @return A list with two elements:
#'  \describe{
#'    \item{stats}{A data frame with two rows and the columns \code{param}
#'                 (\code{mu}, \code{sd}), \code{mu} (posterior mean),
#'                 \code{median} (posterior median), \code{maxlp}, \code{se},
#'                 \code{ci_lower}, \code{ci_upper}, \code{rhat}.}
#'    \item{fit}{A \code{stanfit} object (the result of fitting the model).}
#'  }
#'
#' @export
#'
#' @references
#' \insertRef{zhou2014}{truncnormbayes}
#' \insertRef{stan2022}{truncnormbayes}
#'
#' @examples
#' x <- truncnorm::rtruncnorm(100, a = -1, b = 2, mean = 0.5, sd = 0.5)
#' trunc_est(x, a = -1, b = 2)
trunc_est <- function(x,
                      a,
                      b,
                      mu_start = 0,
                      sigma_start = 1,
                      ci_level = 0.95,
                      ...) {
  stopifnot(a < b)
  stopifnot(sigma_start > 0)
  stopifnot(all(x >= a))
  stopifnot(all(x <= b))

  # set start values for sampler
  init_fcn <- function() list(mean = mu_start, sd = sigma_start)

  stan_fit <- rstan::sampling(stanmodels$trunc_est,
                              cores = 1,
                              init = init_fcn,
                              data = list(n = length(x), a = a, b = b, y = x),
                              ...)

  stan_extract <- rstan::extract(stan_fit)
  stan_summary <- as.data.frame(
    rstan::summary(stan_fit)$summary[c("mu", "sigma"), ]
  )
  means <- stan_summary$mean
  ses <- stan_summary$se_mean
  rhats <- stan_summary$Rhat

  medians <- c(median(stan_extract$mu), median(stan_extract$sigma))

  index_maxlp <- which.max(stan_extract$log_post)
  modes <- c(stan_extract$mu[index_maxlp], stan_extract$sigma[index_maxlp])

  stan_ci <- function(param, q) as.numeric(quantile(stan_extract[[param]], q))
  cil <- c(stan_ci("mu", 1 - ci_level), stan_ci("sigma", 1 - ci_level))
  ciu <- c(stan_ci("mu", ci_level), stan_ci("sigma", ci_level))

  stan_stats <- data.frame(param = c("mu", "sigma"), mode = modes, mean = means,
                           median = medians, se = ses, ci_lower = cil,
                           ci_upper = ciu, rhat = rhats)

  return(list(stats = stan_stats, fit = stan_fit))
}
