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
#' @export
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
#' @export
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

#' Finds the Fisher information matrix contained in n samples from a truncated
#' normal distribution.
#'
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
