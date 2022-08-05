test_that("trunc_est inputs", {

  x <- truncnorm::rtruncnorm(100, a = 0, b = 2, mean = 0.5, sd = 0.5)

  # a must be less than b
  expect_error(trunc_est(x, a = 4, b = 2))

  # sigma_start must be greater than 0
  expect_error(trunc_est(x, a = 4, b = 2, sigma_start = 0))
  expect_error(trunc_est(x, a = 4, b = 2, sigma_start = -1))

  # x values must be between a and b
  expect_error(trunc_est(x, a = 1, b = 2))
  expect_error(trunc_est(x, a = 0, b = 1))

})

check_estimates <- function(mu, sigma, a, b, tolerance = 0.1) {
  x <- truncnorm::rtruncnorm(1000, a = a, b = b, mean = mu, sd = sigma)
  est <- trunc_est(x, a = a, b = b, mu_start = mu, sigma_start = sigma)
  modes <- est$stats$mode
  message(sprintf("estimated: mu = %.2f, sigma = %.2f", modes[1], modes[2]))
  message(sprintf("actual: mu = %.2f, sigma = %.2f", mu, sigma))
  expect_equal(modes[1], mu, tolerance = tolerance)
  expect_equal(modes[2], sigma, tolerance = tolerance)
}

test_that("trunc_est [untruncated] [standard normal]", {
  check_estimates(mu = 0, sigma = 1, a = -10, b = 10)
})

test_that("trunc_est [untruncated] [non-standard normal]", {
  check_estimates(mu = 1, sigma = 0.5, a = -10, b = 10)
})

test_that("trunc_est [truncated symmetrically] [standard normal]", {
  check_estimates(mu = 0, sigma = 1, a = -2, b = 2)
})

test_that("trunc_est [truncated symmetrically] [non-standard normal]", {
  check_estimates(mu = 1, sigma = 0.5, a = -2, b = 2)
})

test_that("trunc_est [truncated left] [standard normal]", {
  check_estimates(mu = 0, sigma = 1, a = -2, b = 10)
})

test_that("trunc_est [truncated right] [standard normal]", {
  check_estimates(mu = 0, sigma = 1, a = -10, b = 2)
})
