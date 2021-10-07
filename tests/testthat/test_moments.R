library(testthat)
library(rstan)

n <- 100
m <- 1.2
s <- 8
ul <- 1.2
uh <- 2
x <- 2

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

test_that("Cumulants evaluated at point matches analytical expression.", {
  expect_equal(get_cumulants("mm", m, s, ul, uh, n), -(n/s**2) + (n/s**2) * ((alphah(m, s, ul, uh) - alphal(m, s, ul, uh)) ** 2 + alphah(m, s, ul, uh) * ush(m, s, ul, uh) - alphal(m, s, ul, uh) * usl(m, s, ul, uh)))
  expect_equal(get_cumulants("ms", m, s, ul, uh, n), (-2 * n/s**3) * Ex_min_m(m,s,ul,uh) + (n/s**2) * (alphal(m, s, ul, uh) - alphah(m, s, ul, uh) + alphah(m, s, ul, uh) * ush(m, s, ul, uh) **2 - alphal(m, s, ul, uh) * usl(m, s, ul, uh)**2 + (alphal(m, s, ul, uh) - alphah(m, s, ul, uh)) * (alphal(m, s, ul, uh)*usl(m, s, ul, uh) - alphah(m, s, ul, uh) * ush(m, s, ul, uh)) ))
  expect_equal(get_cumulants("ss", m, s, ul, uh, n), (n/s**2) - (3 * n/s**4) * Ex_min_m_sq(m, s, ul, uh) + n * ((ush(m, s, ul, uh) * alphah(m, s, ul, uh)/s**2 ) * (ush(m, s, ul, uh) **2 - 2) - (usl(m, s, ul, uh)  * alphal(m, s, ul, uh) / s**2) * (usl(m, s, ul, uh) **2 - 2) + ((alphah(m, s, ul, uh) * ush(m, s, ul, uh) - alphal(m, s, ul, uh) * usl(m, s, ul, uh) )**2)/s**2))
})

test_that("Fisher Information Matrix matches K.", {
  expect_equal(find_K(m, s, ul, uh, n), E_fisher(m, s, n, ul, uh))
})

test_that("Full nlpost is simple plus Jeffreys.", {
  expect_equal(neg_log_post(c(m,s), x, ul, uh ), prior(c(m,s), x, ul, uh) + nlpost_simple(m,s,"var", x,ul,uh))
})

