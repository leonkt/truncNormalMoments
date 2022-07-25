#' @keywords internal
#' @references
#' \insertRef{zhou2014}{truncnormbayes}
#'
#' \insertRef{stan2022}{truncnormbayes}
"_PACKAGE"

## usethis namespace: start
#' @useDynLib truncnormbayes, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom Rdpack reprompt
#' @importFrom stats dnorm pnorm median quantile
## usethis namespace: end
NULL
