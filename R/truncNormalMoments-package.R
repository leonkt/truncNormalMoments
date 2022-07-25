#' @keywords internal
#' @references
#' \insertRef{zhou2014}{truncNormalMoments}
#'
#' \insertRef{stan2022}{truncNormalMoments}
"_PACKAGE"

## usethis namespace: start
#' @useDynLib truncNormalMoments, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom Rdpack reprompt
#' @importFrom stats dnorm pnorm median quantile
## usethis namespace: end
NULL
