#' @keywords internal
#' @references
#' \insertRef{zhou2014}{truncnormbayes}
#'
#' \insertRef{stan2022}{truncnormbayes}
"_PACKAGE"

#' @useDynLib truncnormbayes, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom Rdpack reprompt
#' @importFrom stats dnorm pnorm median quantile
NULL
