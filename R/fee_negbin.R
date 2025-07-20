#' First-Exposure Effect for One-Inflated Zero-Truncated Negative Binomial Model
#'
#' Computes the first-exposure effect (FEE) for a model estimated using the
#' one-inflated zero-truncated negative binomial distribution. The FEE represents
#' the change in expected counts associated with the first exposure to treatment,
#' relative to the counterfactual zero-truncated negative binomial model. This
#' function is used internally by \code{\link{fee}} at given covariate values.
#'
#' @param b A numeric vector of estimated coefficients for the rate.
#' @param g A numeric vector of estimated coefficients for inflation.
#' @param a A dispersion parameter from the negative binomial model.
#' @param X A numeric matrix of covariates for the rate.
#' @param Z A numeric matrix of covariates for inflation.
#'
#' @return A numeric vector of first-exposure effect(s).
#'
#' @keywords internal
#' @seealso \code{\link{fee}}, \code{\link{fee_pois}}

fee_negbin <- function(b, g, a, X, Z) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  th <- l / a
  P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
  L <- -P1 / (1 - P1)
  w <- L + (1 - L) / (1 + t)
  w * (1 - l * (1 - (1 + th) ^ -a) ^ -1)
}