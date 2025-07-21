#' First-Exposure Effect for One-Inflated Positive Poisson Model
#'
#' Computes the first-exposure effect (FEE) for a model estimated using the
#' one-inflated positive Poisson distribution. The FEE represents the change in 
#' expected counts associated with the first exposure to treatment, relative to 
#' the counterfactual positive Poisson model.This function is used internally by
#' \code{\link{fee}} at given covariate values.
#'
#' @param b A numeric vector of estimated coefficients for the rate.
#' @param g A numeric vector of estimated coefficients for inflation.
#' @param X A numeric matrix of covariates for the rate.
#' @param Z A numeric matrix of covariates for inflation.
#'
#' @return A numeric vector of first-exposure effect(s).
#'
#' @keywords internal
#' @seealso \code{\link{fee}}, \code{\link{fee_negbin}}

fee_pois <- function(b, g, X, Z) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
  w * (1 - l * exp(l) / (exp(l) - 1))
}