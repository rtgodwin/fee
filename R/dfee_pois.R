#' Marginal First-Exposure Effects for One-Inflated Poisson Model
#'
#' This function is used internally by \code{\link{dfee}}. The function computes
#' analytical derivatives of the FEE with respect to covariates, and uses discrete
#' differences for binary indicators.
#'
#' @param b A numeric vector of estimated coefficients for the rate.
#' @param g A numeric vector of estimated coefficients for inflation.
#' @param X A numeric matrix of covariates for the rate.
#' @param Z A numeric matrix of covariates for inflation.
#' @param dummies A character vector naming the binary (dummy) variables to compute marginal effects for.
#'
#' @return A matrix of marginal first-exposure effects. Each column corresponds
#' to a covariate, and each row corresponds to an observation in the data.
#'
#' @seealso \code{\link{dfee}}, \code{\link{fee_pois}}, \code{\link{dfee_nb}}
#'
#' @keywords internal

dfee_pois <- function(b, g, X, Z, dummies) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  L <- -l / (exp(l) - l - 1)
  w <- L + (1 - L) / (1 + t)
  
  dldq <- matrix(, nrow(X), (ncol(X) - 1))
  dzdq <- matrix(, nrow(Z), (ncol(Z) - 1))
  colnames(dldq) <- colnames(X)[-1]
  colnames(dzdq) <- colnames(Z)[-1]
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L
  dzdq[, colnames(Z)[-1]] <- -t %*% g[-1]
  dzdq[, -which(colnames(dzdq) %in% colnames(Z)[-1])] <- 0L
  
  dwdq <- dldq * as.vector(((exp(l) - l * exp(l) - 1) / (exp(l) - l - 1) ^ 2) * ((-t) / (1 + t))) - dzdq * as.vector(((exp(l) - 1) / ((exp(l) - l - 1) * (1 + t) ^ 2)))
  dfeedq <- dwdq * as.vector((1 - (l * exp(l)) / (exp(l) - 1))) + dldq * as.vector(((w * exp(l) * (l - exp(l) + 1)) / ((exp(l) - 1) ^ 2)))
  
  for(i in 1:length(dummies)) {
    Xd1 <- Xd0 <- X
    Zd1 <- Zd0 <- Z
    Xd1[ , dummies[i] == colnames(X)] <- 1
    Xd0[ , dummies[i] == colnames(X)] <- 0
    Zd1[ , dummies[i] == colnames(Z)] <- 1
    Zd0[ , dummies[i] == colnames(Z)] <- 0
    dfeedq[, dummies[i]] <- fee_pois(b, g, X=Xd1, Z=Zd1) - fee_pois(b, g, X=Xd0, Z=Zd0)
  }
  dfeedq
}