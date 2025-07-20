#' Marginal First-Exposure Effects for One-Inflated Zero-Truncated Negative Binomial model 
#'
#' This function is used internally by \code{\link{dfee}}. The function computes
#' analytical derivatives of the FEE with respect to covariates, and uses discrete
#' differences for binary indicators.
#'
#' @param b A numeric vector of estimated coefficients for the rate.
#' @param g A numeric vector of estimated coefficients for inflation.
#' @param a A dispersion parameter from the negative binomial model.
#' @param X A numeric matrix of covariates for the rate.
#' @param Z A numeric matrix of covariates for inflation.
#' @param dummies A character vector naming the binary (dummy) variables to compute marginal effects for.
#'
#' @return A matrix of marginal first-exposure effects. Each column corresponds
#' to a covariate, and each row corresponds to an observation in the data.
#'
#' @seealso \code{\link{dfee}}, \code{\link{fee_negbin}}, \code{\link{dfee_pois}}
#'
#' @keywords internal

dfee_nb <- function(b, g, a, X, Z, dummies) {
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  th <- l / a
  P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
  L <- -P1 / (1 - P1)
  w <- L + (1 - L) / (1 + t)
  
  dldq <- matrix(, nrow(X), (ncol(X) - 1))
  dzdq <- matrix(, nrow(Z), (ncol(Z) - 1))
  colnames(dldq) <- colnames(X)[-1]
  colnames(dzdq) <- colnames(Z)[-1]
  dldq[, colnames(X)[-1]] <- l %*% b[-1]
  dldq[, -which(colnames(dldq) %in% colnames(X)[-1])] <- 0L
  dzdq[, colnames(Z)[-1]] <- -t %*% g[-1]
  dzdq[, -which(colnames(dzdq) %in% colnames(Z)[-1])] <- 0L
  
  dfdq <- dldq * as.vector(((1 + l / a) ^ (-a)) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (a ^ 2) * l / ((a + l) ^ 2) - (l / a) * (1 + l / a - (1 + l / a) ^ (1 - a)) ^ (-1) * (1 - (1 - a) * (1 + l / a) ^ (-a))))
  dLdq <- -dfdq / as.vector(((1 - P1) ^ 2))
  dwdq <- dLdq * as.vector((1 - 1 / (1 + t))) - dzdq * as.vector(((1 - P1) / (1 + t) ^ 2))
  dfeedq <- dwdq * as.vector((1 - l / (1 - (1 + l / a) ^ (-a)))) - dldq * as.vector((w / (1 - (1 + l / a) ^ (-a))) * (1 + (l * (1 + l / a) ^ (-a - 1)) / (1 - (1 + l / a) ^ (-a))))
  
  for(i in 1:length(dummies)) {
    Xd1 <- Xd0 <- X
    Zd1 <- Zd0 <- Z
    Xd1[ , dummies[i] == colnames(X)] <- 1
    Xd0[ , dummies[i] == colnames(X)] <- 0
    Zd1[ , dummies[i] == colnames(Z)] <- 1
    Zd0[ , dummies[i] == colnames(Z)] <- 0
    dfeedq[, dummies[i]] <- fee_negbin(b, g, a, X=Xd1, Z=Zd1) - fee_negbin(b, g, a, X=Xd0, Z=Zd0)
  }
  dfeedq
}