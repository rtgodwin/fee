#' Predicted Counts from One-Inflated or Truncated Count Models
#'
#' Computes the predicted count distribution from a fitted model of class \code{"oneinflmodel"},
#' \code{"truncmodel"}, or \code{"basicPoisson"}. The function returns the expected frequency
#' for each count value from 1 up to \code{maxpred}, based on the model's parameters.
#'
#' @param model A fitted model object of class \code{"oneinflmodel"}, \code{"truncmodel"}, or \code{"basicPoisson"}.
#' @param data A data frame containing the covariates used to fit the model.
#' @param maxpred Optional integer specifying the maximum count value for which to compute predicted frequencies.
#'        If not supplied, defaults to the maximum observed count in the data.
#'
#' @return A numeric vector of length \code{maxpred}, giving the predicted expected frequency of each count from 1 to \code{maxpred}.
#'
#' @details
#' The function determines the model type based on its class and the \code{dist} attribute, and applies the appropriate
#' density function:
#' \itemize{
#'   \item For \code{oneinflmodel} (Poisson): one-inflated positive Poisson distribution.
#'   \item For \code{oneinflmodel} (negbin): one-inflated zero-truncated negative binomial.
#'   \item For \code{truncmodel} (Poisson): truncated positive Poisson.
#'   \item For \code{truncmodel} (negbin): zero-truncated negative binomial.
#' }
#'
#' @seealso \code{\link{feeplot}}, \code{\link{fee}}, \code{\link{dfee}}
#'
#' @examples
#' df <- data.frame(x = runif(10,0,10), d = sample(c(0,1), 10, replace=TRUE), y = rpois(10, 3) + 1)
#' model <- oneinfl::oneinfl(formula = y ~ x + d | x + d, df = df, dist = "Poisson")
#' feepred(model, data = df)
#'
#' @export

feepred <- function(model, data, maxpred) {
  
  b <- model$beta
  g <- model$gamma
  
  formula <- model$formula
  cleandata <- oneinfl::makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  y <- cleandata$y
  
  l <- exp(X %*% b)
  t <- exp(-Z %*% g)
  
  if (model$dist == "negbin") {
    a <- model$alpha
    th <- l / a
  }
  
  if (missing(maxpred)) {
    maxpred <- max(y)
  }
  
  if (inherits(model, "truncmodel") && model$dist == "Poisson") {
    pred <- numeric(maxpred)
    for (pp in 1:maxpred) {
      pred[pp] <- sum((l ^ pp) / ((exp(l) - 1) * factorial(pp)))
    }
    return(pred)
    
  } else if (inherits(model, "truncmodel") && model$dist == "negbin") {
    pred <- numeric(maxpred)
    for (pp in 1:maxpred) {
      pred[pp] <- sum((gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
    return(pred)
    
  } else if (inherits(model, "oneinflmodel") && model$dist == "Poisson") {
    pred <- numeric(maxpred)
    w <- -l * (exp(l) - l - 1) ^ -1 + (1 + l * (exp(l) - l - 1) ^ -1) * (1 + t) ^ -1
    pred[1] <- sum(w + (1 - w) * l / (exp(l) - 1))
    for (pp in 2:maxpred) {
      pred[pp] <- sum((1 - w) * (l^pp) / ((exp(l) - 1) * factorial(pp)))
    }
    return(pred)
    
  } else if (inherits(model, "oneinflmodel") && model$dist == "negbin") {
    pred <- numeric(maxpred)
    P1 <- a * ((1 / (1 + th)) ^ a) * th / (1 + th - (1 + th) ^ (1 - a))
    L <- -P1 / (1 - P1)
    w <- L + (1 - L) / (1 + t)
    pred[1] <- sum(w + (1 - w) * a * ((1 / (1 + th)) ^ a) * (th / (1 + th - (1 + th) ^ (1 - a))))
    for (pp in 2:maxpred) {
      pred[pp] <- sum((1 - w) * (gamma(a + pp) / gamma(a) / gamma(pp + 1)) * ((1 / (1 + th)) ^ a) * ((th / (1 + th)) ^ pp) * (1 / (1 - (1 + th) ^ (-a))))
    }
    return(pred)
  }
}