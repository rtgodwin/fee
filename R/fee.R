#' First-Exposure Effect
#'
#' Computes the first-exposure effect (FEE) from a fitted `oneinfl` model object.
#' The FEE measures the difference between the expected count for a first-time exposure
#' and the expected count under the baseline (non-inflated) model. The function supports
#' models estimated using either a one-inflated positive Poisson distribution or a one-inflated
#' zero-truncated negative binomial distribution.
#'
#' The effect can be evaluated in three ways, determined by the `at` argument:
#' \itemize{
#'   \item \code{"AE"}: Average the FEE over all data points (default).
#'   \item \code{"EM"}: Evaluate the FEE at the sample means of the covariates.
#'   \item \code{list}: Evaluate the FEE at a user-specified set of covariate values.
#' }
#'
#' If `at = "AE"`, the returned object also includes the total number of treatment visits
#' implied by the FEE across all observations.
#'
#' @param model A fitted model object of class \code{"oneinfl"}.
#' @param data The original data frame used to fit the model.
#' @param at A character string or list. Specifies how the first-exposure effect should be evaluated.
#'        Options are \code{"AE"} (average effect across the data), \code{"EM"} (effect at means), or
#'        a named list specifying covariate values for evaluating a representative case.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{fee}}{The estimated first-exposure effect.}
#'   \item{\code{sefee}}{The standard error of the estimated effect.}
#'   \item{\code{where}}{A character string describing the evaluation point.}
#'   \item{\code{treatment_visits}}{(Optional) Total implied treatment visits if \code{at = "AE"}.}
#' }
#'
#' @examples
#' # Example usage
#' df <- data.frame(x = runif(10,0,10), d = sample(c(0,1), 10, replace=TRUE), y = rpois(10, 3) + 1)
#' model <- oneinfl::oneinfl(formula = y ~ x + d | x + d, df = df, dist = "Poisson")
#' fee(model, data = df)
#'
#' @export

fee <- function(model, data, at="AE") {
  if (!inherits(model, "oneinflmodel")) {
    stop("Error: 'model' must be of class 'oneinflmodel'.")
  }
  
  q <- list()
  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") {a <- model$alpha}
  formula <- model$formula
  cleandata <- oneinfl::makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z
  if (is.list(at)) {
    q$where <- "First exposure effect evaluated at: "
    for(i in 1:(length(at) - 1)) {
      q$where <- paste(q$where, names(at[i]), " = ", unlist(at[i]), ", ", sep = "")
    }
    q$where <- paste(q$where, names(at[length(at)]), " = ", unlist(at[length(at)]), sep = "")
    if (!all(names(at) %in% names(data))) {stop("variable names in 'at' must match those in the data")}
    for(i in 1:length(at)) {
      if(names(at[i]) %in% colnames(X)) {X <- '[<-'(X, , names(at[i]), unlist(at[i]))}
      if(names(at[i]) %in% colnames(Z)) {Z <- '[<-'(Z, , names(at[i]), unlist(at[i]))}
    }
  } else if (at == "EM") {
    q$where <- "First exposure effect evaluated at the sample means of the data"
    X <- colSums(X) / nrow(X)
    Z <- colSums(Z) / nrow(Z)
  } else if (at == "AE") {q$where <- "First exposure effect averaged over all data points"
  } else {
    stop("'at' must be 'AE' (average effect), 'EM' (effect at means), or a list of representative cases in which to evaluate the effect")
  }
  if (model$dist == "Poisson") {
    q$fee <- mean(fee_pois(b, g, X, Z))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(fee_pois(b, g, X, Z)), c("b", "g")), "gradient")))
    q$sefee <- sqrt(t(J) %*% model$vc %*% J)
    if (at == "AE") {q$treatment_visits <- sum(fee_pois(b, g, X, Z))}
  } else if (model$dist == "negbin") {
    q$fee <- mean(fee_negbin(b, g, a, X, Z))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(fee_negbin(b, g, a, X, Z)), c("b", "g", "a")), "gradient")))
    q$sefee <- sqrt(t(J) %*% model$vc %*% J)
    if (at == "AE") {q$treatment_visits <- sum(fee_negbin(b, g, a, X, Z))}
  }
  q
}