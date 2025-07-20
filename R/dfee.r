#' Marginal First-Exposure Effects
#'
#' Computes marginal first-exposure effects from a fitted \code{oneinfl} model.
#' Dummy variables are automatically detected as those with exactly two unique values
#' in the data, and corresponding marginal effects are instead calculated by differencing
#' the FEE between both values of the dummy.
#' 
#' The marginal effects can be evaluated in three ways, determined by the \code{at} argument:
#' \itemize{
#'   \item \code{"AE"}: Average over all data points (default).
#'   \item \code{"EM"}: Evaluate at the sample means of the covariates.
#'   \item \code{list}: Evaluate at a user-specified set of covariate values.
#' }
#'
#' @param model A fitted model object of class \code{"oneinfl"}.
#' @param data A data frame containing the variables used to fit the model.
#' @param at A character string or list. Specifies where the marginal FEE should be evaluated.
#'        Options are \code{"AE"} (average), \code{"EM"} (means), or a named list of covariate values.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{dfee}}{A named numeric vector of estimated marginal first-exposure effects for each variable.}
#'   \item{\code{sefee}}{A numeric vector of standard errors corresponding to the marginal effects.}
#'   \item{\code{where}}{A character string describing the evaluation point.}
#' }
#'
#' @seealso \code{\link{fee}}, \code{\link{dfee_pois}}, \code{\link{dfee_nb}}
#'
#' @examples
#' df <- data.frame(x = runif(10,0,10), d = sample(c(0,1), 10, replace=TRUE), y = rpois(10, 3) + 1)
#' model <- oneinfl::oneinfl(formula = y ~ x + d | x + d, df = df, dist = "Poisson")
#' dfee(model, data = df)
#'
#' @export

dfee <- function(model, data, at="AE") {
  
  q <- list()

  b <- model$beta
  g <- model$gamma
  if (model$dist == "negbin") {a <- model$alpha}

  names(b) <- substring(names(b), 2)
  names(g) <- substring(names(g), 2)

  formula <- model$formula
  cleandata <- oneinfl::makeXZy(formula, data)
  X <- cleandata$X
  Z <- cleandata$Z

  is.dummy <- function(X) {length(unique(X)) == 2}
  dummies <- colnames(data[-1])[apply(data[-1], 2, is.dummy)]

  if (is.list(at)) {
    q$where <- "Marginal first exposure effect evaluated at: "
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
    q$where <- "Marginal first exposure effect evaluated at the sample means of the data"
    X <- matrix((colSums(X) / nrow(X)), 1, ncol(X), dimnames = list(1, colnames(X)))
    Z <- matrix((colSums(Z) / nrow(Z)), 1, ncol(Z), dimnames = list(1, colnames(Z)))
  } else if (at == "AE") {q$where <- "Marginal first exposure effect averaged over all data points"
  } else {
    stop("'at' must be 'AE' (average effect), 'EM' (effect at means), or a list of representative cases in which to evaluate the effect")
  }

  if (model$dist == "Poisson") {
    q$dfee <- colMeans(dfee_pois(b, g, X, Z, dummies))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(dfee_pois(b, g, X, Z, dummies)), c("b", "g")), "gradient")))
    q$sefee <- sqrt(diag(J %*% model$vc %*% t(J)))
  }

  if (model$dist == "negbin") {
    q$dfee <- colMeans(dfee_nb(b, g, a, X, Z, dummies))
    J <- as.matrix(colMeans(attr(numericDeriv(quote(dfee_nb(b, g, a, X, Z, dummies)), c("b", "g", "a")), "gradient")))
    q$sefee <- sqrt(diag(J %*% model$vc %*% t(J)))
  }
  q
}