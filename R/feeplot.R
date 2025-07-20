#' Plot First-Exposure Effect Distribution
#'
#' Creates a bar plot of the observed count data overlaid with fitted values from
#' a \code{oneinfl} model and its associated counterfactual (non-inflated) model.
#'
#' The factual predictions come from the fitted \code{oneinfl} model, while the 
#' counterfactual distribution is obtained by transforming the model into a non-inflated
#' counterpart using the same estimated parameters.
#'
#' @param model A fitted object of class \code{"oneinflmodel"}.
#' @param data A data frame containing the original data used to fit the model.
#' @param maxpred Optional integer indicating the maximum count value to include on the x-axis.
#'        Defaults to the maximum observed value in \code{data}.
#' @param ylimit Optional numeric value specifying the upper limit of the y-axis.
#'        Defaults to 110\% of the maximum observed count frequency.
#' @param cex A numeric value for point size in the overlay plot. Defaults to \code{1.5}.
#' @param lwd A numeric value for line width in the overlay plot. Defaults to \code{1.5}.
#' @param ... Additional arguments passed to \code{points()} or \code{lines()}.
#'
#' @return A barplot with overlaid predicted values from the factual and counterfactual
#' distributions, including a legend identifying each component.
#'
#' @details
#' The function detects whether the model is based on a Poisson or negative binomial
#' distribution and selects the appropriate counterfactual model.
#'
#' @seealso \code{\link{feepred}}, \code{\link{fee}}, \code{\link{dfee}}
#'
#' @examples
#' df <- data.frame(x = runif(10,0,10), d = sample(c(0,1), 10, replace=TRUE), y = rpois(10, 3) + 1)
#' model <- oneinfl::oneinfl(formula = y ~ x + d | x + d, df = df, dist = "Poisson")
#' feeplot(model, data = df)
#'
#' @export

feeplot <- function(model, data, maxpred, ylimit, cex = 1.5, lwd = 1.5, ...) {
  
  plotpp <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch = 25, col = "darkmagenta", cex = cex)
    lines(x = df.bar[,1], y = preds, col = "darkmagenta", lwd = lwd)
    list(label = "counterfactual distribution", color = "darkmagenta", pch = 25)
  }
  
  plotztnb <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch = 18, col = "red", cex = cex)
    lines(x = df.bar[,1], y = preds, col = "red", lwd = lwd)
    list(label = "counterfactual distribution", color = "red", pch = 18)
  }
  
  plotoipp <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch = 17, col = "green", cex = cex)
    lines(x = df.bar[,1], y = preds, col = "green", lwd = lwd)
    list(label = "factual distribution", color = "green", pch = 17)
  }
  
  plotoiztnb <- function(model, data, maxpred) {
    preds <- feepred(model, data, maxpred)
    points(x = df.bar[,1], y = preds, pch = 16, col = "blue", cex = cex)
    lines(x = df.bar[,1], y = preds, col = "blue", lwd = lwd)
    list(label = "factual distribution", color = "blue", pch = 16)
  }
  
  formula <- model$formula
  cleandata <- oneinfl::makeXZy(formula, data)
  y <- cleandata$y
  
  if (missing(maxpred)) {
    maxpred <- max(y)
  }
  if (missing(ylimit)) {
    ylimit <- max(tabulate(y)) * 1.1
  }
  
  df.bar <- barplot(
    tabulate(y)[1:maxpred],
    names = 1:maxpred,
    xlab = "count",
    ylab = "frequency",
    col = "gray",
    ylim = c(0, ylimit)
  )
  
  leg <- c("actual data")
  cols <- c("gray")
  pchs <- c(15)
  
  if (inherits(model, "oneinflmodel") && model$dist == "Poisson") {
    out1 <- plotoipp(model, data, maxpred)
    modelcounter <- model
    class(modelcounter) <- "truncmodel"
    out2 <- plotpp(modelcounter, data, maxpred)
    leg <- c(leg, out1$label, out2$label)
    cols <- c(cols, out1$color, out2$color)
    pchs <- c(pchs, out1$pch, out2$pch)
  } else if (inherits(model, "oneinflmodel") && model$dist == "negbin") {
    out1 <- plotoiztnb(model, data, maxpred)
    modelcounter <- model
    class(modelcounter) <- "truncmodel"
    out2 <- plotztnb(modelcounter, data, maxpred)
    leg <- c(leg, out1$label, out2$label)
    cols <- c(cols, out1$color, out2$color)
    pchs <- c(pchs, out1$pch, out2$pch)
  }
  
  legend("topright", legend = leg, col = cols, pch = pchs, cex = 1)
}