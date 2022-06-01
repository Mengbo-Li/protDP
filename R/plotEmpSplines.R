#' Plot the fitted empirical spline of detected proportions to average observed
#' intensities
#'
#' @param nuis The list of nuisance parameters.
#' @param X The basis matrix for the natural cubic spline.
#' @param params Fitted coefficients.
#' @param capped Logical. Whether the probabilities are capped by alpha.
#' @param logit Logical. Whether to plot at the logit scale.
#' @param add.jitter Logical. Whether to add jitter to the detected proportion
#' axis.
#' @param jitter.amount Amount of jittering.
#' @param point.cex Size of the points.
#' @param plot.dotts Logical. Whether to keep the raw data points in the plot.
#' @param lineOnly Logical. Whether to plot only fitted lines.
#' @param lty Line type.
#' @param lineCol Line color.
#' @param lwd Line width.
#' @param ylim Limits of the y-axis.
#'
#' @export
#'
#' @examples
#' # See the vignettes.
#' @importFrom graphics legend
#' @importFrom graphics lines
plotEmpSplines <- function(nuis, X, params,
                           capped = TRUE, logit = FALSE,
                           add.jitter = TRUE, jitter.amount = NULL,
                           point.cex = 0.2, plot.dotts = TRUE,
                           lineOnly = FALSE, lty = "solid",
                           lineCol = "darkgreen", lwd = 2,
                           ylim = c(0, 1)) {
  x <- nuis$mu_obs
  y <- nuis$dp
  if (!logit) {
    if (add.jitter) y <- jitter(y, amount = jitter.amount)
    if (plot.dotts)
      plot(x = x, y = y,
           pch = 16, cex = point.cex, ylim = ylim,
           xlab = "Average observed intensity", ylab = "Detection proportion",
           main = "Empirical splines")
    X <- cbind(1, X)
    if (capped) {
      alpha <- params[1]
      betas <- params[-1]
    } else {
      alpha <- 1
      betas <- params
    }
    df <- length(betas) - 1
    eta <- colSums(t(X) * betas)
    fitte.dp <- alpha * plogis(eta)
    lines(x[order(x)], fitte.dp[order(x)], col = lineCol, lwd = lwd, lty = lty)
  } else {
    y.logit <- qlogis(y)
    if (add.jitter) y.logit <- jitter(y.logit, amount = jitter.amount)
    if (plot.dotts)
      plot(x = x, y = y.logit,
           pch = 16, cex = point.cex, ylim = ylim,
           xlab = "Average observed intensity", ylab = "logit(detected proportion)",
           main = "Empirical splines: Logit scale")
    X <- cbind(1, X)
    if (capped) {
      alpha <- params[1]
      betas <- params[-1]
    } else {
      alpha <- 1
      betas <- params
    }
    df <- length(betas) - 1
    eta <- colSums(t(X) * betas)
    fitte.dp <- alpha * plogis(eta)
    fitte.eta <- qlogis(fitte.dp)
    if (lineOnly) plot(x[order(x)], fitte.eta[order(x)], col = lineCol, lwd = lwd, type = "l",
                       ylim = ylim, lty = lty,
                       xlab = "Average observed intensity", ylab = "logit(detected proportion)",
                       main = "Empirical splines: Logit scale")
    else lines(x[order(x)], fitte.eta[order(x)], col = lineCol, lwd = lwd, lty = lty)
  }
}
