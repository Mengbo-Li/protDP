#' Plot the fitted empirical spline of detected proportions to average observed
#' intensities
#'
#' @param nuis The list of nuisance parameters.
#' @param X The basis matrix for the natural cubic spline.
#' @param params Fitted coefficients.
#' @param capped Logical. Whether the probabilities are capped by alpha.
#' @param add.jitter Logical. Whether to add jitter to the detected proportion
#' axis.
#' @param jitter.amount Amount of jittering.
#' @param point.cex Size of the points.
#' @param newPlot Logical. Whether to start a new plot device.
#' @param logit Logical. Whether to plot at the logit scale.
#' @param lty Line type.
#' @param lineCol Line color.
#' @param lwd Line width.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#'
#' @export
#'
#' @examples
#' # See the vignettes.
#' @importFrom graphics legend
#' @importFrom graphics lines
plotEmpSplines <- function(nuis,
                           X,
                           params,
                           capped = FALSE,
                           add.jitter = TRUE,
                           jitter.amount = NULL,
                           point.cex = 0.2,
                           newPlot = TRUE,
                           logit = FALSE,
                           lty = "solid",
                           lineCol = "darkgreen",
                           lwd = 2,
                           xlim = NULL,
                           ylim = c(0, 1)) {

  x <- nuis$mu_obs
  y <- nuis$dp
  if (!logit) {
    if (add.jitter) y <- jitter(y, amount = jitter.amount)
    if (newPlot)
      plot(x = x,
           y = y,
           pch = 16,
           cex = point.cex,
           xlim = xlim,
           ylim = ylim,
           xlab = "Average observed intensity",
           ylab = "Detection proportion",
           main = "Empirical splines")
    X <- cbind(1, X)
    if (capped) {
      alpha <- params[1]
      betas <- params[-1]
    } else {
      alpha <- 1
      betas <- params
    }
    eta <- colSums(t(X) * betas)
    fitte.dp <- alpha * plogis(eta)
    lines(x[order(x)], fitte.dp[order(x)], col = lineCol, lwd = lwd, lty = lty)
  } else {
    # when plotting at the logit scale, there is no capped functionality
    X <- cbind(1, X)
    betas <- params
    eta <- colSums(t(X) * betas)
    if (newPlot) plot(x[order(x)],
                      eta[order(x)],
                      col = lineCol,
                      lwd = lwd,
                      type = "l",
                      xlim = xlim,
                      ylim = ylim,
                      lty = lty,
                      xlab = "Average observed intensity",
                      ylab = "logit(detected proportion)",
                      main = "Empirical splines: Logit scale")
    else lines(x[order(x)], eta[order(x)], col = lineCol, lwd = lwd, lty = lty)
  }
}
