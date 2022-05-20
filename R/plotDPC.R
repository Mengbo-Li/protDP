#' Plot the detection probability curve
#'
#' @param nuis The list of nuisance parameters.
#' @param dpcFit The DPC fit results.
#' @param add.jitter Logical. Whether to add jitter to the detected proportion
#' axis.
#' @param jitter.amount Amount of jittering.
#' @param point.cex Size of the points.
#' @param lwd Line width.
#' @param ylim Limits of the y-axis.
#' @param capped Logical. Whether the probabilities are capped by alpha.
#'
#' @return The plot.
#' @export
#'
#' @examples
#' # See the vignettes.
plotDPC <- function(nuis, dpcFit,
                    add.jitter = TRUE, jitter.amount = NULL,
                    point.cex = 0.2, lwd = 2, ylim = c(0, 1),
                    capped = FALSE) {
  x <- (nuis$mu_obs + dpcFit$mu_mis)/2
  y <- nuis$dp
  if (add.jitter) y <- jitter(y, amount = jitter.amount)
  plot(x = x, y = y,
       pch = 16, cex = point.cex, ylim = ylim,
       xlab = expression(frac(mu[obs] + mu[mis], 2)), ylab = "Detection probability",
       main = "DPC: Average log intensity vs detection probability")
  if (!capped) {
    x <- x[order(x)]
    lines(x = x, y = plogis(nuis$betaStart[1] + nuis$betaStart[2]*x),
          lty = "dashed", lwd = lwd)
    lines(x = x, y = plogis(dpcFit$beta[1] + dpcFit$beta[2]*x),
          col = "blue", lwd = lwd)
    legend("bottomright", legend = c("Start", "Final"), col = c("black", "blue"), lty = c(2, 1), lwd = lwd)
  } else {
    x <- x[order(x)]
    lines(x = x, y = dpcFit$hist["i=0", "alpha"] * plogis(dpcFit$hist["i=0", "b0"] + dpcFit$hist["i=0", "b1"]*x),
          lty = "dashed", lwd = lwd)
    lines(x = x, y = dpcFit$alpha * plogis(dpcFit$beta[1] + dpcFit$beta[2]*x),
          col = "blue", lwd = lwd)
    legend("bottomright", legend = c("Start", "Final"), col = c("black", "blue"), lty = c(2, 1), lwd = lwd)
  }
}
