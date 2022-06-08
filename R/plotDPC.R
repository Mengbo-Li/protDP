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
#'
#' @return The plot.
#' @export
#'
#' @examples
#' # See the vignettes.
plotDPC <- function(nuis,
                    dpcFit,
                    add.jitter = TRUE,
                    jitter.amount = NULL,
                    point.cex = 0.2,
                    lwd = 2,
                    ylim = c(0, 1)) {
  x <- (nuis$mu_obs + dpcFit$mu_mis)/2
  y <- nuis$dp
  if (add.jitter) y <- jitter(y, amount = jitter.amount)
  plot(x = x, y = y,
       pch = 16, cex = point.cex, ylim = ylim,
       xlab = "Intensity", ylab = "Detection probability",
       main = "Detection probability curve")
  x <- x[order(x)]
  lines(x = x, y = plogis(nuis$betaStart[1] + nuis$betaStart[2]*x),
        lty = "dashed", lwd = lwd)
  lines(x = x, y = plogis(dpcFit$beta[1] + dpcFit$beta[2]*x),
        col = "blue", lwd = lwd)
  legend("bottomright", legend = c("Start", "Final"),
         col = c("black", "blue"), lty = c(2, 1), lwd = lwd, cex = 0.8)
}
