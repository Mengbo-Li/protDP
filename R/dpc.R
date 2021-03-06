#' Detection probability curve for label free shotgun proteomics data
#' assuming observed normal intensities
#'
#' @param nuis The list of nuisance parameters.
#' @param maxit Maximum number of iterations.
#' @param eps Convergence tolerance.
#' @param b1.upper Upper bound for beta_1. Typically 1.
#'
#' @return Fitted parameters and the fitting history.
#' @export
#'
#' @examples
#' ## See the vignettes.
dpc <- function(nuis, maxit = 100, eps = 1e-4, b1.upper = 1) {
  dp <- nuis$dp
  wt <- nuis$wt
  s2 <- nuis$s2
  mu_obs <- nuis$mu_obs
  fit0 <- logit_ztbinom(dp = nuis$dp,
                        X = matrix(nuis$mu_obs, ncol = 1),
                        wt = nuis$wt,
                        beta0 = nuis$betaStart,
                        b0.upper = 0,
                        b1.upper = b1.upper)
  betaStart <- fit0$params
  mu_mis <- mu_obs - betaStart[2]*s2
  betas <- betaStart
  betas.hist <- matrix(betas, nrow = 1)
  negLL <- dpc_ztbinom.negLL(betas, dp, wt, mu_obs, mu_mis)
  negLL.hist <- negLL
  for (i in 1:maxit) {
    ztbinomFit <- stats::optim(betas, dpc_ztbinom.negLL,
                               dp = dp, wt = wt,
                               mu_obs = mu_obs, mu_mis = mu_mis,
                               method = "L-BFGS-B", lower = c(-Inf, 0),
                               upper = c(0, b1.upper))
    newBetas <- ztbinomFit$par
    mu_mis <- mu_obs - newBetas[2]*s2
    newNegLL <- dpc_ztbinom.negLL(newBetas, dp, wt, mu_obs, mu_mis)
    if (negLL - newNegLL < eps) break
    betas.hist <- rbind(betas.hist, newBetas)
    negLL.hist <- c(negLL.hist, newNegLL)
    betas <- newBetas
    negLL <- newNegLL
  }
  info <- cbind(betas.hist, negLL.hist)
  colnames(info) <- c("b0", "b1", "neg.ZBLL")
  rownames(info) <- paste("i=", 0:(nrow(info)-1), sep = "")
  return(list(beta = betas, hist = info, betaStart = betaStart, mu_mis = mu_mis))
}
