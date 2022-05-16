#' Detection probability curve for mass spectrometry-based proteomics data
#' assuming observed normal intensities
#'
#' @param Y The log2-transformed intensity matrix where rows are precursors and
#' columns are samples.
#' @param nuis The list of nuisance parameters.
#' @param maxit Maximum number of iterations.
#' @param eps Convergence tolerance.
#' @param b1.upper Upper bound for beta_1. Typically 1.
#' @param capped Logical. Whether to fit a capped logit-linear model.
#' @param alpha0 A start value for alpha if capped is TRUE.
#' @param trace Logical. Whether to print out alpha estimates at each iteration.
#'
#' @return Fitted parameters and the fitting history.
#' @export
#'
#' @examples
#' ## See the vignettes.
dpc <- function(Y, nuis, maxit = 100, eps = 1e-4, b1.upper = 1,
                    capped = FALSE, alpha0 = NA, trace = TRUE) {
  if (!capped) {
    dp <- nuis$dp
    wt <- nuis$wt
    s2 <- nuis$s2
    mu_obs <- nuis$mu_obs
    mu_mis <- mu_obs - nuis$betaStart[2]*s2
    betas <- nuis$betaStart
    betas.hist <- betas
    negLL <- dpc_ztbinom.negLL(betas, dp, wt, mu_obs, mu_mis, alpha = 1)
    negLL.hist <- negLL
    for (i in 1:maxit) {
      ztbinomFit <- optim(betas, dpc_ztbinom.negLL,
                          dp = dp, wt = wt, mu_obs = mu_obs, mu_mis = mu_mis, alpha = 1,
                          method = "L-BFGS-B", lower = c(-Inf, 0), upper = c(0, b1.upper))
      newBetas <- ztbinomFit$par
      mu_mis <- mu_obs - newBetas[2]*s2
      newNegLL <- dpc_ztbinom.negLL(newBetas, dp, wt, mu_obs, mu_mis, alpha = 1)
      if (negLL - newNegLL < eps) break
      betas.hist <- rbind(betas.hist, newBetas)
      negLL.hist <- c(negLL.hist, newNegLL)
      betas <- newBetas
      negLL <- newNegLL
    }
    info <- cbind(betas.hist, negLL.hist)
    colnames(info) <- c("b0", "b1", "neg.ZBLL")
    rownames(info) <- paste("i=", 0:(nrow(info)-1), sep = "")
    return(list(beta = betas, hist = info, mu_mis = mu_mis))
  }
  if (capped) {
    dp <- nuis$dp
    wt <- nuis$wt
    s2 <- nuis$s2
    mu_obs <- nuis$mu_obs
    mu_mis <- mu_obs - nuis$betaStart[2]*s2
    if (is.na(alpha0)) stop("Please supply a starting value for alpha!")
    alpha <- alpha0
    alpha.hist <- alpha
    # get the MLE for betas given alpha0
    ztbinomFit <- optim(nuis$betaStart, dpc_ztbinom.negLL,
                        dp = dp, wt = wt, mu_obs = mu_obs, mu_mis = mu_mis, alpha = alpha,
                        method = "L-BFGS-B", lower = c(-Inf, 0), upper = c(0, b1.upper))
    betas <- ztbinomFit$par
    betas.hist <- betas
    mu_mis <- mu_obs - betas[2]*s2
    # record the start value for negative LL
    negLL <- dpc_ztbinom.negLL(betas, dp, wt, mu_obs, mu_mis, alpha = alpha)
    negLL.hist <- negLL
    for (i in 1:maxit){
      # given MLE of the betas, re-estimate alpha
      eta <- betas[1] + 0.5*betas[2]*(mu_obs + mu_mis)
      p <- plogis(eta)
      mu <- alpha * p
      newAlpha <- sum(wt*dp/(1-mu)) / sum(wt*p/(1-mu) + p*(1-mu)^(wt-1)/iBeta(mu, 1, wt))
      if (newAlpha >= 1) {
        alpha <- 1
        break
      }
      if (abs(newAlpha - alpha) < eps) break
      if (trace) cat(paste("alpha_", i, sep = ""), "=", newAlpha, "\n")
      # re-estimate MLE for betas with the new alpha
      ztbinomFit <- optim(betas, dpc_ztbinom.negLL,
                          dp = dp, wt = wt, mu_obs = mu_obs, mu_mis = mu_mis, alpha = newAlpha,
                          method = "L-BFGS-B", lower = c(-Inf, 0), upper = c(0, b1.upper))
      newBetas <- ztbinomFit$par
      mu_mis <- mu_obs - newBetas[2]*s2
      newNegLL <- dpc_ztbinom.negLL(newBetas, dp, wt, mu_obs, mu_mis, alpha = newAlpha)
      alpha.hist <- c(alpha.hist, newAlpha)
      betas.hist <- rbind(betas.hist, newBetas)
      negLL.hist <- c(negLL.hist, newNegLL)
      alpha <- newAlpha
      betas <- newBetas
    }
    info <- cbind(alpha.hist, betas.hist, negLL.hist)
    colnames(info) <- c("alpha", "b0", "b1", "neg.ZBLL")
    rownames(info) <- paste("i=", 0:(nrow(info)-1), sep = "")
    return(list(alpha = alpha, beta = betas, hist = info, mu_mis = mu_mis))
  }
}
