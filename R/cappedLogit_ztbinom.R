#' Fitting an empirical logistic spline curve for detection proportion with
#' capped probabilities
#'
#' @param dp A vector of detection proportion in all precursors.
#' @param X The basis matrix for a natural cubic spline.
#' @param wt A vector of the numbers of trials (samples size) for zero-truncated
#' binomial distribution.
#' @param alpha0 A start value for alpha.
#' @param beta0 Start values for the beta coefficients.
#' @param maxit Maximum number of iterations.
#' @param trace Logical. Whether to print out alpha estimates at each iteration.
#' @param eps Convergence tolerance.
#'
#' @return Fitted parameters and the fitting history.
#' @export
#'
#' @examples
#' ## See the vignettes.
cappedLogit_ztbinom <- function(dp, X, wt, alpha0 = 0.8, beta0, maxit = 100,
                                trace = TRUE, eps = 1e-4) {

  df <- length(beta0) - 1
  alpha <- alpha0
  alpha.hist <- alpha
  # the order of parameters goes b0, b1, and etc.
  lower.bounds <- c(-Inf, rep(0, df))
  if (df > 0) upper.bounds <- c(0, rep(Inf, df))
  else upper.bounds <- c(Inf)
  # MLE for betas at the start value of alpha0
  ztbinomFit <- optim(beta0, cappedZT_negLL_givenAlpha,
                      dp = dp, wt = wt, X = X, alpha = alpha,
                      method = "L-BFGS-B",
                      lower = lower.bounds, upper = upper.bounds)
  betas <- ztbinomFit$par
  betas.hist <- matrix(betas, nrow = 1)
  negLL <- cappedZT_negLL_givenAlpha(betas, dp, wt, X, alpha)
  negLL.hist <- negLL

  for (i in 1:maxit) {
    if (df > 0)
      eta <- colSums(t(cbind(1, X)) * betas)
    if (df == 0)
      eta <- colSums(t(X) * betas)
    p <- plogis(eta)
    mu <- alpha * p
    newAlpha <- sum(wt * dp/(1 - mu))/sum(wt * p/(1 - mu) +
                                            p * (1 - mu)^(wt - 1)/iBeta(mu, 1, wt))
    if (newAlpha >= 1) {
      alpha <- 1
      if (trace)
        cat(paste("alpha_", i, sep = ""), "=", alpha, "\n")
      ztbinomFit <- optim(betas, cappedZT_negLL_givenAlpha,
                          dp = dp, wt = wt, X = X, alpha = alpha, method = "L-BFGS-B",
                          lower = lower.bounds, upper = upper.bounds)
      newBetas <- ztbinomFit$par
      newNegLL <- cappedZT_negLL_givenAlpha(newBetas, dp, wt, X, alpha)
      alpha.hist <- c(alpha.hist, alpha)
      betas.hist <- rbind(betas.hist, newBetas)
      negLL.hist <- c(negLL.hist, newNegLL)
      betas <- newBetas
      break
    }
    if (abs(newAlpha - alpha) < eps)
      break
    if (trace)
      cat(paste("alpha_", i, sep = ""), "=", newAlpha,
          "\n")
    ztbinomFit <- optim(betas, cappedZT_negLL_givenAlpha,
                        dp = dp, wt = wt, X = X, alpha = newAlpha, method = "L-BFGS-B",
                        lower = lower.bounds, upper = upper.bounds)
    newBetas <- ztbinomFit$par
    newNegLL <- cappedZT_negLL_givenAlpha(newBetas, dp,
                                          wt, X, newAlpha)
    alpha.hist <- c(alpha.hist, newAlpha)
    betas.hist <- rbind(betas.hist, newBetas)
    negLL.hist <- c(negLL.hist, newNegLL)
    alpha <- newAlpha
    betas <- newBetas
  }

  # clean up results
  info <- cbind(alpha.hist, betas.hist, negLL.hist)
  colnames(info) <- c("alpha", names(beta0), "neg.ZBLL")
  rownames(info) <- paste("i=", 0:(nrow(info)-1), sep = "")
  final.params <- c(alpha, betas)
  names(final.params)[1] <- "alpha"

  list(params = final.params, info = info)
}
