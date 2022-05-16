#' Log-likelihood under zero-truncated binomial
#' This is the objective function to minimise in `cappedLogit_ztbinom`
#' @noRd
cappedZT_negLL_givenAlpha <- function(params, dp, wt, X, alpha) {
   # df - 1 because there is an intercept
   df <- length(params) - 1
   betas <- params
   if (df > 0) X <- cbind(1, X)
   eta <- colSums(t(X) * betas)
   p <- alpha * plogis(eta)
   # negative log-liklihood under zt-binom, to minimise in optim
   -sum(dbinom(dp*wt, size = wt, p = p, log = TRUE) -
           pbinom(0.5, size = wt, p = p, log.p = TRUE, lower.tail = FALSE))
}
