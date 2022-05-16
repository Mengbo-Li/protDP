#' Negative log-likelihood under the zero-truncated binomial distribution
#' This is the objective function to minimise in `logit_ztbinom`.
#' @noRd
logit_ztbinom.ZT_negLL <- function(params, dp, wt, X) {
  df <- length(params) - 1
  if (df > 0) X <- cbind(1, X)
  eta <- colSums(t(X) * params)
  p <- plogis(eta)
  -sum(dbinom(dp*wt, size = wt, p = p, log = TRUE) -
         pbinom(0.5, size = wt, p = p, log.p = TRUE, lower.tail = FALSE))
}
