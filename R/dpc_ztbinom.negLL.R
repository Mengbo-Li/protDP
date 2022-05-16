#' Negative log-likelihood for zero-truncated binomial distribution
#' The objective function for `optim` to minimise
#' alpha = 1 for non-capped DPC
#' alpha not default when capped = TRUE in `dpc`
#' @noRd
dpc_ztbinom.negLL <- function(params, dp, wt, mu_obs, mu_mis, alpha = 1) {
  b0 <- params[1]
  b1 <- params[2]
  eta <- b0 + 0.5*b1*(mu_obs + mu_mis)
  p <- alpha * plogis(eta)
  -sum(dbinom(dp*wt, size = wt, p = p, log = TRUE) -
         pbinom(0.5, size = wt, p = p, log.p = TRUE, lower.tail = FALSE))
}
