#' Estimation of hyperparameters for the empirical Bayes method
#'
#' @param y Log2 transformed precursor-level intensity data.
#'
#' @return List of hyperparameters and the posterior means and variances.
#' @importFrom stats var
#' @export
#'
#' @examples
#' # See vignettes.
hyperparams <- function(y) {
  mu <- rowMeans(y, na.rm = TRUE)
  s2 <- apply(y, 1, var, na.rm = TRUE)
  df.residual <- rowSums(!is.na(y)) - 1
  # Estimate variance hyperparameters
  sv <- limma::squeezeVar(s2, df = df.residual)
  df.prior <- sv$df.prior
  df.total <- df.prior + df.residual
  s2.prior <- sv$var.prior
  s2.post <- sv$var.post
  # Estimate mu0
  mu0 <- mean(mu)
  if (is.infinite(df.prior)) {
    n0 <- Inf
    mu_obs.post <- mu
    s2_obs.post <- s2
  } else {
    # Estimate n0
    n <- df.residual + 1
    modt <- (mu - mu0) / sqrt(s2.post/n)
    n0 <- n0EstimateFromModT(tstat = modt,
                             df = df.total,
                             n = n,
                             niter = 20L,
                             eps = 1e-5,
                             trace = FALSE)
    # posterior mean and variance
    mu_obs.post <- (n*mu + n0*mu0) / (n + n0)
    s2_obs.post <- (n*n0*(mu-mu0)^2/(n+n0) + (n-1)*s2 +
                      df.prior*s2.prior) /(n+df.prior)
    s2_obs.post[is.na(s2_obs.post)] <- s2.prior
  }
  list(mu0 = mu0,
       n0 = n0,
       df.prior = df.prior,
       s2.prior = s2.prior,
       mu_obs.post = mu_obs.post,
       s2_obs.post = s2_obs.post)
}


