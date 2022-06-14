#' Fitting an empirical logistic spline curve for detection proportion
#'
#' @param dp A vector of detection proportion in all precursors.
#' @param X The basis matrix for a natural cubic spline.
#' @param wt A vector of the numbers of trials (samples size) for zero-truncated
#' binomial distribution.
#' @param beta0 Start values for the beta coefficients.
#' @param b0.upper Upper bound for b0.
#' @param b1.upper Upper bound for b1.
#'
#' @return Fitted beta coefficients and the fitting history.
#' @examples
#' ## See the vignettes.
#' @importFrom stats optim pbinom dbinom binomial coef glm pbeta plogis qlogis
#' @export
logit_ztbinom <- function(dp, X, wt, beta0, b0.upper = 0, b1.upper = Inf) {

  df <- length(beta0) - 1
  params <- beta0
  params.hist <- matrix(params, nrow = 1)
  negLL <- logit_ztbinom.ZT_negLL(params, dp, wt, X)
  negLL.hist <- negLL


  # the order of parameters goes alpha, b0, b1, and etc.
  lower.bounds <- c(-Inf, rep(0, df))
  if (df > 0) upper.bounds <- c(b0.upper, b1.upper, rep(Inf, df-1))
  else upper.bounds <- c(Inf)
  ztbinomFit <- stats::optim(params,
                             logit_ztbinom.ZT_negLL,
                             dp = dp, wt = wt, X = X,
                             method = "L-BFGS-B",
                             lower = lower.bounds,
                             upper = upper.bounds)
  newParams <- ztbinomFit$par
  newnegLL <- logit_ztbinom.ZT_negLL(newParams, dp, wt, X)
  params.hist <- rbind(params.hist, newParams)
  negLL.hist <- c(negLL.hist, newnegLL)


  # clean up results
  info <- cbind(params.hist, negLL.hist)
  colnames(info) <- c(names(beta0), "neg.LL")
  rownames(info) <- paste("iter", 0:(nrow(info)-1), sep = " ")

  list(params = newParams, info = info)

}
