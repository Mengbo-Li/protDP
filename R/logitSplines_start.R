#' Get the initial estimates for a logistic spline fit by assuming the
#' binomial detection probability
#'
#' @param dp A vector of detection proportion in all precursors.
#' @param mu A vector of average observed intensities in all precursors.
#' @param wt A vector of the numbers of trials (samples size).
#' @param df Degrees of freedom for the natural cubic spline.
#'
#' @return A list of fitted coefficients and the basis matrix for a spline with
#' some degrees of freedom.
#'
#' @examples
#' # See the vignettes.
#' @importFrom splines ns
#' @export
logitSplines_start <- function(dp, mu, wt, df = 1) {
  if (df == 0) {
    fit <- glm(dp ~ 1, weights = wt, family = binomial)
    X <- matrix(1, nrow = length(dp), ncol = 1)
  }
  else {
    X <- splines::ns(mu, df = df)
    fit <- glm(dp ~ X, weights = wt, family = binomial)
  }
  names(fit$coefficients) <- paste("b", 0:df, sep = "")
  list(betas_start = fit$coefficients, X = X)
}
