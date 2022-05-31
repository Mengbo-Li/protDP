#' Get nuisance parameters and an initial estimation of the detection
#' probability curve (DPC)
#'
#' The `getNuisance` function fits a logistic-linear model for detection
#' proportion to the observed average intensity in each precursor. It also
#' returns a list of nuisance parameters for convenience of use when estimating
#' DPC, such as the vector of observed mean intensity in each precursor, the
#' total samples size and the eBayes-trended observed variances in each
#' precursor.
#'
#' @param Y The log2-transformed intensity matrix where rows are precursors and
#' columns are samples.
#'
#' @return A list of nuisance parameters for use to estimate the detection
#' probability curve.
#'
#' @examples
#' ## NOT RUN
#'
#' @export

getNuisance <- function(Y) {
  wt <- rep(ncol(Y), nrow(Y))
  dp <- rowMeans(!is.na(Y))
  mu_obs <- rowMeans(Y, na.rm = TRUE)
  require(limma)
  fit0 <- eBayes(lmFit(Y, design = matrix(1, nrow = ncol(Y))), trend = TRUE)
  s2 <- fit0$s2.prior
  glmfit0 <- glm(dp ~ mu_obs, weight = wt, family = binomial)
  betaStart <- as.numeric(coef(glmfit0))
  return(list(wt = wt,
              dp = dp,
              mu_obs = mu_obs,
              s2 = s2,
              betaStart = betaStart,
              fit0 = fit0,
              glmfit0 = glmfit0))
}
