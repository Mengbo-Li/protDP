#' Wrapper function of all presented results
#'
#' @param data The log2-transformed intensity matrix where rows are precursors
#' and columns are samples.
#'
#' @return List of results including the logistic-spline fits, the deviance
#' partitioning plot, the capped logit-linear fit and the fitted detection
#' probability curve assuming normal observed intensities.
#'
#' @examples
#' ## See the vignettes.
#'
#' @export
gatherResults <- function(data) {
  nuis <- getNuisance(data)
  # first get all logistic-splines fits
  splineFits <- list()
  splineFits_params0 <- list()
  for (df in c(1, 3, 5)) {
    params0 <- logitSplines_start(dp = nuis$dp, mu = nuis$mu_obs,
                                  wt = nuis$wt, df = df)
    fit <- logit_ztbinom(dp = nuis$dp, X = params0$X,
                         wt = nuis$wt, beta0 = params0$betas_start)
    splineFits_params0[[(df+1)/2]] <- params0
    splineFits[[(df+1)/2]] <- fit
  }
  # next the deviance reduced
  baseline_params0 <- logitSplines_start(dp = nuis$dp, mu = nuis$mu_obs,
                                         wt = nuis$wt, df = 0)
  baselineFit <- logit_ztbinom(dp = nuis$dp, X = baseline_params0$X,
                               wt = nuis$wt, beta0 = baseline_params0$betas_start)
  baselineDev <- 2*min(baselineFit$info[, "neg.LL"])
  devs <- data.frame(df = c(0, 1, 3, 5),
                     dev = c(baselineDev,
                             sapply(splineFits, function(i)
                               2*min(i$info[, "neg.LL"]))))
  devs$dev.decre <- baselineDev - devs$dev
  devs$percDevReduced = devs$dev.decre / max(devs$dev.decre)
  # next the capped logit-linear fit
  cappedLinear_params0 <- splineFits_params0[[1]]
  cappedLinearFit <- cappedLogit_ztbinom(dp = nuis$dp, X = cappedLinear_params0$X,
                                         wt = nuis$wt, alpha0 = 0.9,
                                         beta0 = cappedLinear_params0$betas_start,
                                         trace = FALSE)
  # Detection probability curve assuming normal observed intensities
  dpcFit <- dpc(nuis)
  return(list(nuis = nuis,
              splineFits_params0 = splineFits_params0,
              splineFits = splineFits,
              devs = devs,
              cappedLinearFit = cappedLinearFit,
              dpcFit = dpcFit))
}
