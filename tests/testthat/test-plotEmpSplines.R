test_that("Plotting an empirical spline works", {
  data("datasetA")
  eg.nuis <- getNuisance(log2(datasetA$prot))
  # df=1
  eg.start1 <- logitSplines_start(dp = eg.nuis$dp,
                                  mu = eg.nuis$mu_obs,
                                  wt = eg.nuis$wt,
                                  df = 1)
  eg.fit1 <- logit_ztbinom(dp = eg.nuis$dp,
                           wt = eg.nuis$wt,
                           X = eg.start1$X,
                           beta0 = eg.start1$betas_start)
  # df=5
  eg.start5 <- logitSplines_start(dp = eg.nuis$dp,
                                  mu = eg.nuis$mu_obs,
                                  wt = eg.nuis$wt,
                                  df = 5)
  eg.fit5 <- logit_ztbinom(dp = eg.nuis$dp,
                           wt = eg.nuis$wt,
                           X = eg.start5$X,
                           beta0 = eg.start5$betas_start)
  theplot <- function() {
    # raw scale
    set.seed(104)
    par(mfrow = c(1, 2))
    plotEmpSplines(eg.nuis,
                   X = eg.start1$X,
                   params = eg.fit1$params,
                   jitter.amount = 0.1,
                   ylim = c(0, 1.04))
    plotEmpSplines(eg.nuis,
                   X = eg.start5$X,
                   params = eg.fit5$params,
                   ylim = c(0, 1.04),
                   lineCol = "red",
                   newPlot = FALSE)
    # logit scale for df = 1
    plotEmpSplines(eg.nuis,
                   X = eg.start1$X,
                   params = eg.start1$betas_start,
                   ylim = c(-2, 12),
                   logit = TRUE,
                   lty = "dotted")
    plotEmpSplines(eg.nuis,
                   X = eg.start1$X,
                   params = eg.fit1$params,
                   ylim = c(-2, 12),
                   logit = TRUE,
                   newPlot = FALSE)
  }

  vdiffr::expect_doppelganger("Empirical splines", theplot)
})
