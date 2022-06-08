test_that("Fitting a capped logistic spline for detection proportion works", {
  data("datasetA")
  eg.nuis <- getNuisance(log2(datasetA$prot))
  eg.start <- logitSplines_start(dp = eg.nuis$dp,
                                 mu = eg.nuis$mu_obs,
                                 wt = eg.nuis$wt,
                                 df = 5)
  eg.fit <- cappedLogit_ztbinom(dp = eg.nuis$dp,
                                wt = eg.nuis$wt,
                                X = eg.start$X,
                                beta0 = eg.start$betas_start,
                                alpha0 = 0.9,
                                trace = FALSE)
  expect_type(eg.fit, "list")
  expect_length(eg.fit, 2)
  expect_equal(ncol(eg.fit[[2]]), 8)
  expect_equal(colnames(eg.fit[[2]])[1], "alpha")
  expect_identical(length(eg.fit[[1]]), length(eg.start[[1]])+1L)
  expect_identical(length(eg.fit[[1]]), ncol(eg.start[[2]])+2L)
})
