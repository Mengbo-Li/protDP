test_that("Fitting a logistic spline for detection proportion works", {
  data("datasetA")
  eg.nuis <- getNuisance(log2(datasetA$prot))
  eg.start <- logitSplines_start(dp = eg.nuis$dp,
                                 mu = eg.nuis$mu_obs,
                                 wt = eg.nuis$wt,
                                 df = 5)
  eg.fit <- logit_ztbinom(dp = eg.nuis$dp,
                          wt = eg.nuis$wt,
                          X = eg.start$X,
                          beta0 = eg.start$betas_start)
  expect_type(eg.fit, "list")
  expect_length(eg.fit, 2)
  expect_equal(dim(eg.fit[[2]]), c(2, 7))
  expect_identical(length(eg.fit[[1]]), length(eg.start[[1]]))
  expect_identical(length(eg.fit[[1]]), ncol(eg.start[[2]])+1L)
})
