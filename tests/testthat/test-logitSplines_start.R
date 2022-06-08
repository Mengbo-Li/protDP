test_that("Getting start values for parameters of regression splines work", {
  data("datasetA")
  eg.nuis <- getNuisance(log2(datasetA$prot))
  eg <- logitSplines_start(dp = eg.nuis$dp,
                           mu = eg.nuis$mu_obs,
                           wt = eg.nuis$wt,
                           df = 5)
  expect_type(eg, "list")
  expect_length(eg, 2)
  expect_length(eg[[1]], 6)
  expect_equal(ncol(eg[[2]]), 5)
})
