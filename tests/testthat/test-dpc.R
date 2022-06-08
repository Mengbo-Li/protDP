test_that("Fitting a detection probability curve works", {
  data("datasetA")
  eg.nuis <- getNuisance(log2(datasetA$prot))
  eg.dpc <- dpc(eg.nuis)

  expect_type(eg.dpc, "list")
  expect_length(eg.dpc, 3)
  expect_length(eg.dpc[[1]], 2)
  expect_equal(ncol(eg.dpc[[2]]), 3)
  expect_type(eg.dpc[[3]], "double")
  expect_equal(length(eg.dpc[[3]]), nrow(datasetA$prot))
  expect_false(any(eg.dpc[[1]] == eg.nuis[[5]]))

})
