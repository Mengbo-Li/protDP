test_that("Fitting a detection probability curve works", {
  data("datasetA")
  eg.dpc <- dpc(log2(datasetA$prot))

  expect_type(eg.dpc, "list")
  expect_length(eg.dpc, 10)
  expect_length(eg.dpc[[1]], 2)
  expect_equal(ncol(eg.dpc[[2]]), 3)
  expect_type(eg.dpc[[3]], "double")
  expect_equal(length(eg.dpc[[4]]), nrow(datasetA$prot))

})
