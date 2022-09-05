test_that("Plotting a detection probability curve works", {
  testthat::skip_on_os(c("linux", "windows"))
  data("datasetA")
  eg.dpc <- dpc(log2(datasetA$prot))
  theplot <- function() {
    set.seed(104)
    plotDPC(eg.dpc, ylim = c(0, 1.04), add.jitter = FALSE)
  }

  vdiffr::expect_doppelganger("DPC plot", theplot)
})
