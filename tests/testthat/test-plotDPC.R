test_that("Plotting a detection probability curve works", {
  data("datasetA")
  eg.nuis <- getNuisance(log2(datasetA$prot))
  eg.dpc <- dpc(eg.nuis)
  theplot <- function() {
    set.seed(104)
    plotDPC(eg.nuis, eg.dpc, ylim = c(0, 1.04))
  }

  vdiffr::expect_doppelganger("DPC plot", theplot)
})
