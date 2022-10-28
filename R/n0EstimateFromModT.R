#' Estimate the effective sample size hyperparamter n0
#' This function is called in the `hyperparams` function.
#' Created 4 Aug 2022. Last modified 26 Oct 2022.
#' @noRd
n0EstimateFromModT <- function(tstat,
                               df,
                               n,
                               niter = 20L,
                               eps = 1e-5,
                               trace = FALSE) {
  tstat <- as.numeric(tstat)
  ExpectedLogF <- mean(logmdigamma(df/2) - logmdigamma(1/2))
  RHS <- 2*mean(log(abs(tstat))) - ExpectedLogF

  # Monotonically convergent Newton iteration for v0 = 1/n0
  # Solve mean(log(1+n*v0)) = RHS
  if(RHS <= 0) return(Inf)
  v0 <- 0
  iter <- 0L
  repeat {
    iter <- iter+1L
    Diff <- RHS - mean(log(1+n*v0))
    Deriv <- mean( n/(1+n*v0) )
    Step <- Diff/Deriv
    v0 <- v0 + Step
    if(trace) cat("iter=",iter, " v0=",v0," Step=",Step,"\n")
    if(Step < eps) break
    if(iter >= 20L) {
      warning("Iteration limit reached")
      break
    }
  }
  1/v0
}
