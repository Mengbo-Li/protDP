# get the design matrix from observed mu without groups
defaultDesign <- function(mu) {
  model.matrix(~ 1 + mu, data = data.frame(1, mu))
}

# calculate pi's
logistic.p <- function(X, beta) {
  # X = design matrix, Nx2, only intercept and mu_obs
  # beta = 2x1, beta0 and beta1
  X <- as.matrix(X)
  beta <- as.vector(beta)
  p <- plogis(X %*% beta)
  return(p)
}

# negative log-likelihood of zero-truncated binomial distribution
ZBloglik <- function(y, w, p) {
  # input: y = detected counts = rowSums(!is.na(Y))
  # w = samples size = ncol(Y) 
  # p = probabilities
  # sum(dztbinom(y, size = w, p = p, log = TRUE))
  dens <- dbinom(y, size = w, p = p, log = TRUE) - 
    pbinom(0.5, size = w, p = p, log.p = TRUE, lower.tail = FALSE)
  -sum(dens)
}

# incomplete beta function B(z, a, b)
# this is equivalent to actuar:::betaint(z, a, b)/gamma(a+b)
iBeta <- function(z, a, b) {
  pbeta(z, a, b) * beta(a, b)
}


# Newton-Raphson for betas given alpha
NR_ZBLL_betas <- function(X, dp, wt, betaStart, alpha = 0.9,
                          eps = 1e-4, maxit = 100) {

  # calculate X matrix
  # X <- defaultDesign(mu_obs)
  y <- dp*wt
  
  betas <- betaStart
  names(betas) <- paste("b", 0:(length(betas)-1), sep = "")
  loglikStart <- ZBloglik(y, wt, alpha*logistic.p(X, betas))
  loglik <- loglikStart
  logliks <- loglik
  diff.loglik <- -Inf
  
  ## initialise a few output values
  hist <- data.frame(i = 0, alpha = alpha, 
                     matrix(betas, nrow = 1), 
                     loglik = loglik, 
                     diff.loglik, 
                     step.size = NA)
  colnames(hist)[3:(2+length(betas))] <- names(betas)
  ## line search step
  search.steps <- setdiff(seq(-2, 2, by = 0.1), 0)
  i <- 1
  
  while ((i <= maxit) & (diff.loglik < 0) & (-diff.loglik > eps)) {
    # calculate increment
    ## first-order derivative: Dloglik
    p <- as.numeric(logistic.p(X, betas))
    mu <- alpha * p
    the1minusMu <- 1-mu
    iBeta_mu_i <- iBeta(mu, 1, wt)
    dlogIpi <- the1minusMu^(wt-1) / iBeta_mu_i # D_i: first order derivative for logI_mu_i wrt mu_i
    Dloglik <- t(X) %*% as.matrix((wt*(dp-mu)/the1minusMu - mu*dlogIpi)*(1-p))
    
    ## negative of second-order derivative: nDDloglik
    the1m2mu <- 1-2*mu
    the1mp1m2p <- (1-p)*(1-2*p)
    vi <- -wt * (mu^2+the1m2mu*dp) *(1-p)^2 / the1minusMu^2 + 
      wt*(dp-mu)*the1mp1m2p / mu / the1minusMu + 
      ((wt-1)*dlogIpi/the1minusMu + dlogIpi^2) * mu^2 * (1-p)^2 - 
      dlogIpi*mu*the1mp1m2p
    V <- diag(-vi)
    nDDloglik <- t(X) %*% V %*% X
    # solve for increment
    increm <- solve(nDDloglik, Dloglik)
    
    # line search for improved step size
    loglik.searches <- rep(NA, length(search.steps))
    for (this_step in 1:length(search.steps)) {
      loglik.searches[this_step] <- ZBloglik(y, wt, alpha*logistic.p(X, betas+increm*search.steps[this_step]))
    }
    # step size index for maximum decrease in negative log-likelihood (if tie, use the first appeared)
    ind.best.step <- which(loglik.searches == min(loglik.searches, na.rm = TRUE))[1]
    
    # update beta
    betaNew <- as.numeric(betas + increm * search.steps[ind.best.step])
    loglikNew <- ZBloglik(y, wt, alpha*logistic.p(X, betaNew))
    diff.loglik <- loglikNew - loglik
    
    # iteration history
    hist <- rbind(hist, 
                  c(i, alpha = alpha, 
                    matrix(betaNew, nrow = 1), 
                    loglik = loglikNew, 
                    diff.loglik, 
                    search.steps[ind.best.step]))
    betas <- as.vector(betaNew)
    loglik <- loglikNew
    i <- i + 1
  }
  # results
  names(betas) <- paste("b", 0:(length(betas)-1), sep = "")
  list(betas = betas, negLL = loglik, hist = hist)
}



# update alpha
cappedLogitLinear_ZBNR <- function(X, dp, wt, alpha0 = 0.9, beta0, 
                                   maxit = 100, trace = TRUE, eps = 1e-04) {
  alpha <- alpha0
  alpha.hist <- alpha
  
  beta_byNR <- NR_ZBLL_betas(X, dp, wt, beta0, alpha = alpha0,
                             eps = 1e-4, maxit = 100)
  betas <- beta_byNR$betas
  betas.hist <- betas
  negLL <- beta_byNR$negLL
  negLL.hist <- negLL
  
  for (i in 1:maxit) {
    p <- as.numeric(logistic.p(X, betas))
    mu <- alpha * p
    newAlpha <- sum(wt * dp/(1-mu))/sum(wt * p/(1-mu) + 
                                            p * (1-mu)^(wt-1)/iBeta(mu, 1, wt))
    if (newAlpha >= 1) {
      alpha <- 1
      if (trace) 
        cat(paste("alpha_", i, sep = ""), "=", alpha, "\n")
      beta_byNR <- NR_ZBLL_betas(X, dp, wt, betas, alpha = alpha,
                                 eps = 1e-4, maxit = 100)
      newBetas <- beta_byNR$beta
      newNegLL <- beta_byNR$negLL
      alpha.hist <- c(alpha.hist, alpha)
      betas.hist <- rbind(betas.hist, newBetas)
      negLL.hist <- c(negLL.hist, newNegLL)
      betas <- newBetas
      break
    }
    if (abs(newAlpha - alpha) < eps) break
    if (trace) 
      cat(paste("alpha_", i, sep = ""), "=", newAlpha, "\n")
    beta_byNR <- NR_ZBLL_betas(X, dp, wt, betas, alpha = newAlpha,
                               eps = 1e-4, maxit = 100)
    newBetas <- beta_byNR$beta
    newNegLL <- beta_byNR$negLL
    alpha.hist <- c(alpha.hist, newAlpha)
    betas.hist <- rbind(betas.hist, newBetas)
    negLL.hist <- c(negLL.hist, newNegLL)
    alpha <- newAlpha
    betas <- newBetas
  }
  info <- cbind(alpha.hist, betas.hist, negLL.hist)
  colnames(info) <- c("alpha", names(betas), "neg.ZBLL")
  rownames(info) <- paste("i=", 0:(nrow(info) - 1), sep = "")
  final.params <- c(alpha, betas)
  names(final.params) <- c("alpha", names(betas))
  list(params = final.params, info = info)
}













