

impCov <- function(fit) {

  # basic info
  nfactors <- fit$factors

  # number of variables
  p <- nrow(fit$loadings)

  # factor correlations
  # psi <- fit$r.scores # correlations between factors
  if(nfactors==1) psi <- matrix(1,1,1) else psi <- fit$Phi

  # residual variances
  theta <- matrix(0, p, p)
  diag(theta) <- diag(fit$residual)

  # factor loadings
  load <- matrix(0, p, nfactors)
  for(i in 1:nfactors) load[, i] <- fit$loadings[, i]

  implied_cov <- load %*% psi %*% t(load) + theta
  implied_cor <- cov2cor(implied_cov)

  return(implied_cor)

} # eoF

