# A simple, non-optimized Gibbs sampler for a mixture of binomials.
# Used mainly to check the results from the variational algorithm.
# As the variational model, assumes the total number of trials is known.

# Parameters:
# y: integer vector of successes
# Nt: integer scalar or vector of total number of trials
# alpha, beta: positive real scalars of prior hyperparameters for theta
# eta: positive real vector of prior hyperparameter for phi
# warmup: number of initian iterations to discard
# iters: total number of iterations

# Returns:
# A names list with:
# phi: vector of cluster proportions
# theta: binomial parameter for each cluster
# Z: cluster attribution for each datapoint

library(gtools) # Needed for Rdirichlet initialization
gibbsMixBern <- function(y, Nt, K, alpha=1, beta=1, eta=1/K, warmup=1000, iters=2000){
  N <- length(y)  
  if (length(Nt)==1){
    Nt <- rep(Nt, N)
  }
  # Random initialization
  phi <- rdirichlet(1, rep(eta, K))
  theta <- rbeta(K, alpha, beta)
  Z <- sample(1:K, N, replace = T, prob=phi)
  out <- list(phi=phi, theta=theta, Z=Z)
  
  for (i in 1:iters){
    if ((((i/iters)*100) %% 10) == 0) {cat(paste(' ', (i/iters)*100, '%', collapse='', sep=''))}
    # Uses a loop to generate count instead of 'table' function
    # in case one cluster remains with zero counts
    tableZ <- rep(NA, K)
    for (k in 1:K){
      tableZ[k] <- sum(Z==k)
    }
    
    # Update phi parameter based on full conditional
    phi <- rdirichlet(1, eta+tableZ)
    # Update Z cluster attribution vector
    Z <- rep(NA, N)
    for (z in 1:N){
      Z[z] <- sample(1:K, 1, prob=phi*dbinom(y[z], Nt[z], theta))
    }
    #Update theta parameters
    theta <- rep(NA, K)
    for (k in 1:K){
      theta[k] <- rbeta(1, alpha + sum(y[Z==k]), beta + sum(Nt[Z==k] - y[Z==k]))
    }
    # Write random samples to output
    if (i > warmup){
      out <- list(phi=rbind(out$phi, phi), theta=rbind(out$theta, theta), Z=rbind(out$Z, Z))
    }
  }
  out <- list(phi=out$phi[-1,], theta=out$theta[-1,], Z=out$Z[-1,])
  out
}
