gibbsMixBern <- function(y, Nt, K, alpha=1, beta=1, eta=1/K, warmup=2000, iters=5000){
  N <- length(y)  
  phi <- rdirichlet(1, rep(eta, K))
  theta <- rbeta(K, alpha, beta)
  Z <- sample(1:K, N, replace = T, prob=phi)
  out <- list(phi=phi, theta=theta, Z=Z)
  
  for (i in 1:iters){
    tableZ <- rep(NA, K)
    for (k in 1:K){
      tableZ[k] <- sum(Z==k)
    }
    phi <- rdirichlet(1, eta+tableZ)
    Z <- rep(NA, N)
    for (z in 1:N){
      Z[z] <- sample(1:K, 1, prob=phi*dbinom(y[z], Nt, theta))
    }
    theta <- rep(NA, K)
    for (k in 1:K){
      theta[k] <- rbeta(1, alpha + sum(y[Z==k]), beta + sum(Nt - y[Z==k]))
    }
    if (i > warmup){
      if (((i/iters) %% 0.1) ==0){
        cat((i/iters) * 100)
        cat('\n')
      }
      out <- list(phi=rbind(out$phi, phi), theta=rbind(out$theta, theta), Z=rbind(out$Z, Z))
    }
  }
  out <- list(phi=out$phi[-1,], theta=out$theta[-1,], Z=out$Z[-1,])
  out
}
