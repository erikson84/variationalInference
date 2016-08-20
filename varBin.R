library(gtools) # For rdirichlet
varMixBern <- function(y, Nt, K, alpha=1, beta=1, eta=1/K, epsilon=0.0001,
                       init=rdirichlet(length(y), rep(eta, K))){
  # Some basic quantities
  N <- length(y)
  
  # Initialize g_nk, each data point responsability
  gnk <- init
  
  # Compute some quantities we need...
  Nk <- apply(gnk, 2, sum)
  
  # Define variational parameters
  etaVar <- eta + Nk
  alphaVar <- alpha + colSums(matrix(y, nrow=N, ncol=K)*gnk) 
  betaVar <- beta + colSums(matrix(Nt-y, nrow=N, ncol=K)*gnk)
  
  # Define some quantities needed to compute ELBO.
  lnPhi <- digamma(etaVar) - digamma(sum(etaVar))
  lnTheta <- digamma(alphaVar) - digamma(alphaVar + betaVar)
  ln1mTheta <- digamma(betaVar) - digamma(alphaVar + betaVar)
  
  # 
  ElogpYZT <- function(){
    sum(gnk * (y %*% t(lnTheta) + (Nt-y) %*% t(ln1mTheta) + lchoose(Nt, y)))
  }
  
  ELBO <- sum((eta-1) * lnPhi) + ((alpha - 1)*sum(lnTheta) + (beta - 1)*sum(ln1mTheta)) +
    sum(gnk * matrix(lnPhi, N, K, byrow=T)) + ElogpYZT() - sum(gnk*log(gnk)) -
    (sum((etaVar-1)*lnPhi) + lgamma(sum(etaVar)) - sum(lgamma(etaVar))) - 
    sum(((alphaVar-1)*lnTheta) + ((betaVar-1)*ln1mTheta) -
       lgamma(alphaVar) - lgamma(betaVar) + lgamma(alphaVar + betaVar))
  print(ELBO)
  
  newELBO <- ELBO+1
  iter <- 1
  while(abs(newELBO-ELBO)>epsilon) {
    gnk <- matrix(lnPhi, N, K, byrow=T) + y %*% t(lnTheta) +
      (Nt-y) %*% t(ln1mTheta) + matrix(lchoose(Nt, y), N, K)
    gnk <- exp(gnk)
    gnk <- t(apply(gnk, 1, function(x) x/sum(x)))
    
    

    Nk <- apply(gnk, 2, sum)
    
    etaVar <- eta + Nk
    alphaVar <- alpha + colSums(matrix(y, nrow=N, ncol=K)*gnk) 
    betaVar <- beta + colSums(matrix(Nt-y, nrow=N, ncol=K)*gnk)
    
    lnPhi <- digamma(etaVar) - digamma(sum(etaVar))
    lnTheta <- digamma(alphaVar) - digamma(alphaVar + betaVar)
    ln1mTheta <- digamma(betaVar) - digamma(alphaVar + betaVar) 
    
    ELBO <- newELBO
    newELBO <- sum((eta-1) * lnPhi) + ((alpha - 1)*sum(lnTheta) + (beta - 1)*sum(ln1mTheta)) +
      sum(gnk * matrix(lnPhi, N, K, byrow=T)) + ElogpYZT() - sum(gnk*log(gnk)) -
      (sum((etaVar-1)*lnPhi) + lgamma(sum(etaVar)) - sum(lgamma(etaVar))) - 
      sum(((alphaVar-1)*lnTheta) + ((betaVar-1)*ln1mTheta) -
            lgamma(alphaVar) - lgamma(betaVar) + lgamma(alphaVar + betaVar))
    print(iter)
    iter <- iter+1
    print(newELBO)
    
  }
  
  list(etaVar=etaVar, alphaBetaVar=cbind(alphaVar, betaVar), catVar=gnk, ELBO=newELBO)
  
}

