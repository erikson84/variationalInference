library(gtools)
source('varBin.R')
source('gibbsBin.R')

set.seed(300)
# Set the number of clusters
# The variational algorithm is pretty fast
# even with hundred of clusters. But the Gibbs
# sampler is very slow, so beware!
K <- 3
phi <- rdirichlet(1, rep(1/K, K))
theta <- rbeta(K, 1, 1)
# Assume different number of trials for each datapoint
Nt <- sample(10:30, 2000, replace=T)
Z <- sample(1:K, 2000, replace=T, phi)
y <- rbinom(2000, Nt, theta[Z])

# Use the variational algorith to approximate the posterior
teste <- varMixBern(y, Nt=Nt, K=3)
apply(teste$alphaBetaVar, 1, function(x) x[1]/sum(x))
teste$etaVar/sum(teste$etaVar)
teste$ELBO

# Use the Gibbs sampling algorithm to approximate the posterior
testeGibbs <- gibbsMixBern(y, Nt=Nt, K=3)
apply(testeGibbs$phi, 2, mean)
apply(testeGibbs$theta, 2, mean)
