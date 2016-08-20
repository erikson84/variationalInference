library(gtools)
source('varBin.R')
source('gibbsBin.R')

set.seed(2010)
# Set the number of clusters
# The variational algorithm is pretty fast
# even with hundred of clusters. But the Gibbs
# sampler is very slow, so beware!
N <- 2000
K <- 3
phi <- rdirichlet(1, rep(1/K, K))
theta <- rbeta(K, 1, 1)
# Assume different number of trials for each datapoint
Nt <- sample(10:30, N, replace=T)
Z <- sample(1:K, N, replace=T, phi)
y <- rbinom(N, Nt, theta[Z])

# Use the variational algorith to approximate the posterior
testVar <- varMixBern(y, Nt=Nt, K=3)
apply(testVar$alphaBetaVar, 1, function(x) x[1]/sum(x))
testVar$etaVar/sum(testVar$etaVar)
testVar$ELBO

# Use the Gibbs sampling algorithm to approximate the posterior
testGibbs <- gibbsMixBern(y, Nt=Nt, K=3)
apply(testGibbs$phi, 2, mean)
apply(testGibbs$theta, 2, mean)

par(mfrow=c(2, 3))
plot(density(testGibbs$theta[,1], bw=.005), xlab='theta[1]', main='Posterior for theta[1]', ylim=c(0, 240))
curve(dbeta(x, testVar$alphaBetaVar[3, 1], testVar$alphaBetaVar[3, 2]), add=T, col='blue')
abline(v=theta[1], lty=2)

plot(density(testGibbs$theta[,2], bw=.005), xlab='theta[2]', main='Posterior for theta[2]', ylim=c(0, 240))
curve(dbeta(x, testVar$alphaBetaVar[1, 1], testVar$alphaBetaVar[1, 2]), add=T, col='blue')
abline(v=theta[3], lty=2)

plot(density(testGibbs$theta[,3], bw=.005), xlab='theta[3]', main='Posterior for theta[3]', ylim=c(0, 240))
curve(dbeta(x, testVar$alphaBetaVar[2, 1], testVar$alphaBetaVar[2, 2]), add=T, col='blue')
abline(v=theta[2], lty=2)

plot(density(testGibbs$phi[,1], bw=.005), xlab='phi[1]', main='Posterior for phi[1]', ylim=c(0, 60))
curve(dbeta(x, testVar$etaVar[3], sum(testVar$etaVar[c(2, 1)])), add=T, col='blue')
abline(v=phi[1], lty=2)

plot(density(testGibbs$phi[,2], bw=.005), xlab='phi[2]', main='Posterior for phi[2]', ylim=c(0, 60))
curve(dbeta(x, testVar$etaVar[1], sum(testVar$etaVar[c(2, 3)])), add=T, col='blue')
abline(v=phi[3], lty=2)

plot(density(testGibbs$phi[,3], bw=.005), xlab='phi[3]', main='Posterior for phi[3]', ylim=c(0, 60))
curve(dbeta(x, testVar$etaVar[2], sum(testVar$etaVar[c(1, 3)])), add=T, col='blue')
abline(v=phi[2], lty=2)
