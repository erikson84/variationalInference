source('varBin.R')
source('gibbsBin.R')

K <- 3
phi <- rdirichlet(1, rep(1/K, K))
theta <- rbeta(K, 1, 1)
Nt <- sample(10:30, 2000, replace=T)
Z <- sample(1:K, 2000, replace=T, phi)
y <- rbinom(2000, Nt, theta[Z])

teste <- varMixBern(y, Nt=Nt, K=3)
(apply(teste$alphaBetaVar, 1, function(x) x[1]/sum(x)))
#theta
teste$etaVar/sum(teste$etaVar)
#phi
teste$ELBO

testeGibbs <- gibbsMixBern(y, Nt=30, K=3)
apply(testeGibbs$phi, 2, mean)
apply(testeGibbs$theta, 2, mean)