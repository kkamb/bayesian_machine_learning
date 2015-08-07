library(DPpackage)
data(galaxy)
galaxy <- data.frame(galaxy,speeds=galaxy$speed/1000)
attach(galaxy)
prior1 <- list(alpha=1,m1=rep(0,1),psiinv1=diag(0.5,1),nu1=4,tau1=1,tau2=100)
# Fixing alpha and m1
prior2 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.5,1)),nu1=4,nu2=4, tau1=1,tau2=100)
# Fixing only alpha
prior3 <-
  list(alpha=1,m2=rep(0,1),s2=diag(100000,1),psiinv2=solve(diag(0.5,1)),nu1=4,nu2=4,tau1=1,tau2=100)
#Everything is random
prior4 <- list(a0=2,b0=1,m2=rep(0,1),s2=diag(100000,1), psiinv2=solve(diag(0.5,1)),
               nu1=4,nu2=4,tau1=1,tau2=100)
fit1.1 <- DPdensity(y=speeds,prior=prior1,mcmc=mcmc,state=state,status=TRUE)
fit1.2 <- DPdensity(y=speeds,prior=prior2,mcmc=mcmc,state=state,status=TRUE)
fit1.3 <- DPdensity(y=speeds,prior=prior3,mcmc=mcmc,state=state,status=TRUE)
fit1.4 <- DPdensity(y=speeds,prior=prior4,mcmc=mcmc,state=state,status=TRUE)
plot(fit1.1,ask=FALSE)
plot(fit1.2,ask=FALSE)
plot(fit1.3,ask=FALSE)
plot(fit1.4,ask=FALSE)