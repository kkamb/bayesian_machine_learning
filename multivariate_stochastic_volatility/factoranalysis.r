##factoranalysis.r
##adapted from mike west's latent factor analysis code

n=length(Yfactor[,1])
M<-colno # Number of manifest variables: instruments
tau<-rep(1,M) # Error precision
beta<-matrix(0,M,K)
diag(beta)<-1 # Ensure first K diagonal elements > 0
beta.nrow=M*K

#### Store #######
Beta<-matrix(0,nsim,beta.nrow) # Loadings
Beta[1,]<-c(t(beta)) # Starting Values
Sigma<-matrix(0,nsim,M) # Error Variances
Sigma[1,]<-1/tau # Starting Values
F1<-F2<-F3<-matrix(0,nsim,n)
###### MCMC ########
for (i in 2:nsim){
# Update F
varf<-solve(diag(K)+t(beta)%*%diag(tau)%*%beta)
mf<-varf%*%t(beta)%*%diag(tau)%*%t(Yfactor)
v<-t(chol(varf))
z<-matrix(rnorm(n*K),K,n)
F<-t(v%*%z+mf)
# Update tau_j
for (j in 1:M){
d<-crossprod(Yfactor[,j]-F%*%beta[j,])
tau[j]<-rgamma(1,a+n/2,(b+d)/2)
}
# Update beta[1,1]
v<-1/(T0+tau[1]*crossprod(F[,1]))
m<-v*(T0*mu0+tau[1]*crossprod(F[,1],Yfactor[,1]))
beta[1,1]<-rtnorm(1,m,sqrt(v))
# Update rows 2:K of beta
for (j in 2:K) {
v<-solve(T0+tau[j]*crossprod(F[,1:j]))
m<-v%*%(T0*rep(mu0,j)+tau[j]*crossprod(F[,1:j],Yfactor[,j]))
beta[j,1:j]<-rtmvnorm(1,c(m),sigma=v,lower=c(rep(-Inf,j-1),0))
}
# Update remaining rows of beta
for (j in (K+1):M) {
v<-solve(T0+tau[j]*crossprod(F))
m<-v%*%(T0*rep(mu0,K)+tau[j]*crossprod(F,Yfactor[,j]))
beta[j,]<-rmvnorm(1,m,v)
}
# Store
Beta[i,]<-c(t(beta))
Sigma[i,]<-1/tau
F1[i,]<-F[,1]
F2[i,]<-F[,2]
F3[i,]<-F[,3]
}
##### Results#####
end.pt=length(Beta[,1])
beg.pt=10000
mbeta<-apply(Beta[beg.pt:end.pt,],2,mean)
Mbeta<-matrix(mbeta,M,K,byrow=T) #matrix of the factor loadings
f1=colMeans(F1[beg.pt:nsim,]) #f1 is the vector that contains the first factor for all times
f2=colMeans(F2[beg.pt:nsim,]) #f2 is the vector containing the second factor for all times
f3=colMeans(F3[beg.pt:nsim,]) #f3 is the vector containing the third factor for all times
Fmat=rbind(f1,f2,f3) #matrix of 3 factors
msig<-apply(Sigma[beg.pt:nsim,],2,mean)
Betas.uncertainty=matrix(NA,nrow=beta.nrow,ncol=3)
for (i in 1:beta.nrow){
Betas.uncertainty[i,1]=quantile(Beta[beg.pt:nsim,i],c(.025,.5,.975))[1]
Betas.uncertainty[i,2]=quantile(Beta[beg.pt:nsim,i],c(.025,.5,.975))[2]
Betas.uncertainty[i,3]=quantile(Beta[beg.pt:nsim,i],c(.025,.5,.975))[3]
}
#DIAGNOSTICS
Beta.mcmc=mcmc(Beta[beg.pt:end.pt,])
Sigma.mcmc=mcmc(Sigma[beg.pt:end.pt,])
pdf("betaloadings.pdf")
par(mfrow=c(5,3))
traceplot(Beta.mcmc[,1],main="Factor1 S&P Loadings")
traceplot(Beta.mcmc[,4],main="Factor1 10-Year Loadings")
traceplot(Beta.mcmc[,5],main="Factor2 10-Year Loadings")
traceplot(Beta.mcmc[,7],main="Factor1 Oil Loadings")
traceplot(Beta.mcmc[,9],main="Factor3 Oil Loadings")
traceplot(Beta.mcmc[,10],main="Factor1 Bond Loadings")
traceplot(Beta.mcmc[,11],main="Factor2 Bond Loadings")
traceplot(Beta.mcmc[,13],main="Factor1 Emerging Market Loadings")
traceplot(Beta.mcmc[,14],main="Factor2 Emerging Market Loadings")
traceplot(Beta.mcmc[,15],main="Factor3 Emerging Market Loadings")
traceplot(Beta.mcmc[,16],main="Factor1 JPN Yen Loadings")
traceplot(Beta.mcmc[,17],main="Factor2 JPN Yen Loadings")
traceplot(Beta.mcmc[,18],main="Factor3 JPN Yen Loadings")
traceplot(Beta.mcmc[,20],main="Factor2 Gold Loadings")
traceplot(Beta.mcmc[,21],main="Factor3 Gold Loadings")
dev.off()
pdf("factorsigmas.pdf")
par(mfrow=c(4,2))
traceplot(Sigma.mcmc[,1],main="S&P Idiosyncratic Variance")
traceplot(Sigma.mcmc[,2],main="10 Year Idiosyncratic Variance")
traceplot(Sigma.mcmc[,3],main="Oil Commodity Idiosyncratic Variance")
traceplot(Sigma.mcmc[,4],main="Bond Market Idiosyncratic Variance")
traceplot(Sigma.mcmc[,5],main="Emerging Markets Idiosyncratic Variance")
traceplot(Sigma.mcmc[,6],main="JPN Yen Idiosyncratic Variance")
traceplot(Sigma.mcmc[,7],main="Gold Futures Idiosyncratic Variance")
