n=dim(Y)[1]
p=dim(Y)[2]
nu0=p+2
Sigma=cov(Y)
S0=(nu0-p-1)*Sigma0
ybar<-apply(Y,2,mean)

set.seed(12)
THETA=SIGMA=NULL
for (s in 1:nsamples) {
    ##update theta
    Ln=solve(solve(L0)+n*solve(Sigma))
    mun=Ln%*%(solve(L0)%*%mu0+n*solve(Sigma)%*%ybar)
    theta=rmvnorm(1,mun,Ln)

    ##update Sigma
    Sn=S0+(t(Y)-c(theta))%*%t(t(Y)-c(theta))
    Sigma=solve(rwish(1,nu0+n,solve(Sn)))

    #save results
    THETA=rbind(THETA,theta)
    SIGMA=rbind(SIGMA,c(Sigma))
}
THETA2=THETA[1001:nsamples,1:2]
SIGMA2=SIGMA[1001:nsamples,1:4]
