##ffbsandgibbs.r file
##adapted from heibert lopes' stochastic volatility R code

ffbsu = function(y,F,alpha,V,mu,phi,tau,m0,C0){
  n = length(y)
  if (length(F)==1){F = rep(F,n)}
  if (length(alpha)==1){alpha = rep(alpha,n)}
  if (length(V)==1){V = rep(V,n)}
  a = rep(0,n)
  R = rep(0,n)
  m = rep(0,n)
  C = rep(0,n)
  B = rep(0,n-1)
  H = rep(0,n-1)
  a1 = mu + phi*m0
  R1 = phi^2*C0 + tau
  # time t=1
  a[1] = a1
  R[1] = R1
  f    = alpha[1]+F[1]*a[1]
  Q    = R[1]*F[1]**2+V[1]
  A    = R[1]*F[1]/Q
  m[1] = a[1]+A*(y[1]-f)
  C[1] = R[1]-Q*A**2
  # forward filtering
  for (t in 2:n){
    a[t] = mu + phi*m[t-1]
    R[t] = C[t-1]*phi**2 + tau
    f    = alpha[t]+F[t]*a[t]
    Q    = R[t]*F[t]**2+V[t]
    A    = R[t]*F[t]/Q
    m[t] = a[t]+A*(y[t]-f)
    C[t] = R[t]-Q*A**2
    B[t-1] = C[t-1]*phi/R[t]
    H[t-1] = sqrt(C[t-1]-R[t]*B[t-1]**2)
  }
  theta = rep(0,n)
  # backward sampling
  theta[n] = rnorm(1,m[n],sqrt(C[n]))
  for (t in (n-1):1)
    theta[t] = rnorm(1,m[t]+B[t]*(theta[t+1]-a[t+1]),H[t])
  return(theta)
}

gibbs.step = function(y,X,theta0,A,nu0,s0){
  n     = length(y)
  k     = ncol(X)
  par1  = (nu0+n)/2
  var   = solve(crossprod(X,X)+A)
  mean  = matrix(var%*%(crossprod(X,y)+crossprod(A,theta0)),k,1)
  par2  = nu0*s0 + sum((y-crossprod(t(X),mean))^2)
  par2  = (par2 + crossprod(t(crossprod(mean-theta0,A)),mean-theta0))/2
  sig2  = 1/rgamma(1,par1,par2)
  var   = var*sig2
  mean  = mean + crossprod(t(chol(var)),rnorm(2))
  return(c(mean,sig2))
}

#----------------------------------------------------------------------------
# Sample Z from 1,2,...,k, with P(Z=i) proportional to q[i]N(mu[i],sig2[i])
#----------------------------------------------------------------------------
ncind = function(y,mu,sig,q){
  w = dnorm(y,mu,sig)*q
  return(sample(1:length(q),size=1,prob=w/sum(w)))
}

#------------------------------------------------
# Sampling the log-volatilities in a 
# standard univariate stochastic volatility model
#------------------------------------------------
svu = function(logy2,h,mu,phi,tau,m0,C0){
  mu.est   = c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859)
  sig2.est = c(5.795960,2.613690,5.179500,0.167350,0.640090,0.340230,1.262610)
  q.prob    = c(0.007300,0.105560,0.000020,0.043950,0.340010,0.245660,0.257500)
  sig.est  = sqrt(sig2.est)
  z    = sapply(logy2-h,ncind,mu.est,sig.est,q.prob)
  ffbsu(logy2,1.0,mu.est[z],sig2.est[z],mu,phi,tau,m0,C0)
}

ps=MU=PHI=TAU=NULL
h1s=h2s=h3s=h4s=h5s=h6s=h7s=NULL
fh1s=fh2s=fh3s=NULL
X=NULL

for (iter in 1:S){

for(icols in 1:ncol){
  h[,icols]=svu(ymat[,icols],h[,icols],mu[icols],phi[icols],tau2[icols],m0,C0) #ffbs step to update implicit volatility vector h
  var=1/(iC0+phi[icols]^2/tau2[icols])
  mean=var*(iC0m0+phi[icols]*(h[1,icols]-mu[icols])/tau2[icols])
  h0=rnorm(ncol,mean,sqrt(var))
  X=cbind(1,c(h0[icols],h[1:(n-1),icols]))
  par=gibbs.step(h[1:n,icols],X,theta0,iV0,nu0,s0) #gibbs sampler step
  mu[icols]=par[1]
  phi[icols]=par[2]
  tau[icols]=par[3]
}
  MU=rbind(MU,mu)
  PHI=rbind(PHI,phi)
  TAU=rbind(TAU,tau)
h1s=rbind(h1s,h[,1]);h2s=rbind(h2s,h[,2]);h3s=rbind(h3s,h[,3]);h4s=rbind(h4s,h[,4])
h5s=rbind(h5s,h[,3]);h6s=rbind(h6s,h[,6]);h7s=rbind(h7s,h[,7])
fh1s=rbind(fh1s,h[,8]);fh2s=rbind(fh2s,h[,9]);fh3s=rbind(fh3s,h[,10])
}

end.pt=length(MU[,1])
beg.pt=end.pt-2000

hrows=length(h1s[1,])
hmat.mean=hmat.upper=hmat.lower=matrix(NA,nrow=hrows,ncol=7)
fhmat.mean=fhmat.upper=fhmat.lower=matrix(NA,nrow=hrows,ncol=3)

for (i in 1:hrows){
hmat.lower[i,1]=quantile(h1s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
hmat.mean[i,1]=quantile(h1s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
hmat.upper[i,1]=quantile(h1s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
hmat.lower[i,2]=quantile(h2s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
hmat.mean[i,2]=quantile(h2s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
hmat.upper[i,2]=quantile(h2s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
hmat.lower[i,3]=quantile(h3s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
hmat.mean[i,3]=quantile(h3s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
hmat.upper[i,3]=quantile(h3s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
hmat.lower[i,4]=quantile(h4s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
hmat.mean[i,4]=quantile(h4s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
hmat.upper[i,4]=quantile(h4s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
hmat.lower[i,5]=quantile(h5s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
hmat.mean[i,5]=quantile(h5s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
hmat.upper[i,5]=quantile(h5s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
hmat.lower[i,6]=quantile(h6s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
hmat.mean[i,6]=quantile(h6s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
hmat.upper[i,6]=quantile(h6s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
hmat.lower[i,7]=quantile(h7s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
hmat.mean[i,7]=quantile(h7s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
hmat.upper[i,7]=quantile(h7s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
fhmat.lower[i,1]=quantile(fh1s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
fhmat.mean[i,1]=quantile(fh1s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
fhmat.upper[i,1]=quantile(fh1s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
fhmat.lower[i,2]=quantile(fh2s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
fhmat.mean[i,2]=quantile(fh2s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
fhmat.upper[i,2]=quantile(fh2s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
fhmat.lower[i,3]=quantile(fh3s[beg.pt:end.pt,i],c(.025,.5,.975))[1]
fhmat.mean[i,3]=quantile(fh3s[beg.pt:end.pt,i],c(.025,.5,.975))[2]
fhmat.upper[i,3]=quantile(fh3s[beg.pt:end.pt,i],c(.025,.5,.975))[3]
}

timeperiod=c(1:707)
pdf("f1vol.pdf")
plot(timeperiod,fhmat.mean[,1],type="n",main="Equities Factor Volatilities")
lines(fhmat.mean[,1])
lines(fhmat.lower[,1],lty=2)
lines(fhmat.upper[,1],lty=2)
dev.off()

pdf("f2vol.pdf")
plot(timeperiod,fhmat.mean[,2],type="n",main="Fixed Income Factor Volatilities")
lines(fhmat.mean[,2])
lines(fhmat.lower[,2],lty=2)
lines(fhmat.upper[,2],lty=2)
dev.off()

pdf("f3vol.pdf")
plot(timeperiod,fhmat.mean[,3],type="n",main="Commodities Factor Volatilities")
lines(fhmat.mean[,3])
lines(fhmat.lower[,3],lty=2)
lines(fhmat.upper[,3],lty=2)
dev.off()

pdf("idiosyncraticvols.pdf")
par(mfrow=c(2,2))
plot(timeperiod,hmat.mean[,1],type="n",main="SPX Idiosyncratic Volatilities")
lines(hmat.mean[,1])
lines(hmat.lower[,1],lty=2)
lines(hmat.upper[,1],lty=2)

plot(timeperiod,hmat.mean[,2],type="n",main="10yr Idiosyncratic Volatilities")
lines(hmat.mean[,2])
lines(hmat.lower[,2],lty=2)
lines(hmat.upper[,2],lty=2)

plot(timeperiod,hmat.mean[,3],type="n",main="Oil Idiosyncratic Volatilities")
lines(hmat.mean[,3])
lines(hmat.lower[,3],lty=2)
lines(hmat.upper[,3],lty=2)

plot(timeperiod,hmat.mean[,4],type="n",main="Bond Idiosyncratic Volatilities")
lines(hmat.mean[,4])
lines(hmat.lower[,4],lty=2)
lines(hmat.upper[,4],lty=2)

pdf("idiosyncraticvols2.pdf")
par(mfrow=c(2,2))
plot(timeperiod,hmat.mean[,5],type="n",main="MXEF Idiosyncratic Volatilities")
lines(hmat.mean[,5])
lines(hmat.lower[,5],lty=2)
lines(hmat.upper[,5],lty=2)

plot(timeperiod,hmat.mean[,6],type="n",main="Yen Idiosyncratic Volatilities")
lines(hmat.mean[,6])
lines(hmat.lower[,6],lty=2)
lines(hmat.upper[,6],lty=2)

plot(timeperiod,hmat.mean[,7],type="n",main="Gold Idiosyncratic Volatilities")
lines(hmat.mean[,7])
lines(hmat.lower[,7],lty=2)
lines(hmat.upper[,7],lty=2)
dev.off()
