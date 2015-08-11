library(MASS);library(coda);library(xtable);library(car);library(mvtnorm);library(msm);library(moments);library(tmvtnorm)

#################data input and cleaning##################
findata=read.table("weeklydata.csv",header=T,sep=",")
attach(findata)
y=findata[2:8]
y=y[1:708,1:7] #get rid of NAs
xno=length(y[,1])
colno=length(y)
ylog=log(y) #ylog are log of prices
yret=ylog
yret[1,]=0
j=2
for (j in 2:xno){
    for (i in 1:colno){
yret[j,i]=ylog[j,i]-ylog[j-1,i]
if (yret[j,i]==0) yret[j,i]=.0000001
}
}
yret[,2]=-8*yret[,2] #since this is the US10 year, quoted in yields
yret=yret[2:xno,]
#############################################################

yvol=findata$SPX.Volume
yvol=yvol[1:708]
ma <- function(x,n){filter(x,rep(1/n,n), sides=1)}
yvol.50avg=ma(yvol,50)
yvol=yvol[51:708]
yvol.50avg=yvol.50avg[51:708]
yvol.X=yvol/yvol.50avg #creating a yvol.X factor, that normalizes the S&P volume
yret=yret[51:708,]

########MCMC ALGORITHM###############
set.seed(12)

#factor analysis inputs#######
Yfactor=yret
Yfactor=scale(Yfactor)
nsim<-20000
K<-3 #no of latent factors
mu0<-0	# Prior mean of beta_ij 
C0<-1 	# Prior variance of beta_ij (note: beta_ii (i=1,..,k) trunc. >0)
T0<-1/C0     # Prior precision
a<-b<-.01	# Hyperparameters for Sigma
source("factoranalysis.r")
####needs above inputs, outputs Mbeta, matrix of factor loadings, and Fmat, matrix of factors

y.ind=Yfactor-t(Mbeta%*%Fmat) #taking out the factors to get the idiosyncratic volatilities

########Analyzing data################
pdf("preandpostcorr.pdf")
par(mfrow=c(1,2))
z=cor(Yfactor)
require(lattice)
levelplot(z)
z=cor(y.ind)
require(lattice)
levelplot(z)
dev.off()

#########ffbs and gibbs step#####################

y.ind=log(y.ind^2)
F.logsquare=t(log(Fmat^2)) #factor volatities
ymat=cbind(y.ind,F.logsquare) #creates a matrix of both independent series and factor volatities, whose implicit volatities we are trying to derive
n=length(ymat[,1])
ncol=length(ymat[1,])

####priors for theta###############################
mu=rep(0,ncol)
phi=rep(.99,ncol)
theta0=c(mu[1],phi[1]) #prior beta vector
V0=diag(100,2) #prior covariance
iV0=solve(V0)

####priors for sigma (inverse gamma distribution)##
nu0=10
s0=.02
tau2=1/rgamma(ncol,nu0/2,nu0/2*s0)
tau=sqrt(tau2)

####priors for latent volatility vector h##########
m0=0 #h is normally distributed with prior parameters m0 and c0
C0=100
iC0=1/C0
iC0m0=iC0*m0
h0=rep(0,ncol) #setting initial value to zero
h=matrix(0,nrow=n,ncol=ncol)
h[1,]=rnorm(ncol,mu+phi*(h0-mu),tau)
for (t in 2:n) #creating prior volatilities for each time period t
  h[t,] = rnorm(ncol,mu+phi*(h[t-1,]-mu),tau)

S=4000 #number of iterations
source("ffbsandgibbs.r")
####outputs MU,PHI,TAU2, h1s to h7s, fh1s to fh3s

#DIAGNOSTICS
end.pt=length(MU[,1])
beg.pt=end.pt-2000

mu.mcmc=mcmc(MU[beg.pt:end.pt,])
phi.mcmc=mcmc(PHI[beg.pt:end.pt,])
tau.mcmc=mcmc(TAU[beg.pt:end.pt,])
h1.mcmc=mcmc(h1s[beg.pt:end.pt,])
fh1.mcmc=mcmc(fh1s[beg.pt:end.pt,])
fh2.mcmc=mcmc(fh2s[beg.pt:end.pt,])
fh3.mcmc=mcmc(fh3s[beg.pt:end.pt,])

pdf("idiosyncraticfactor1.pdf")
par(mfrow=c(3,3))
traceplot(mu.mcmc[,1],main="S&P Idiosyncratic Volatility mu")
traceplot(phi.mcmc[,1],main="S&P Idiosyncratic Volatility phi")
traceplot(tau.mcmc[,1],main="S&P Idiosyncratic Volatility tau")
traceplot(mu.mcmc[,2],main="10yr Idiosyncratic Volatility mu")
traceplot(phi.mcmc[,2],main="10yr Idiosyncratic Volatility phi")
traceplot(tau.mcmc[,2],main="10yr Idiosyncratic Volatility tau")
traceplot(mu.mcmc[,3],main="Oil Idiosyncratic Volatility mu")
traceplot(phi.mcmc[,3],main="Oil Idiosyncratic Volatility phi")
traceplot(tau.mcmc[,3],main="Oil Idiosyncratic Volatility tau")
dev.off()

pdf("idiosyncraticfactor2.pdf")
par(mfrow=c(4,3))
traceplot(mu.mcmc[,4],main="Bond Idiosyncratic Volatility mu")
traceplot(phi.mcmc[,4],main="Bond Idiosyncratic Volatility phi")
traceplot(tau.mcmc[,4],main="Bond Idiosyncratic Volatility tau")
traceplot(mu.mcmc[,5],main="MXEF Idiosyncratic Volatility mu")
traceplot(phi.mcmc[,5],main="MXEF Idiosyncratic Volatility phi")
traceplot(tau.mcmc[,5],main="MXEF Idiosyncratic Volatility tau")
traceplot(mu.mcmc[,6],main="Yen Idiosyncratic Volatility mu")
traceplot(phi.mcmc[,6],main="Yen Idiosyncratic Volatility phi")
traceplot(tau.mcmc[,6],main="Yen Idiosyncratic Volatility tau")
traceplot(mu.mcmc[,7],main="Gold Idiosyncratic Volatility mu")
traceplot(phi.mcmc[,7],main="Gold Idiosyncratic Volatility phi")
traceplot(tau.mcmc[,7],main="Gold Idiosyncratic Volatility tau")
dev.off()

pdf("commonfactors.pdf")
par(mfrow=c(3,3))
traceplot(mu.mcmc[,8],main="Equities Factor Volatility mu")
traceplot(phi.mcmc[,8],main="Equities Factor Volatility phi")
traceplot(tau.mcmc[,8],main="Equities Factor Volatility tau")
traceplot(mu.mcmc[,9],main="Bond Factor Volatility mu")
traceplot(phi.mcmc[,9],main="Bond Factor Volatility phi")
traceplot(tau.mcmc[,9],main="Bond Factor Volatility tau")
traceplot(mu.mcmc[,10],main="Commodities Factor Volatility mu")
traceplot(phi.mcmc[,10],main="Commodities Factor Volatility phi")
traceplot(tau.mcmc[,10],main="Commodities Factor Volatility tau")
dev.off()

fname="mcmc1state.Rda"
fname="mcmc1h.Rda"
mcmc1state=data.frame(MU,PHI,TAU)
mcmc1h=data.frame(fh1s,fh2s,hf3s,h1s,h2s,h3s,h4s,h5s,h6s,h7s)
save(mcmc1state,file="mcmc1state.Rda")
save(mcmc1h,file="mcmc1h.Rda")
