library(mvtnorm);library(MASS);library(coda);library(hdrcde)


mu0=c(0,0)
Sigmatrue=matrix(c(1,.8,.8,1),nrow=2)
L0=matrix(c(.5,0,0,.5),nrow=2)
Sigma0 = matrix(c(1,0.5,0.5,1),nrow=2)
n=100
pmiss=c(.05,.1,.15,.2,.25,.3,.35,.4,.45,.5)
x.prob=length(pmiss)
nsamples=3000

setClass(Class="discardMat",
         representation(
            Y="matrix",
            index="vector"
          )
)

rwish<-function(n,nu0,S0)
    {
        sS0<-chol(S0)
        S<-array( dim=c(dim(S0),n))
        for (i in 1:n)
            {
                Z<-matrix(rnorm(nu0*dim(S0)[1]),nu0,dim(S0)[1])%*% sS0
                S[,,i]<- t(Z)%*%Z
            }
        S[,,1:n]
    }
discard<-function(n,ratio)
    {
    n=100
    discard.vec<-sample(n*2,n*2*ratio)
    set.seed(12)
    Y=rmvnorm(n,mu0,Sigmatrue)
    Y[discard.vec]<-NA
    return(new("discardMat",
          Y=Y,
          index=discard.vec))
    }

####full data analysis
set.seed(12)
Y = rmvnorm(n, mu0, Sigmatrue)
source("gibbsrun.r")
full.data.theta=THETA2
full.data.sigma=SIGMA2

theta.mcmc=mcmc(full.data.theta)
sigma.mcmc=mcmc(full.data.sigma)
summary(theta.mcmc)
summary(sigma.mcmc)

plot(theta.mcmc)
plot(sigma.mcmc)
autocorr.plot(theta.mcmc)
autocorr.plot(sigma.mcmc)

####complete case analysis
theta.comp.array=array(data=NA,dim=c(nsamples-1000,2,x.prob))
sigma.comp.array=array(data=NA,dim=c(nsamples-1000,4,x.prob))
for (i in 1:x.prob)
    {
        Ymat=discard(n,pmiss[i])
        Y=Ymat@Y
        Y<-Y[complete.cases(Y)==TRUE,]
        source("gibbsrun.r")
        theta.comp.array[,,i]=THETA2
        sigma.comp.array[,,i]=SIGMA2
        theta.mcmc=theta.comp.array[,,i]
        sigma.mcmc=sigma.comp.array[,,i]
    }
theta.mcmc=mcmc(theta.comp.array[,,1])
sigma.mcmc=mcmc(sigma.comp.array[,,1])
summary(theta.mcmc)
summary(sigma.mcmc)
theta.mcmc=mcmc(theta.comp.array[,,2])
sigma.mcmc=mcmc(sigma.comp.array[,,2])
summary(theta.mcmc)
summary(sigma.mcmc)
theta.mcmc=mcmc(theta.comp.array[,,3])
sigma.mcmc=mcmc(sigma.comp.array[,,3])
summary(theta.mcmc)
summary(sigma.mcmc)


####missing analysis
theta.miss.array=array(data=NA,dim=c(nsamples-1000,2,x.prob))
sigma.miss.array=array(data=NA,dim=c(nsamples-1000,4,x.prob))
for (l in 1:x.prob)
    {
        Ymat=discard(n,pmiss[l])
        Y=Ymat@Y
        index=Ymat@index
        source("gibbsmissing.r")
        theta.miss.array[,,l]=THETA2
        sigma.miss.array[,,l]=SIGMA2
    }


plot(theta.mcmc)
plot(sigma.mcmc)
autocorr.plot(theta.mcmc)
autocorr.plot(sigma.mcmc)

#plot
library(hdrcde)
    mu.missing1=mu.missing2=NULL
    mu.missing1.higher=mu.missing1.lower=mu.missing2.higher=mu.missing2.lower=NULL
    var.missing1=var.missing2=cov.missing=NULL
    var.missing1.higher=var.missing1.lower=var.missing2.higher=var.missing2.lower=NULL
    cov.missing.higher=cov.missing.lower=NULL
for(i in 1:x.prob){    
    mu.missing1[i]=mean(theta.miss.array[,1,i])
    mu.missing1.lower[i]=hdr(theta.miss.array[,1,i])$hdr[2,1]
    mu.missing1.higher[i]=hdr(theta.miss.array[,1,i])$hdr[2,2]
    mu.missing2[i]=mean(theta.miss.array[,2,i])
    mu.missing2.lower[i]=hdr(theta.miss.array[,2,i])$hdr[2,1]
    mu.missing2.higher[i]=hdr(theta.miss.array[,2,i])$hdr[2,2]
    var.missing1[i]=mean(sigma.miss.array[,1,i])
    var.missing1.lower[i]=hdr(sigma.miss.array[,1,i])$hdr[2,1]
    var.missing1.higher[i]=hdr(sigma.miss.array[,1,i])$hdr[2,2]
    cov.missing[i]=mean(sigma.miss.array[,2,i])
    cov.missing.lower[i]=hdr(sigma.miss.array[,2,i])$hdr[2,1]
    cov.missing.higher[i]=hdr(sigma.miss.array[,2,i])$hdr[2,2]
    var.missing2[i]=mean(sigma.miss.array[,4,i])
    var.missing2.lower[i]=hdr(sigma.miss.array[,4,i])$hdr[2,1]
    var.missing2.higher[i]=hdr(sigma.miss.array[,4,i])$hdr[2,2]
}

plot(probs,parestimates,type="n")
lines(mu.missing1)
lines(mu.missing2,col="red")
lines(mu.missing1.higher,lty=2)
lines(mu.missing1.lower,lty=2)
lines(mu.missing2.lower,col="red",lty=2)
lines(mu.missing2.higher,col="red",lty=2)

parestimates=range(var.missing1.higher,var.missing1.lower,var.missing2.higher,var.missing2.lower,cov.missing.higher,cov.missing.lower)
plot(probs,parestimates,type="n")
lines(var.missing1)
lines(var.missing2,col="red")
lines(var.missing1.higher,lty=2)
lines(var.missing1.lower,lty=2)
lines(var.missing2.lower,col="red",lty=2)
lines(var.missing2.higher,col="red",lty=2)
lines(cov.missing.lower,col="blue",lty=2)
lines(cov.missing.higher,col="blue",lty=2)
lines(cov.missing,col="blue")

library(hdrcde)
    mu.comp1=mu.comp2=NULL
    mu.comp1.higher=mu.comp1.lower=mu.comp2.higher=mu.comp2.lower=NULL
    var.comp1=var.comp2=cov.comp=NULL
    var.comp1.higher=var.comp1.lower=var.comp2.higher=var.comp2.lower=NULL
    cov.comp.higher=cov.comp.lower=NULL
for(i in 1:x.prob){    
    mu.comp1[i]=mean(theta.miss.array[,1,i])
    mu.comp1.lower[i]=hdr(theta.miss.array[,1,i])$hdr[2,1]
    mu.comp1.higher[i]=hdr(theta.miss.array[,1,i])$hdr[2,2]
    mu.comp2[i]=mean(theta.miss.array[,2,i])
    mu.comp2.lower[i]=hdr(theta.miss.array[,2,i])$hdr[2,1]
    mu.comp2.higher[i]=hdr(theta.miss.array[,2,i])$hdr[2,2]
    var.comp1[i]=mean(sigma.miss.array[,1,i])
    var.comp1.lower[i]=hdr(sigma.miss.array[,1,i])$hdr[2,1]
    var.comp1.higher[i]=hdr(sigma.miss.array[,1,i])$hdr[2,2]
    cov.comp[i]=mean(sigma.miss.array[,2,i])
    cov.comp.lower[i]=hdr(sigma.miss.array[,2,i])$hdr[2,1]
    cov.comp.higher[i]=hdr(sigma.miss.array[,2,i])$hdr[2,2]
    var.comp2[i]=mean(sigma.miss.array[,4,i])
    var.comp2.lower[i]=hdr(sigma.miss.array[,4,i])$hdr[2,1]
    var.comp2.higher[i]=hdr(sigma.miss.array[,4,i])$hdr[2,2]
}

parestimates=range(mu.comp1.higher,mu.comp1.lower,mu.comp2.higher,mu.comp2.lower)

plot(probs,parestimates,type="n")
lines(mu.comp1)
lines(mu.comp2,col="red")
lines(mu.comp1.higher,lty=2)
lines(mu.comp1.lower,lty=2)
lines(mu.comp2.lower,col="red",lty=2)
lines(mu.comp2.higher,col="red",lty=2)

plot(probs,parestimates,type="n")
lines(var.comp1)
lines(var.comp2,col="red")
lines(var.comp1.higher,lty=2)
lines(var.comp1.lower,lty=2)
lines(var.comp2.lower,col="red",lty=2)
lines(var.comp2.higher,col="red",lty=2)
lines(cov.comp.lower,col="blue",lty=2)
lines(cov.comp.higher,col="blue",lty=2)
lines(cov.comp,col="blue")
