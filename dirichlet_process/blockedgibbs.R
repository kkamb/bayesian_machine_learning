# from David Dunson's code

# Fit Simple DP Model to 3-component normal mixture
# Sample Data from y=.1xN(-3,.5) + .5xN(0,.75)+.4xN(3,1) 
# Compares with DPpackage
###################
library(mixtools)
library(Hmisc)
library(DPpackage)

#################
# Generate Data #
#################
set.seed(102923)
n<-500
pi<-c(.1,.5,.4)
mu<-c(-3,0,3)
s2<-c(.5,.75,1)
y<-rnormmix(n,pi,mu,sqrt(s2))
H<-20    # Max Number of clusters in block DP

plot(density(y),type="l",col="darkgreen",lwd=2)

hist(y,breaks=20,freq=F)
lines(density(y,bw=.75),type="l",col="darkgreen",lwd=2)

##########
# Priors #
##########
alpha<-1  # DP Scale Paramteter  (Could try gamma hyperprior)
m0<-0		    # mean of mu all clusters
t0<-.001     # Prec of mu all clusters	
a<-b<-.01	  # Hyperparms for sigma all clusters

#########
# Inits #
#########
pi<-ns<-rep(0,H) 	# Mixing weights and number of subjects per cluster
v<-rep(1/H,H)		# Conditional weights -- pr(c_i=h|c_i not in l<h)
v[H]<-1			# Apply DP truncation to H classes
mu<-rep(0,H)		# Cluster-specific means
tau<-sigma2<-rep(1,H)	
p<-tmp2<-matrix(0,n,H) # p[i,h] = posterior prob that subject i belongs to cluster h

#########
# Store #
#########
nsim<-1000
V<-Mu<-Sigma<-N<-Pi<-matrix(0,nsim,H)
C<-matrix(0,nsim,n)
grid<-seq(min(y),max(y),length=500)
Y<-array(0,dim=c(nsim,length(grid),H))

#########
# GIBBS #
#########
for (i in 1:nsim) {
  # Update c, cluster indicator
  cumv<-cumprod(1-v)
  pi[1]<-v[1]
  for (h in 2:H) pi[h]<-v[h]*cumv[h-1]
  for (h in 1:H) tmp2[,h]<-pi[h]*dnorm(y,mu[h],sqrt(sigma2[h]))
  p<-tmp2/apply(tmp2,1,sum)
  C[i,]<-c<-rMultinom(p,1)
  Pi[i,]<-pi
  for (h in 1:H) ns[h]<-length(c[c==h])  # Must allow zeros for empty clusters
  
  # Update v
  for (h in 1:(H-1)) v[h]<-rbeta(1,1+ns[h],alpha+sum(ns[(h+1):H]))
  V[i,]<-v
  
  # Update mu and sigma2 and Yhat (density estimate)
  for (h in 1:H) {
    var<-1/(t0+tau[h]*ns[h])
    m<-var*(t0*m0+tau[h]*sum(y[c==h])) 
    Mu[i,h]<-mu[h]<-rnorm(1,m,sqrt(var))
    tau[h]<-rgamma(1,a+ns[h]/2,b+sum((y[c==h]-mu[h])^2)/2)
    Sigma[i,h]<-sigma2[h]<-1/tau[h]
    Y[i,,h]<-pi[h]*dnorm(grid,mu[h],sqrt(sigma2[h]))
  }
  N[i,]<-ns  		# Number per cluster
  if (i%%100==0) print(i)
  
}

nh<-table(round(apply(C[501:nsim,],2,mean)))			# Avg Number per class (crude est.)	
mprop<-round(apply(N[501:nsim,],2,mean)[which(nh>0)]/n,2) 	# Avg Proportion per class (another way)
mmu<-round(apply(Mu[501:nsim,],2,mean)[which(nh>0)],2)
msigma<-round(apply(Sigma[501:nsim,],2,mean)[which(nh>0)],2)
cat(mprop,"\n",mmu,"\n",msigma)

Ytmp<-matrix(0,nsim,length(grid))
for (i in 1:nsim) {
  for (j in 1:length(grid)) Ytmp[i,j]<-sum(Y[i,j,])
}

yhat<-apply(Ytmp,2,mean)

###################
# Using DPpackage #
###################
nburn <- 500
nsave <- 500
nskip <- 1
ndisplay <- 100
mcmc <- list(nburn=nburn,nsave=nsave,nskip=nskip,ndisplay=ndisplay)

prior1 <- list(alpha=1,m1=rep(0,1),psiinv1=diag(0.5,1),nu1=4,
               tau1=1,tau2=100)
fit <- DPdensity(y=y,prior=prior1,mcmc=mcmc,
                 state=state,status=TRUE)

################
# Plot Density #
################
#plot(fit,ask=F)
hist(y,breaks=20,freq=F,main="Dirichlet Process Fit to Four-Component Mixture",ylim=c(0,.25))
lines(grid,yhat,type="l",lwd=2,col="darkgreen")
lines(fit$x1,fit$dens,type="l",col="darkred",lwd=2,lty=2)
#lines(density(y),type="l",col="blue4",lwd=2,lty=2)
lines(density(rnormmix(1000000,c(.1,.5,.4),c(-3,0,3),sqrt(s2))),  
      type="l",col="blue4", lty=3,lwd=2)
legend("topright",col=c("darkgreen","darkred","blue4"),lty=c(1,2,3),
       legend=c("Estimated","DPpackage","True Density"),lwd=2,cex=.9)
