## A Dirichlet process model ##
## test to see how well a DP model can estimate the simulated y below
y=.75*rnorm(100,0,1)+.15*rnorm(100,5,10)+.1*rnorm(100,-1,1) #simulate y
n=length(y)
nsim=200; nburn=100; nstore=100 #MCMC run parameters
pars=2 #number of parameters to be estimated
psi=array(0,dim=c(nsim,pars,n)); psi[,2,]=1
mu.star=list(); sigma.star=list()
n.star=rep(0,nsim)
alpha=rep(0,nsim)
munot=rep(0,nsim); taunot=rep(0,nsim); sigmanot=rep(0,nsim)
alpha[1]=1
munot[1]=0; taunot[1]=1; sigmanot[1]=1
m=0; v=3; k=.01; r=.01; c=.01; anot=1; bnot=1; t=1

##########to be deleted after final#############
tableindex=chinese_restaurant(n,alpha[1])
count=0;mut=NULL;sigmat=NULL
for(kcr in 1:n){
  if (count<tableindex[kcr]){
    count=count+1;
    mut[count]=rnorm(1,munot[1],taunot[1])
    sigmat[count]=rgamma(1,v/2,v*sigmanot[1]^2/2)
  }
  psi[1,1,kcr]=mut[tableindex[kcr]]
  psi[1,2,kcr]=sigmat[tableindex[kcr]]
}
psi.star=matrix(union(psi[1,,],NULL),nrow=2)
n.star=length(psi.star[1,])
w=list()
for(i in 1:n.star){
  w[[i]]=seq(1,n)[psi[1,1,]==psi.star[1,i]]
}

sno=2
psi=psi[sno-1,,]
alpha=alpha[sno-1]
munot=munot[sno-1]
taunot=taunot[sno-1]
sigmanot=sigmanot[sno-1]

psi.wout.i=psi[,2:n]
y.i=y[1]
#################################################

for(sno in 2:nsim){
  psi[sno,,] = sample.psi(psi[sno-1,,], n, alpha[sno-1], munot[sno-1], taunot[sno-1], v, sigmanot[sno-1], y)
  mu.star[[sno]]=union(psi[sno,1,],NULL)
  sigma.star[[sno]]=union(psi[sno,2,],NULL)
  n.star[sno]=length(mu.star[[sno]]) #number of tables
  w=list() #observations that fall within that particular c
  #i.e., w[[1]] gives list of observations at the 1st table, w[[2]] list of obs at the 2nd table, etc
  for(i in 1:(n.star[sno])){
    w[[i]]=seq(1,n)[psi[sno,1,]==(mu.star[[sno]][i])]
  }
  #####update cluster locations (i.e., m_js and sigma_j's)
  mu.star[[sno]]=sample.mu.star(w,n.star[sno],sigma.star[[sno]],munot[sno-1],taunot[sno-1],y)
  sigma.star[[sno]]=sample.sigma.star(w,n.star[sno],mu.star[[sno]],sigmanot[sno-1],v,y)
  for(i in 1:(n.star[sno])){
    psi[sno,1,w[[i]]]=(mu.star[[sno]])[i]
    psi[sno,2,w[[i]]]=(sigma.star[[sno]])[i]    
  }
  #####update hyperparameters
  alpha[sno]=sample.alpha(alpha[sno-1],n,n.star[sno],anot,bnot)
  munot[sno]=sample.mu(mu.star[[sno]],n.star[sno],taunot[t-1],m,k)
  taunot[sno]=sample.tau(mu.star[[sno]],n.star[sno],taunot[sno],munot[sno],r,k,t)
  sigmanot[sno]=sample.sigma(sigma.star[[sno]],n.star[sno],v,c)
}

#########FUNCTIONS####################

sample.psi=function(psi,n,alpha,munot,taunot,v,sigmanot,y){
  psi[,1]=sample.psi.i(psi[,2:n],alpha,munot,taunot,v,sigmanot,y[1])
  for(i in 2:(n-1)){
    psi[,i]=sample.psi.i(psi[,c(1:(i-1),(i+1):n)],alpha,munot,taunot,v,sigmanot,y[i])
  }
  psi[,n]=sample.psi.i(psi[,1:(n-1)],alpha,munot,taunot,v,sigmanot,y[n])
  return(psi)
}

sample.psi.i=function(psi.wout.i,alpha,munot,taunot,v,sigmanot,y.i){
  psi.star.wout.i=matrix(union(psi.wout.i,NULL),nrow=2) #unique table indices, with mean in 1st row and sigma in 2nd
  n.star.wout.i=length(psi.star.wout.i[1,]) #length of unique table index
  #number of people sitting at each table is n.j.minus
  n.j.minus=apply(matrix(psi.star.wout.i[1,],nrow=n.star.wout.i,ncol=1),1,function(a)sum(psi.wout.i[1,]==a))
  #auxilary sampling method
  q.0=alpha
  q.j=n.j.minus
  ddd=q.0+sum(q.j)
  uu=runif(1,0,1)
  if(uu<(q.0/ddd)){
    mu=rnorm(1,munot[sno-1],taunot[sno-1])
    sigma=1/rgamma(1,v/2,v*sigmanot[sno-1]/2)
    mean.par=(munot/taunot + y.i*sigma^2)/(1/taunot+sigma^2)
    sigma.par=1/(1/taunot + sigma^2)
    mean.h=rnorm(1,mean=mean.par,sd=sigma.par)
    v.par = v+1
    s.par = v.par^(-1) * (v*sigma^2 + (y.i-mean.h)^2)
    sigma.h=rgamma(1,v.par/2,v.par*s.par/2)
    return(matrix(c(mean.h,sigma.h),nrow=2))
  }  
  else{
    return(matrix(c(sample(psi.star.wout.i[1,],1,replace=TRUE,prob=q.j),sample(psi.star.wout.i[2,],1,replace=TRUE,prob=q.j)),nrow=2))
  }
}

sample.sigmanot=function(y,n,v,c){
  a=c-n*v/2
  b=c-n*v*y/2
  return(rgamma(1,a,b))
}

sample.taunot=function(mu.star,n.star,taunot,munot,r,k,t){
  s2=var(mu.star)
  vn=r+n.star
  s2n=(r*t^2+(n.star-1)*s2+k*n.star*(mean(mu.star)-munot)^2/(k+n.star))/vn
  return(1/rgamma(1,vn/2,vn*s2n/2))
}
sample.munot=function(mu.star,n.star,taunot,m,k){
  munot.mean=(k*m+sum(mu.star))/(k+n.star)
  return(rnorm(1,munot.mean,taunot/(k+nstar))
}

sample.alpha=function(alpha.old,n,n.star,anot,bnot){
  eta=rbeta(1,alpha.old+1,n)
  tt1=anot+n.star-1
  tt2=bnot-log(eta)
  prob=tt1/((n*tt2)+anot+n.star-1)
  uu=runif(1,0,1)
  if (uu<prob){
    return(rgamma(1,tt1+1,tt2))
  }
  else{
    return(rgamma(1,tt1,tt2))
  }
}

sample.mu.star=function(w,n.star,sigmaold,munot,taunot,y){
  mu.star=rep(0,nstar)
  for(i in 1:nstar){
    y.star=y[w[[i]]]
    sigmaold.star=sigmaold[w[[i]]]
    n.j=length(y.star)
    mean.j=(munot*taunot^2+n.j*mean(y.star)/(sigmaold.star^2))/(taunot^2+n.j/sigmaold.star^2)
    var.j=(taunot^2+n/sigmaold.star^2)^(-1)
    mu.star[i]=rnorm(1,mean=mean.j,sd=sqrt(var.j))
  }
  return(mu.star)
}
sample.sigma.star=function(w,n.star,muold,sigmanot,v,y){
  sigma.star=rep(1,nstar)
  for(i in 1:nstar){
    y.star=y[w[[i]]]
    muold.star=muold[w[[i]]]
    n.j=length(y.star)
    vn=v+n.j
    sigman=(1/vn)*(v*sigmanot^2+sum((y.star-muold.star)^2))
    sigma.star[i]=1/rgamma(1,vn/2,vn*sigman/2)
  }
  return(sigma.star)
}