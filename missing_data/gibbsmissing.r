n=dim(Y)[1]
p=dim(Y)[2]
nu0=p+2
S0=(nu0-p-1)*Sigma0
nmissing=length(index)

Y.full<-Y
O<-1*(!is.na(Y))

for(j in 1:p)
       {
          Y.full[is.na(Y.full[,j]),j]<-mean(Y.full[,j],na.rm=TRUE)
       }

Sigma=cov(Y.full)

set.seed(12)
THETA=SIGMA=Y.MISS=NULL
for (s in 1:nsamples) {
    ##update theta
    ybar<-apply(Y.full,2,mean)
    Ln=solve(solve(L0)+n*solve(Sigma))
    mun=Ln%*%(solve(L0)%*%mu0+n*solve(Sigma)%*%ybar)
    theta=rmvnorm(1,mun,Ln)

    ##update Sigma
    Sn=S0+(t(Y.full)-c(theta))%*%t(t(Y.full)-c(theta))
    Sigma=solve(rwish(1,nu0+n,solve(Sn)))

    ##update missing data
    for(k in 1:nmissing)
        {
            col=ifelse(index[k]>n,2,1)
            row=ifelse(index[k]>n,index[k]-n,index[k])
            b<- (O[row,]==0)
            a<- (O[row,]==1)
	    if(sum(b)==1){
            iSa<- solve(Sigma[a,a])
            beta.j<- Sigma[b,a]%*%iSa
            Sigma.j <- Sigma[b,b]-Sigma[b,a]%*%iSa%*%Sigma[a,b]
            theta.j <- theta[b] + beta.j%*%((Y.full[row,a])-theta[a])
            Y.full[row,b]<-rmvnorm(1,theta.j,Sigma.j)
	    }
        }
    
    #save results
    THETA=rbind(THETA,theta)
    SIGMA=rbind(SIGMA,c(Sigma))
    Y.MISS=rbind(Y.MISS,Y.full[O==0])
}
THETA2=THETA[1001:nsamples,1:2]
SIGMA2=SIGMA[1001:nsamples,1:4]
