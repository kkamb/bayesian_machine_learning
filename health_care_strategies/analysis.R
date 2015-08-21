library(MASS);library(coda);library(xtable);library(car);library(mvtnorm);library(moments);library(tmvtnorm); library(monomvn); library(lars); library(pls);library(leaps)

cost.fulldata=read.table("costs.txt", header=TRUE)
attach(cost.fulldata)
df=data.frame(COST, RXPM, GS,RI,COPAY,AGE,F,MM)

###Exploratory analysis
pdf("dataplot.pdf")
plot(df)
dev.off()

predictors=df[-1]
mat.pred=as.matrix(predictors)
yf = df$COST
Xf=scale(mat.pred)
n=length(yf)

##vif correlation plot
model1=lm(yf~RXPM+GS+RI+COPAY+AGE+F+MM)
vif(model1)

###BAYESIAN OUTLIERS
cost.lm<-lm(COST~RXPM+GS+RI+COPAY+AGE+F+MM,qr=T,data=cost.fulldata)
##added variable plot
pdf("avplot.pdf")
avPlots(cost.lm)
dev.off()

#residuals-outlier plots
pdf("residuals.pdf")
par(mfrow=c(2,2))
plot(cost.lm)
dev.off()

source("bayes-outliers.R")
k = qnorm(.5 + .5*.95^(1/21))
Bout <- Bayes.outlier.prob(cost.lm,k=k)

plot(Bout$prob.outlier, ylab="Posterior Probability of Outlier", xlab="Case", type="h")

indices = outer(1:21, 1:21, FUN=paste)
cbind(indices[Bout$prob.pair.outlier > .0027^2],
      round(Bout$prob.pair.outlier[Bout$prob.pair.outlier > .0027^2],
            digits=6))

install.packages("BMA")
library(BMA)
help(package=BMA)
help(MC3.REG)
cost.MC3= MC3.REG(COST, mat.pred
  ,num.its=10000, outliers=TRUE, M0.out=rep(FALSE, n), outs.list=1:n)

xtable(summary(cost.MC3))

# other diagnostics
im = influence.measures(stack.lm)
plot(rstudent(stack.lm) ~ hatvalues(stack.lm), ylab="Externally Studentized Residual", xlab="Leverage")
identify(rstudent(stack.lm) ~hatvalues(stack.lm) )
# Prob that observation with largest studentized residual is an outlier
2*(1- pt(max(abs(rstudent(stack.lm))), stack.lm$df - 1))
# Bonferonni 
.05/21
max(abs(rstudent(stack.lm)))
qt(1 - .025/21, 16)

#taking out outliers
outlier.vector=c(10,19)
Xfout=mat.pred[-outlier.vector,]
yfout=yf[-outlier.vector]
nout=length(yfout)

##separating into testing and training cases
n.train=as.integer(nout/1.8)
n.test=nout-n.train
in.test = sample(1:nout, n.test)  # random test cases
in.train = (1:nout)[-in.test]  # remaining training data
y = yfout[in.train]
x = Xfout[in.train,]
p=dim(Xfout)[2]
dx.scale=matrix(c(rep(1,n.test),scale(Xfout[in.test,])),ncol=p+1,nrow=n.test)
dx.unscaled=matrix(c(rep(1,n.test),Xfout[in.test,]),ncol=p+1,nrow=n.test)

pdf("GStest.pdf")
par(mfrow=c(1,2))
plot(Xfout[in.test,2],yfout[in.test],main="Test vs Training Data for GS",col="blue",type="p",pch=16)
abline(lm(yfout[in.test]~Xfout[in.test,2]),col="blue")
lines(yfout,Xfout[,2],col="red",type="p")
       abline(lm(yfout~Xfout[,2]),col="red")
plot(Xfout[in.test,3],yfout[in.test],main="Test vs Training Data for RI",col="blue",type="p",pch=16)
abline(lm(yfout[in.test]~Xfout[in.test,3]),col="blue")
lines(yfout,Xfout[,3],col="red",type="p")
       abline(lm(yfout~Xfout[,3]),col="red")
dev.off()

######modeling using g-priors##########
####################################
library(BMS)
library(hdrcde)#to compute uncertainty
S<-1000
g=length(y)
nu0=1
s20=var(y)
p=dim(X)[2]
X=matrix(c(rep(1,n.train),Xfout[in.train,]),ncol=p+1,nrow=n.train)

lm.gprior<-function(y,X,g=dim(X)[1],nu0,s20,S=1000)
{

  n<-dim(X)[1] ; p<-dim(X)[2]
  Hg<- (g/(g+1)) * X%*%solve(t(X)%*%X)%*%t(X)
  SSRg<- t(y)%*%( diag(1,nrow=n)  - Hg ) %*%y

  s2<-1/rgamma(S, (nu0+n)/2, (nu0*s20+SSRg)/2 )

  Vb<- g*solve(t(X)%*%X)/(g+1)
  Eb<- Vb%*%t(X)%*%y

  E<-matrix(rnorm(S*p,0,sqrt(s2)),S,p)
  beta<-t(  t(E%*%chol(Vb)) +c(Eb))

  list(beta=beta,s2=s2)                                
}   

greg=lm.gprior(y,X,g,nu0,s20,S)
BETA=greg$beta[501:1000,]
SIGMA=greg$s2[501:1000]

#uncertainty using HDR
mean(BETA[,1]);hdr(BETA[,1],95)
mean(BETA[,2]);hdr(BETA[,2],95)
mean(BETA[,3]);hdr(BETA[,3],95)
mean(BETA[,4]);hdr(BETA[,4],95)
mean(BETA[,5]);hdr(BETA[,5],95)
mean(BETA[,6]);hdr(BETA[,6],95)
mean(BETA[,7]);hdr(BETA[,7],95)
mean(BETA[,8]);hdr(BETA[,8],95)
mean(SIGMA);hdr(SIGMA,95)

#uncertainty using quantile
quantile(BETA[,1],c(.025,.5,.975))
quantile(BETA[,2],c(.025,.5,.975))
quantile(BETA[,3],c(.025,.5,.975))
quantile(BETA[,4],c(.025,.5,.975))
quantile(BETA[,5],c(.025,.5,.975))
quantile(BETA[,6],c(.025,.5,.975))
quantile(BETA[,7],c(.025,.5,.975))
quantile(SIGMA,c(.025,.5,.975))

beta.mcmc=mcmc(BETA)

ztrial=zlm(df)
summary(ztrial)

##############modeling using LASSO#################
################################################

scalex=scale(x)
cost.lars = lars(scalex,y, type="lasso")
Cp = summary(db.lars)$Cp
pdf("lassograph.pdf")
plot(cost.lars)
dev.off()
xtable(summary(cost.lars))
round(coef(cost.lars),3)

best = (1:length(Cp))[Cp == min(Cp)]     # step with smallest Cp
LASSO.Est=coef(cost.lars,s=best)
y.pred = predict(cost.lars, s=best, newx=scale(Xfout[in.test,]))
mse.lasso =  sum((yfout[in.test] - y.pred$fit)^2)/n.test

cost2.lm<-lm(y~scalex)
##added variable plot
pdf("avplot2.pdf")
avPlots(cost2.lm)
dev.off()

##############modeling using Ridge Regression############
######################################################

cost.ridge=lm.ridge(y ~ scalex, lambda=seq(0,.1,.0001))
best=which.min(cost.ridge$GCV)
y.pred=dx.scale%*%coef(cost.ridge)[best,]
RIDGE.Est=coef(cost.ridge)[best,]
mse.ridge=sum((yfout[in.test]-y.pred)^2)/n.test
pdf("ridgegraph.pdf")
plot(cost.ridge)
dev.off()

cost.ols=lm(y~x)
confint(cost.ols)
dx2.unscaled=matrix(c(rep(1,n.test),Xfout[in.test,]),ncol=p+1,nrow=n.test)
y.pred=dx2.unscaled%*%coef(cost.ols)
OLS.Est=coef(cost.ols)
mse.ols=sum(yfout[in.test]-y.pred)^2/n.test


############modeling using Bayesian Model Averaging################
###########################################
install.packages("BAS")
library(BAS)
costtrain.bma=bas.lm(y~x, data=cost.fulldata, prior="ZS-null", n.models=2^7, initprobs="Uniform")
#set alpha at n, number of data points. should be 67 for training data.
#set n.models at the number of response variables
#update function not necessary - optional

pdf("bmadiagnostics.pdf")
par(mfrow=c(2,2))
plot(costtrain.bma)
dev.off()

pdf("bmacolor.pdf")
image(costtrain.bma)
dev.off()

pdf("bmacoefs.pdf")
beta=coef(costtrain.bma)
par(mfrow=c(3,3))
plot(beta,subset=2:9,ask=F)
dev.off()

which.mat=list2matrix.which(costtrain.bma,1:(2^8))
    GS.poll.in=(which.mat[,3])==0
    post.prob=1-sum(GS.poll.in*costtrain.bma$postprob)
    odds=(1-post.prob)/(post.prob)
    prior.odds=(1-.2^1)/.2^1
    GS.bayes.factor=odds/prior.odds
    RI.poll.in=(which.mat[,4])==0
    post.prob=1-sum(RI.poll.in*costtrain.bma$postprob)
    odds=(1-post.prob)/(post.prob)
    prior.odds=(1-.2^1)/.2^1
    RI.bayes.factor=odds/prior.odds

############other variable selection procedures#########

##forward-backward selection - frequentist methodology
library(leaps)
LEAP=leaps(Xfout[in.train,],yfout[in.train],int=TRUE,method="adjr2",nbest=1)
LEAP.CP=leaps(Xfout[in.train,],yfout[in.train],int=TRUE,method="Cp",nbest=1)
LEAP
LEAP.CP

##BIC
check=regsubsets(yfout[in.train]~Xfout[in.train,],data=df,nbest=1)
pdf("BICplot-bic.pdf")
plot(check,scale="bic")
dev.off()
pdf("BICplot-adjr2.pdf")
plot(check,scale="adjr2")
dev.off()

##adjusted AIC
null=lm(yfout[in.train]~1)
fullmodel=lm(yfout[in.train]~Xfout[in.train,])
STEP=stepAIC(null,scope=list(lower=null, upper=fullmodel),direction="forward")
summary(STEP)

##########fitting final models#######
BIC.model1.out=lm(yfout[in.test]~Xfout[in.test,1]+Xfout[in.test,2]+Xfout[in.test,5]+Xfout[in.test,7])
BIC.model2.out=lm(yfout[in.test]~Xfout[in.test,1]+Xfout[in.test,2]+Xfout[in.test,5])
BIC.model3.out=lm(yfout[in.test]~Xfout[in.test,2]+Xfout[in.test,5])
AIC.model.out=lm(yfout[in.test]~Xfout[in.test,1]+Xfout[in.test,2])

xtable(summary(BIC.model1.out))
xtable(summary(BIC.model2.out))
xtable(summary(BIC.model3.out))
xtable(summary(AIC.model.out))
xtable(anova(BIC.model3.out,BIC.model2.out,BIC.model1.out))

##testing the potential outliers in our chosen model##
pdf("avplotfinal.pdf")
par(mfrow=c(1,2))
avPlot((BIC.model3.out),Xfout[in.test,2])
avPlot((BIC.model3.out),Xfout[in.test,5])
dev.off()
