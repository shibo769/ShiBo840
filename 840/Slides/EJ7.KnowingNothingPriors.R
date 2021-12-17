# When designing priors ... 
# We can not not nothing about everything !!
# Isn't that nice ?

#
# Rt = rho Rt-1  + noise		Autoregression
#

# 1) Flat prior on rho - Effect on R-square

rho<-runif(10000,-1,1)
Rsquare<-rho^2
par(mfrow=c(2,1),mgp=c(1.5,0.5,0),mar=c(2.5,2.5,2,0.1))
hist(rho,nclass=50,prob=T);box()
hist(Rsquare,nclass=100,prob=T);box()
 

# Note we don't need to simulate, we can compute p(Rsq)
# by change of variable from rho to Rsq
# Make sure you can do that 
xx <-seq(0.001,1,length=1000)
lines(xx,1/(2*sqrt(xx)),col="blue",lwd=2)

# 2) Flat prior on Rsquare - Effect on rho.

Rsquare<-runif(1000000,0,1)		# Simulate prio R-squares
hist(Rsquare,nclass=100,prob=T);box()

# We can allow rho in [-1,1]

rho<-c(sqrt(Rsquare),-sqrt(Rsquare))	 # This trick is called antihetic
hist(rho,nclass=200,prob=T);box()

# Or maybe we don't and force rho>0
rho<-sqrt(Rsquare)
hist(rho,nclass=100,prob=T);box()

#
# Uniform prior on two regression Beta coefficient 
#

bet1<-runif(100000,0,2)
bet2<-runif(100000,0,2)

# Prior for the difference ?
hist(bet2-bet1,nclass=100,prob=T);box() 


#################
# A flexible prior for bounded parameters (like rho or Rsquare)
# The Beta density
# It allows to model "knowing a little" or "knowing a lot"


theta <-seq(0,1,length=1000)

aa<- 1
bb<- 1
ptheta<-dbeta(theta,aa,bb)
parsim<-rbeta(100000,aa,bb)

hist(parsim,nclass=100,prob=T,main="");box()
lines(theta,ptheta)
title(paste("Beta Distribution: a =",aa," b=",bb),line=0.1)

# By linear transformation the Beta domain [0,1]
# can be changed to any [LB,UB], like [-1,1] for a correlation
# coefficient






