###############################################################
4.1 Basic Simulation

set.seed(1)
runif(5)
runif(5)
set.seed(1)
runif(5)


4.2 Inverse transform 

# Exponential by inverse transform
unifdraw <-   runif(100000)

# ANTIHETHIC ACCELERATION
c(1-unifdraw,unifdraw)

x1       <- - log(unifdraw)		# Exponential by Inverse transform
x2		 <-  rexp(100000)
par(mfrow=c(1,3))
hist(x1,nclass=50,prob=T);box()
hist(x2,nclass=50,prob=T);box()
qqplot(x1,x2);abline(0,1)

#  P.13 NORMAL BY BOX-MULLER (from uniforms)

u1 <-runif(10000)
u2 <-runif(10000)
x1 <- sqrt(-2*log(u1))*cos(2*pi*u2)
x2 <- sqrt(-2*log(u1))*sin(2*pi*u2)
par(mfrow=c(1,3))
plot(x1,x2)
qqnorm(x1)
qqnorm(x2)


# P 14  CORRELATED NORMAL BY CHOL (from uncorrelated)
#	Chol(V) gives you directly the P* in V=P*P*'

x1<-rnorm(20000)
x2<-0.5*x1 + sqrt(1-0.5^2)*rnorm(20000)
cor(x1,x2);sd(x1);sd(x2)

covmat<-matrix(c(1,0.5,0.4,0.5,1,0.3,0.4,0.3,1),ncol=3)
t(chol(covmat))%*%chol(covmat)
simx <-matrix(rnorm(3*200000),ncol=3)
round(var(simx),2)
simx <-simx%*%chol(covmat)
round(var(simx),2)


# P8  Student-T by Data Augmentation
#		aka "mixture representation"

par(mfrow=c(1,2))
dof      <- 8
lambdavec<- dof / rchisq(100000,dof) # an inverse gamma with s = 1
epsvec   <- rnorm(100000)
uvec     <- epsvec * sqrt(lambdavec)
hist(uvec,nclass=80,prob=T);box()
round(mean(uvec),2);round(var(uvec)/(dof/(dof-2)),3)
qqnorm(sort(uvec),type="l");abline(0,1)

# 5.3 ACCEPT REJECT P 16  ###########################
# Need package VGAM for laplace
# Draw Normal by A/R on Laplace blanket

library(VGAM)
par(mfrow=c(1,1))
xx<-seq(-4,4,by=0.001)plot(xx,dlaplace(xx),type="l")lines(xx,dnorm(xx),lty=2)# NEED TO FIND c:  c q(x) > p(x)
# Need to find the maximum q(x)/p(x). It will be our c# Lazy way: Just plot q/p to find the maximum

plot(xx,dnorm(xx)/dlaplace(xx),type="l",ylab="p/q")
max(dnorm(xx)/dlaplace(xx))

# Analytically, find the maximum of q/p

cc<-sqrt(2/pi)*exp(0.5);cc

# Does it work?plot(xx,dlaplace(xx)*cc,type="l",ylab="c * q  and p")lines(xx,dnorm(xx),lty=2,col="red")
legend("topleft",c("c*q","p"),lty=c(1,2),col=c("black","red"))
# Simulate from Q, from Uniform, keep “good” draws

nsim  <- 100000qdraw <- rlaplace(nsim)udraw <- runif(nsim)pdraw  <- qdraw[udraw < dnorm(qdraw)/( cc* 0.5 *exp(-abs(qdraw)))]

# Check Results
# rejection rate?

1-length(pdraw)/nsim

xx<-seq(-3,3,by=0.001)
par(mfrow=c(1,2))hist(qdraw,nclass=150,freq=F);box();lines(xx,dnorm(xx),col="red") hist(pdraw,nclass=100,freq=F);box();lines(xx,dnorm(xx),col="red") 

#	A/R Normal by Laplace, 
#   Not knowing the normalization constant of p.
#	We only use the kernel p*

xx<-seq(-3,3,by=0.001)
plot(xx,exp(-xx^2/2)/dlaplace(xx),type="l",ylab="p*/q")

max(exp(-xx^2/2)/dlaplace(xx))		# MAX ( p* / q)

#Analytically, the max of p*/q is

cc<-2*exp(0.5);cc

nsim  <- 50000
qdraw <- rlaplace(nsim)
udraw <- runif(nsim)
pdraw <- qdraw[udraw < exp(-0.5*qdraw^2)/(dlaplace(qdraw)*cc)]
1-length(pdraw)/nsim
hist(qdraw,nclass=150,freq=F,xlab="qdraws",col="green");box()
lines(xx,dnorm(xx),col="red",lwd=2) 
hist(pdraw,nclass=150,freq=F,xlab="pdraws",col="green");box();
lines(xx,dnorm(xx),col="red") 

par(mfrow=c(1,2))
qqnorm(qdraw,ylab="qdraw",pch=20);abline(0,1)
qqnorm(pdraw,ylab="pdraw",pch=20);abline(0,1)


###############################################################
# Draw a Beta by q ~ uniform with 
# INDEPENDENCE METROPOLIS
# We could do it by A/R also

x<-seq(0,1,by=0.01)
plot(x,dbeta(x,3.1,7.5),type="l")		#3.1  7.5
lines(x,dunif(x),lty=2)

nsim	<-100000
xx      <-rep(0,nsim) 	# accepted draws will be here
xcandid <-runif(nsim) 	# candidate draws from q
xx[1]	<-0.99			# Arbitrary start in the domain
uvec    <-runif(nsim) 	# draw for accept/repeat step
accept  <-rep(0,nsim)	# and acceptance decision at each step


# You should be able to tell why q disappeared from the acceptance
# probability

for (i in 2:nsim){
pcan  <- xcandid[i]^(3.1-1)*(1-xcandid[i])^(7.5-1)	# p*(candidate)
pmin1 <- xx[i-1]^(3.1-1)*(1-xx[i-1])^(7.5-1)   		# p*(previous draw)
prob  <- pcan / pmin1			# acceptance probability
accept[i]   <- (uvec[i]<prob)	# acceptance decision
xx[i]  		<- xx[i-1] * (1-accept[i]) + xcandid[i]* accept[i]
}

round((nsim-length(unique(xx)))/nsim,2)		# % repeats


library(forecast)	# to use adult Acf

par(mfrow=c(1,2),mgp=c(1.5,0.5,0))
ts.plot(xx[1:1000],ylab="Draw",main="")
title("Metropolis: time series of accepted draws",line=0.2)
Acf(xx,lag.max=20,main="")
title("Metropolis: ACF of accepted draws",line=0.2)

par(mfrow=c(1,1))
hist(xx,nclass=150,freq=F);box()
lines(x,dbeta(x,3.1,7.5),col="blue",lwd=2)

##############################################################
# GIBBS SAMPLER:
# Bivariate Normal by the two conditionals x1|x2  and x2|x1

rho  <- 0.99; nsim <-100000; 
mu1  <- 0.5; mu2  <- 2; 
condsd <-sqrt(1-rho^2)		# Conditional standard deviation
t1draw<-rnorm(nsim); 
t2draw<-rnorm(nsim) 	# Draw N(0,1)s ahead of time

t1draw[1]<- 7		# 1st draw, Any acceptable value in domain
t2draw[1]<- t2draw[1]*condsd + rho*t1draw[1]  + mu1-rho*mu2

for (i in 2:nsim){
t1draw[i]<- t1draw[i]*condsd + rho*t2draw[i-1] + mu1-rho*mu2
t2draw[i]<- t2draw[i]*condsd + rho*t1draw[i]   + mu2-rho*mu1 }

cor(t1draw,t2draw)
tab<-matrix(c(mean(t1draw),sd(t1draw), mean(t2draw),sd(t2draw)),ncol=2)
colnames(tab)<-c("X1","X2")
rownames(tab)<-c("Mean","Std.dev")
round(tab,3)

par(mfrow=c(2,1))
ts.plot(t1draw[1:1000]);title("X1 | X2 Draws",line=-1) 
ts.plot(t2draw[1:1000]);title("X2 | X1 Draws",line=-1) 

Acf(t1draw[101:100000],lag.max=100,main="X1|X2");box()
Acf(t2draw[101:100000],lag.max=100,main="X2|X1");box()

# Numerical Efficiency 
# Can use package bayesm
library(bayesm)
numEff(t1draw)
numEff(t2draw)

xx<-seq(-3.5,4.5,length=500)hist(t1draw,nclass=150,col="green",prob=T);box();abline(v=mu1,lwd=3) 
lines(xx,dnorm(xx,mu1,1),lwd=2,col="red")
hist(t2draw,nclass=100,col="green",prob=T);box();abline(v=mu2,lwd=3)lines(xx,dnorm(xx,mu2,1),lwd=2,col="red")


################################################################
# Nasty Bimodal Target 
# by INDEPENDENCE Metropolis
# 
par(mfrow=c(1,1))
pkern<-function(x){
pkern<- exp(-(x-1)^2/(2*0.45^2))*(2.8*exp(-2.8*(-x+6)) + 
		exp(-(x+2)^2/(2*0.35^2)))
pkern
}
# We don't know the normalization constant of p, we only have p*
# Plot standardizes everything by its max so to show both on same plot

x<-seq(-2,3.5,by=0.01)
plot(x,pkern(x)/max(pkern(x)),type="l",ylab="p* unnormalized")
lines(x,dnorm(x,1.5,1.5)/max(dnorm(x,1.5,1.5)),lty=2) # Use q ~ Normal

# Only know the density kernel: 

nsim	<-50000
xx     	<-rep(0,nsim) 			# accepted draws will be here
xcandid <-rnorm(nsim,1.5,1.5) 	# draws from NORMAL q proposal
uvec    <-runif(nsim) 			# draw for accept/repeat step
prob    <- rep(0,nsim)
xx[1]   <- xcandid[1]

for (i in 2:nsim){ 
pcan  	<- pkern(xcandid[i]) 	# candidate draw kernel
plast 	<- pkern(xx[i-1]) 		# Previous draw kernel
prob[i] <- min((pcan / dnorm(xcandid[i],1.5,1.5)) / (plast / dnorm(xx[i-1],1.5,1.5) ),1)
accept <- ( uvec[i] < prob[i]) 
xx[i]  	<- xx[i-1] * (1-accept) + xcandid[i]* accept
}

1-length(unique(xx))/nsim			# Fraction of repeats

par(mfrow=c(3,1),mar=c(3.5,3,2,0.1),mgp=c(1.5,0.5,0))Acf(xx[1:10000],lag.max=200) ts.plot(xx[1:10000],ylab="draw")hist(xx,nclass=100,freq=F);box()lines(x,pkern(x)*  max(hist(xx,plot=F)$density)/max(pkern(x)))

#
# Now a random walk METROPOLIS
# See effect of step size on results

nsim	<- 500000
xx      <- rep(0,nsim) 	# accepted draws will be here
sigma   <- 1 				# Random Walk Std. Dev.
step    <- rnorm(nsim,0,sigma)	#	Random Walk Step Size 
uvec    <- runif(nsim) 	# draw for accept/repeat step
prob    <- rep(0,nsim)
accept  <- rep(0,nsim)
xx[1]   <- 3				# Start with anything

for (i in 2:nsim){ 
xcandid	<- xx[i-1] + step[i]		# The RW candidate draw
pcan   	<- pkern(xcandid) 
plast  	<- pkern(xx[i-1]) 
prob[i] <- min(  pcan  / plast  , 1 ) 
accept[i]   <- ( uvec[i] < prob[i]) 
xx[i]     	<- xx[i-1] * (1-accept[i]) + xcandid * accept[i]
}

1-length(unique(xx))/nsim	# Fraction of Repeats

par(mfrow=c(3,1),mar=c(2,3,1,0))Acf(xx,lag.max=100) ts.plot(xx[1:50000])hist(xx,nclass=100,freq=F);box()lines(x,pkern(x)*max(hist(xx,plot=F)$density)/max(pkern(x)))




