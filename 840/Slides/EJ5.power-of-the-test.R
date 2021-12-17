#
#	COMPUTING POWER OF THE TEST, TYPE I and TYPE II ERRORS
#

###############################################################
#  Estimation, property of the sample statistic under the null
#

sig <-  1 
TT  <- 200 
mu0 <-  0 
lb5 <- mu0 - 1.96 *sig/sqrt(TT)	# 5% bounds used to accept/reject under H0
ub5 <- mu0 + 1.96 *sig/sqrt(TT) 
c(lb5,ub5)				# Rejection bounds at 5%

##############################################################
#
# What is the distribution of muhat given the true mu?
# Null H0: 			mu = 0,   
# Alternative H1: 	mu takes other values
# Need to pick one alternative value: trumu 
# Under H1: 			muhat ~ N(trumu, sig/sqrt(T))
# Plot the distribution of muhat under the alternative mu=trumu

par(mgp=c(1.5,0.5,0))
trumu    <- 0.08					# Alternative values, mu1 = 0 is H0
muvals   <-seq(-0.4,0.8,len=400)		# values to draw density
denmuhat <-dnorm(muvals,trumu,sig/sqrt(TT))	# Compute density

plot(muvals,denmuhat,type="l",ylab="Density",xlab="mu")
abline(v=c(lb5,ub5),col="red")
abline(v=trumu,lty=2,lwd=2)
legend("topright",c("95% bounds for H0",expression(paste("True ",mu))),lty=c(1,2),lwd=c(1,2),col=c("red","black"),bty="n")
title(bquote(paste("Sampling Distribution of ",hat(mu),"  under alternative H1: ",mu == .(trumu))), line=1)

# Where is the area of "rejection"
# Shade the areas of rejection 
# ... when muhat is less than lb or more than ub

# If muhat > upper bound under H0
xcoord<-c(ub5,ub5,muvals[muvals > ub5],max(muvals))    
ycoord<-c(0,dnorm(ub5),dnorm(muvals[muvals > ub5],trumu,sig/sqrt(TT)),0)
polygon(xcoord,ycoord,col="red")

# If muhat < lower bound under HO
xcoord<-c(muvals[1],muvals[muvals< lb5],lb5,lb5)
ycoord<-c(0,dnorm(c(muvals[muvals<lb5],lb5),trumu,sig/sqrt(TT)),0)
polygon(xcoord,ycoord,col="red")

# Compute the power !

leftpower <-  pnorm(lb5,trumu,sig/sqrt(TT))	# Prob(muhat < lb | mu=mu1)
rightpower<- 1-pnorm(ub5,trumu,sig/sqrt(TT))	# Prob(muhat > ub | mu=m1)
power<-  leftpower + rightpower

round(rbind(leftpower,rightpower,power),4)

legend("left",legend= bquote(paste("p(",hat(mu),"<LB) = ",.(round(leftpower,3)))),
bty="n",cex=0.8)
legend("right",legend=bquote(paste("p(",hat(mu),">UB) = ",.(round(rightpower,3)))),
bty="n",cex=0.8)

legend("topleft",legend=bquote(paste("Power = ",.(round(power,3)))),bty="n",col="red",cex=0.8)

###########################################################
#
# Now consider the entire alternative space for the "True" mu
# We can compute the power for every possible value of mu
# Plot the POWER (curve) of the test vs mu


sig <-  1 
mu0 <-  0 
mualt<-seq(-0.6*sig,0.6*sig,length=400)	# Possible alternative values of mu

# Compute Bounds for a 5% and a 1% test

TT<-2000
lb5 <- mu0 - 1.96 *sig/sqrt(TT)		# 5% bounds
ub5 <- mu0 + 1.96 *sig/sqrt(TT)
lb1 <- mu0 - 2.576*sig/sqrt(TT)		# 1% bounds 
ub1 <- mu0 + 2.576*sig/sqrt(TT) 

# Power: Probability of rejecting H0 when H0 is false

pow5<-pnorm(lb5,mualt,sig/sqrt(TT))+    1-pnorm(ub5,mualt,sig/sqrt(TT))
pow1<-pnorm(lb1,mualt,sig/sqrt(TT))+    1-pnorm(ub1,mualt,sig/sqrt(TT))

# What is power when mu = mu0 ?
# Prob(reject null when the null is TRUE) = Prob(Type 1 error)

par(mgp=c(1.5,0.5,0))
plot(mualt,pow5,type="l",ylim=c(0,1),xlab=expression(paste("Unknown True ",mu)),
ylab="Prob.(Rejecting H0)")
title(expression(paste("Power of classical test vs True  ",mu)),line=0.5)
lines(seq(-0.1,0.1,len=2),rep(0.05,2),lty=2,lwd=2)

lines(mualt,pow1,lty=2,col="red")			# Power for the 1% test
lines(seq(-0.1,0.1,len=2),rep(0.01,2),lty=2,col="red",lwd=2)
legend("right",legend=c("5% test","1% test"),col=c("black","red"),lwd=1,bty="n")

#
# Pick a large TT, 
#	Compare P(Type I) and P(Type II) 
#	What is your probability of failing to reject the null if wrong
#	What is your prob. of rejecting the null if true: 5%


# type II error: failing to reject the null when the null is false
#	p(type II error) = 1 - power
#   At the null, failing to reject is the right thing to do!

plot(mualt,1-pow5,type="l",ylab="Pr.(Failing to reject H0)",
xlab=expression(paste("Unknown True ",mu)),ylim=c(0,1))
title(line=0.2,"Probability of Type II error: failing to reject Null when Null true",cex.main=1.1)
abline(h=0.05,lty=2,col="red")
axis(4,at=0.05,cex=0.5,col="red")
legend("topright",legend="Prob(Type I)=0.05",lty=2,cool="red",bty="n")




