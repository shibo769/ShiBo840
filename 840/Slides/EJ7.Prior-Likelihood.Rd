#
# Prior versus likelihood
# Different Priors may look flat in the "relevant parameter range"
# What matters is how peaked the prior is relative to the likelihood

trumu  <- 0.1
trusig <- 0.2	# Sigma known
TT  <- 60

# Pick values of mu for plotting
muvals<-seq(muhat-5*trusig/sqrt(TT),muhat+5*trusig/sqrt(TT),length=201)


# Simulate Sample and plot likelihood vs mu values
# Likelihood computed at MLE forfor sigma

ydat  <- rnorm(TT,trumu,trusig)
muhat <- mean(ydat);muhat
sighat<- sd(ydat)*sqrt(1-1/TT);sighat

# For this sample ydat, we can compute the likelihood
# for various value of the parameter mu
# Pick values of mu for plotting
muvals<-seq(muhat-4*trusig/sqrt(TT),muhat+4*trusig/sqrt(TT),length=201)

# Compute the likelihood for muvals and for sigma=sigma_MLE
lik <- rep(0,201)
for (i in 1:201){
lik[i] <- exp(-sum((ydat-muvals[i])^2)/(2*sighat^2))
}
lik<-lik/max(lik) # Plotting trick: Only the shape counts

# Plot it
par(mgp=c(1.5,0.5,0),mar=c(2.5,2.5,2,0.2))
plot(muvals,lik,xlab="mu",ylab="Likelihood",
type="l",main="Likelihood p(y|mu) vs mu",cex.axis=0.9)
abline(v=c(mu,muhat),col=c("blue","red"),lwd=2)

#
# Look at two different priors for the mean only.
# For example, two normal priors with different mean & std.dev

parvals<-seq(-1,1.5,length=1000)	# pick parameter values
plot(parvals,dnorm(parvals,0.3,0.5),type="l",col="red",
ylim=c(0,1.35),xlab="mean",ylab="density")
lines(parvals,dnorm(parvals,0.05,0.3),col="blue")
abline(v=c(0.3,0.05),col=c("red","blue"))
mtext("Two different priors on the mean")

#
# Plot Priors and likelihood together
#

par(mgp=c(1.5,0.5,0),mar=c(2.5,2.5,2,0.2))
plot(muvals,lik,type="l",ylab="Likelihood and Prior",xlab="Parameter Mu")
axis(side=1,at=muhat,labels="MLE",mgp=c(-2,-1,-0))
text(muvals[130],0.6,"l(m|y)")
title("Likelihood and 2 different Priors",line=0.1)

pteta1<-dnorm(muvals,0.3,0.5)
pteta2<-dnorm(muvals,0.05,0.4)

lines(muvals,pteta1,col="red")
text(muvals[180],0.8,"Prior(0.3,0.5)",col="red")

lines(muvals,pteta2,col="blue")
text(muvals[20],0.95,"Prior(0.05,0.4)",col="blue")


#********************###########################
# Compare lhood(sigma)to improper prior 1/sigma
# Likelihood drawn at the MLE for mu
#

ydat  <-rnorm(TT,trumu,trusig)		# Get a sample
muhat <-mean(ydat)
sighat<-sd(ydat)*sqrt((TT-1)/TT)
sigsig<-sighat/sqrt(2*TT)

# Range of sigma values for likelihood computation
sigvals<-seq(sighat-4*sigsig,sighat+4*sigsig,length=101)

lik <- rep(0,101)
for (i in 1:101){
lik[i]<- sigvals[i]^(-TT)*exp(-sum((ydat-muhat)^2)/(2*sigvals[i]^2))
}
lik<-lik/max(lik) # Plotting trick: Only the shape counts

par(mgp=c(1.5,0.5,0),mar=c(2.5,2.5,2,0.2))
plot(sigvals,lik,type="l",ylab="Likelihood and Prior",xlab="Sigma")
axis(side=1,at=sighat,labels="MLE",mgp=c(-2,-1,-0))
text(teta[75],0.6,"l(m|y)")

# Add the diffuse prior 1/sigma 
# ... rescaled, only shape matters

lines(sigvals, min(sigvals)/sigvals,col="red") 



