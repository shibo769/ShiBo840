#
# Sorting on Estimates
#
# Simulate np portfolio returns over nyear years
# Sort on basis of winners to make long-short hedge portfolio
# Can we predict winners in the next period?

#########################################################
#	Estimate mean vector every years with 5 year RW
#	Make a portfolio of top and bottom 10%

# Sorting on Estimates
#
# Simulate np portfolio returns over nyear years

freq  <- 12				# data frequency
year  <- 5				# length of an estimation period
lenn  <- 8				# Number of estimation periods	
nyear <- lenn*year		# Total number of years of backtesting
np    <- 40				# number of portfolios				

# Actual means are different

mulow <- 0.06 / freq  
muhigh<- 0.18 / freq
muvec <-seq(mulow,muhigh,length=np)		# vector of np TRUE means

trutop<-mean(muvec[(9*np/10):np])*freq
trubot<-mean(muvec[1:(1+np/10)])*freq
trutop;trubot

# Covariance matrix 
rho   <- 0.3
sig   <- 0.25 / sqrt(freq)
cormat<- matrix(rho,ncol=np,nrow=np)+ diag(rep(1-rho,np))
covmat<- diag(rep(sig,np)) %*%cormat %*% diag(rep(sig,np))

# Ready? Simulate the data

library(MASS)	# Mass has mvrnorm for multivariate normal
mrets <- mvrnorm(nyear*freq,muvec,covmat)

library(zoo)
muhats  <-rollapply(mrets,width=year*freq,mean,by=year*freq)*freq
win 	<- matrix(0,nrow=lenn-1,ncol=2)	# will contain winners
loss	<- win							# 			.... losers
for(ii in 1:(lenn-1)){
	sorts 	 <- muhats[c(ii,(ii+1)),order(muhats[ii,])]	# Sort on basis of per. 1
														# keep track of per. 2
	 win[ii,]<- apply(sorts[,(9*np/10):np],1,mean)
	loss[ii,]<- apply(sorts[,1:(1+np/10)],1,mean)
	}
plot(c(win[,1],loss[,1]),c(win[,2],loss[,2]), xlab="Period 1 mean",
ylab="Period 2 mean",xlim=range(muhats),ylim=range(muhats))
title("Predictive power of winner vs loser portolio formed on 5-year estimates",line=0.1,cex.main=1)
abline(lsfit(c(win[,1],loss[,1]),c(win[,2],loss[,2])))
abline(v=c(trubot,trutop),lty=2)
abline(h=c(trubot,trutop),lty=2)


plot(c(win[,1]-loss[,1]),c(win[,2]-loss[,2]), xlab="predicted mean difference",
ylab="realized mean difference",xlim=range(muhats),ylim=range(muhats))
title("Predictive power of W-L portolio formed on 5-year estimates",line=0.1,cex.main=1)

