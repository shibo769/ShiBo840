nn<-50	# n assets
TT<-250	# T periods

# Simulate uncorrelated returns
totvar<-seq(0.18,0.45,length=nn)^2
retmat<-matrix(rnorm(nn*TT,0,rep(sqrt(totvar),rep(TT,nn))),nrow=TT,ncol=nn) # Uncorrelated Data

# Simulate an index model
mum <- 0.08; varm <- 0.15^2
totvar<- seq(0.18,0.45,length=nn)^2
sigeps<-sqrt(totvar-varm)	# residual stdev's
Rsq   <- varm/totvar;Rsq		# all betas=1
rm	  <- rnorm(TT,mum,sqrt(varm))
reps  <- t(matrix(rnorm(nn*TT,0,sigeps),ncol=TT))
retmat<- reps + rm

# Simulate an equal correlation model
library(mnormt)
rho  <-0.5
diagsdmat<-diag(seq(0.18,0.45,length=nn))
cormat<- matrix(rho,ncol=nn,nrow=nn)+ diag(1-rho, ncol=nn,nrow=nn)
covmat<- diagsdmat %*% cormat %*% diagsdmat
retmat<-rmnorm(TT,varcov=covmat)

# correlations
corr<-cor(retmat);corr<-corr[corr<1]
plot(density(corr,adjust=1.5),xlab="Rho",main="Correlations") # Correlations
sd(corr)/sqrt(1/TT)

pcs1   <-prcomp(retmat)		
summary(pcs1)
pcs1$rotation		# The B's
apply(pcs1$rotation^2,2,sum)	# For each component the ei^2 sum to 1
screeplot(pcs1)

pcs2   <-princomp(retmat)	# gets you the scores		
summary(pcs2)
pcs2$loadings		# The B's
pcs2$sdev
pcs2$scores[,50]	# The scores
screeplot(pcs2)

facs<- factanal(retmat,factors=1,scores="Bartlett")	# Joreskog's MLE
facs$loadings		# The Bs
facs$scores				# The scores, the Ft's



indus<-read.csv("47indus-day.csv")
nind <-length(indus[1,])-3;nind 
rets <-cbind(indus[,2:(nind+1)])
pcs<-princomp(rets)


(retmat%*%pcs$loadings[,2])-pcs$scores[,2]

