STUDENT-T CODELikelihood Function#  Function to compute the likelihood of data(TT) ~ sd * Student(df)t.loglik <- function(data,TT,df,sd){loglik <- TT*df*log(sd) + sum(log(dt(data,df)))loglik <- loglik - 0.5*(df+1)*sum(log((df*sd^2+data^2)/(df+data^2)))loglik	}

# Functions R should have !

colmax<- function(X){apply(X,2,max)}
rowmax<- function(X){apply(X,1,max)}


Simulating Student-t Data# Simulate TT data Student-t(nudraw,sdraw)TT 	 <- 500 nu	 <- 10				# Try low and high nu.sig  <- 0.2/sqrt(52)	# Weekly US index standard deviationrm(xx)xx    <- rt(TT,nu)*sig		# simulated data##Plotting the Likelihood## Grid of values for likelihood plot# Making sure we cover large enough range, this is overkillshat <- sd(xx);shat/sig				# sample standard deviation too highshigh<- shat * sqrt(qchisq(0.99,TT)/TT) # bounds for gridslow <- shat * sqrt(qchisq(0.01,TT)/TT) / sqrt(4/(4-2))  nuvec	<- seq(4, 53,length=50)#sigvec	<- seq(slow,shat, length=35)		# bottom half#sigvec	<- c(sigvec,seq(shat+(shigh-shat)/(35-1),shigh,length=35)) # add top
sigvec <-seq(slow,shigh,length=70)
## Plot the likelihood vs. Nu and sigma#

rm(loglikmat); loglikmat<- matrix(0,nrow=50,ncol=70)inu<-1;  isig<-1;  i<-1for (inu in 1:50){  for (isig in 1:70){loglikmat[inu,isig]<-t.loglik(xx,TT,nuvec[inu],sigvec[isig])   }}contour(nuvec,sigvec,loglikmat,nlevels=50,xlab="nu",ylab="sigma")abline(v= nuvec[which(loglikmat==max(loglikmat),arr.ind=T)[1]])abline(h= sigvec[which(loglikmat==max(loglikmat),arr.ind=T)[2]])persp(nuvec,sigvec,loglikmat,theta=-20,phi=25)

par(mfrow=c(2,1))plot(sigvec*sqrt(52),colmax(loglikmat),type="l",xlab="sigma",ylab="Log-Lik")plot(nuvec,rowmax(loglikmat),type="l",xlab="nu",ylab="Log-Lik")#persp3d(nuvec,sigvec,loglikmat,col="green")  # Need rgl##Estimate parameters by MLE#
t.loglik <- function(data,TT,df,sd){
loglik <- TT*df*log(sd) + sum(log(dt(data,df)))
loglik <- loglik - 
0.5*(df+1)*sum(log((df*sd^2+data^2)/(df+data^2)))
loglik	}

t.loglikmle <- function(p,xx,TT){		# routine computes minus the likelihoodloglik <- TT*p[2]*log(p[1]) + sum(log(dt(xx,p[2])))loglik <- loglik - 0.5*(p[2]+1)*sum(log((p[2]*p[1]^2+xx^2)/(p[2]+xx^2)))loglik <- -loglikloglik   }# MLE in R: 	nlminb (no numerical Hessian)#			optim  (   numerical Hessian)#			fitdistr in MASS uses optim# ind<-which(loglikvec==max(loglikvec),arr.ind=T)xx    <- rt(TT,nu)*sig		# simulated dataout1<-nlminb(c(shat,30),t.loglikmle,xx=xx,TT=TT,lower=c(0,3),upper=c(Inf,100))out1$parout2<-optim(c(shat,20),xx=xx,TT=TT,t.loglikmle, lower=c(slow,4), upper=c(10,50), hessian=T,method = "L-BFGS-B")out2$value; out2$par; sqrt(diag(solve(out2$hessian)))###############################################################
fitdistr(rt(500,5)*0.2/sqrt(52),"t",list(m=0,s=1,df=20))?Prior on Sigma		Better use proper conjugate prior: Inverted Gamma		Calibrate so that it is uninformative if desiredR: 	Gamma(shape, scale) 	~ xshape-1 exp(-x / scale)	Then ? = 1/sqrt(x) 	~ Inverted Gamma  ~  exp(-1/scale ?2) / ? 2 shape+1	shape  = nu/2	scale    = 2 / (nu*ssq)TT<-1000; 	mu<-0.1/12; 	sig<-0.2/sqrt(12)rt         <- rnorm(TT, mu,sig)muhat<-mean(rt) nu    <- TT-1ssq   <-sum(scale(rt,scale=F)^2)sh     <- nu/2sc     <- 2/ (nu * ssq)sigdraw <- rgamma(1000,shape=sh,scale=sc)sigdraw <- 1/sqrt(sigdraw)mean(sigdraw); quantile(sigdraw,c(0.025,0.975))2