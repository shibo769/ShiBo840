#
# Hotelling T2 (Wald) interval  
# Q1: Compare cutoffs from the exact F(m,T-m) distribution
#		vs approximate Chisquare(m) distribution
#

TT <-120
mm <-seq(1,20)
fac<-(TT-1)*mm/(TT-mm)

par(mgp=c(1.5,0.5,0))
plot( fac*qf(0.95,mm,TT-mm),qchisq(0.95,mm), ylim=c(2,41),
ylab=bquote(paste(chi^2,"(m) 95%")),
xlab=bquote(paste("(",.(TT),"-1)/(",.(TT),"-m) * m F(m,",.(TT),"-m) 95%")),type="o")
abline(0,1,col="red")
title(bquote(paste("Wald test: Exact F based vs Approx. ",chi^2,"cutoff, m = 1, .. , 20")),line=1)

# One-at-a-time (Wrong) vs simultaneous Hotelling T2 cutoffs
# Very large intervals 
# Discussed in class in details

TT     <-120
mm     <- seq(1,20)	# extent of data snooping
fac	   <-(TT-1)*mm/(TT-mm)

tcut   <- qt(1-0.05/2,TT-1)
Tsqcut <- qf(0.95,mm,TT-mm)*mm*(TT-1)/(TT-mm)
plot(mm,sqrt(Tsqcut),type="l",xlab="dimension",
ylab="Tsquare cutoff")
title(bquote(paste("Hotelling F-based cutoff vs dimension m, T=",.(TT))),line=0.5)


# Assuming independent components
# Actual probability of one rejection (actual test alpha)
# if using one-at-a-time (wrong) t-stats

reject<- 1-(1-0.05)^mm		# much larger than 0.05
plot(mm,reject,type="l",ylab="Pr(rejection)")
title(bquote(paste("Independent components, Pr(reject) under H0, nominal ",alpha,
"=0.05")),line=0.6)

# Assuming independent components
# Correcting indidual t alphas for simultaneous 
# (data mining) test to have size alpha=0.05

alphcorr <- 1-0.95^(1/mm)
indepcut <- qt(1-alphcorr/2,TT-1)
plot(mm,indepcut,type="l",ylab="cutoff tstat with corrected alpha")
title("Correcting alpha with independent components 
to keep Pr(reject)=0.05")


# Bonferroni cutoffs based on the number of dimensions
# no correlation assumptions 

alphbonf<- 0.05/mm
bonfcut <- qt(1 - alphbonf/2,TT-1)
plot(mm,indepcut,xlab="dimension",ylab="Bonferroni",type="l")
title("Bonferroni t cutoff vs number of tests",line=0.2)

lines(mm,bonfcut,col="red")
legend("topleft",c("independent","Bonferroni"),lty=c(1,1),col=c("black","red"),bty="n")



