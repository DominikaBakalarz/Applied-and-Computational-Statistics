library(MASS)

#Part1
#Step 1 - load data and get basic information
blood <- read.table('http://www.stats.ox.ac.uk/~nicholls/CompStats/home.txt',header=T)
attach(blood)
par(mfrow=c(1,1))
plot(Subject, home ,xlab='Subject',ylab='blood pressure',pch =16, ylim=c(65,110))
par(new=TRUE)
plot(Subject, hospital,ylab='blood pressure' ,col='black', ylim=c(65,110))
legend(1, 110, legend=c("home", "hospital"),
       col=c("black", "black"), pch=c(16,1), cex=0.8)
#add some boxplots and histograms:
par(mfrow=c(1,3))
boxplot(home,xlab='home',ylab='blood pressure', ylim = c(73,107))
boxplot(hospital,xlab='hospital',ylab='blood pressure', ylim = c(73,107))
boxplot(hospital - home,xlab='hospital - home',ylab='difference in blood pressure')

#SWAP to get hospital - home
#Step 2 - test to see if the home and hospital measurements have the same distribution 
diff <- hospital - home
mean(diff) #5.18
median(diff) #7
#fitdistr(diff, "normal") #mean 5.8
sort(diff) # we can see that there are no ties, so we can use wilcoxon.test
wilcox.test(hospital, home, paired = TRUE ,mu=0,conf.int=TRUE) #two-sided alternative, because we are given no expert knowledge on what to expect
#p-value = 0.01855, estimate = 6, reject H0

#Step 3 - test to see if there is an offset of 6 units
home2 <- home + 6
diff2 <- hospital - home2
sort(diff2) #no ties, but there's 0 in the data, and ties in the abs values. so can't use wilcoxon.test
#wilcox.test(home2, hospital, paired = TRUE ,mu=0,conf.int=TRUE)
diff <- diff2[diff2!=0] #remove 0

srt <- data.frame(diff)
srt$adiff <- abs(srt$diff)
srt$rank <- rank(srt$adiff)
srt$srank <- srt$rank*sign(srt$diff)
W <- sum(srt$srank[srt$srank>0])
W #26

get_pvalue_srt <- function(w_observed, n, nsim=100000){
  w_rand <- numeric(nsim)
  rnk<-1:n
  for (sim in 1:nsim){
    sign <- sample(c(-1,1), replace=TRUE, size=n)
    srank <- rnk*sign
    w_rand[sim] <- sum(srank[srank>0])
  }
  
  prob_larger <- mean(w_observed <= w_rand)
  prob_smaller <- mean(w_observed >= w_rand)
  p_value <- 2*min(prob_smaller, prob_larger)
  return(list(prob_smaller = prob_smaller, prob_larger = prob_larger, p_value = p_value, w = w_rand))
}

sim <- get_pvalue_srt(W,n=length(diff2!=0)) 
sim$p_value #0.77, no evidence to reject H0, but collecting more data would be good


#Part2
#Step 1 - load data and get basic information
clotting <- read.table('http://www.stats.ox.ac.uk/~nicholls/CompStats/clotting.txt',header=T)
attach(clotting)
par(mfrow=c(1,1))
plot(Subject, old ,xlab='Subject',ylab='blood clotting time',pch =16 , ylim=c(0,650))
par(new=TRUE)
plot(Subject, new, ylab='blood clotting time' ,col='black', ylim=c(0,650))
legend(1, 610, legend=c("old", "new"),
       col=c("black", "black"), pch=c(16,1), cex=0.8)

par(mfrow=c(1,3))
boxplot(old,xlab='old drug',ylab='blood clotting time', ylim=c(60,610))
boxplot(new,xlab='new drug',ylab='blood clotting time', ylim=c(60,610))
boxplot(old - new,xlab='old - new',ylab='difference in the blood clotting time')

#Step 2 - test to see if the blood clotting times for the old and new drugs have the same distribution 
diff <- old - new
sort(diff) # we can see that there are no ties, BUT THERE ARE TIES IN THE ABS VALUES, so can't use wilcoxon.test

srt <- data.frame(diff)
srt$adiff <- abs(srt$diff)
srt$rank <- rank(srt$adiff)
srt$srank <- srt$rank*sign(srt$diff)
W <- sum(srt$srank[srt$srank>0])
W #97.5

sim <- get_pvalue_srt(W,n=length(old))
sim$prob_larger #0.01481, one-sided test, reject  H0 

#Step 3 - confidence interval for the change
X=srt$diff
n=length(X)
A=matrix(NA,n,n,dimnames=list(paste(X),paste(X)));
for (i in 1:n) {for (j in i:n) {A[i,j]=(X[i]+X[j])/2}}
#the n(n+1)/2=120 sorted Walsh averages
v=sort(A[upper.tri(A,diag=TRUE)]); v
round(psignrank(1:(n*(n+1)/2),n),4)
# these are alpha/2 values available, W<=25 is 0.0240
#take alpha/2=0.0240 so alpha=0.048 <= 0.05
psignrank(25,n);  #so the lower end of the CI is at W=26
n*(n+1)/2-26 #94
psignrank(116,n); #upper end at W=116,  
#We want the 26th and 94th of the Walsh averages
v[c(26,94)] #( 7.5, 210.0)
