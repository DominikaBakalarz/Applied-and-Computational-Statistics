install.packages('resampledata')  
library('resampledata')
# Load Earthquakes data
x = Quakes[,2]
theta_MLE <- 1/mean(x)
theta_MLE

#1. bootstrap estimate of V(theta_MLE)
bootstrap_variance <- function(X, B) {
  n <- length(X)
  mhat.boot <- numeric(B) 
  for (j in 1:B){
    X.boot <- rexp(n, theta_MLE)
    mhat.boot[j] <- 1/mean(X.boot) # Median bootstrap samples 
    }
    var.boot <- var(mhat.boot) # Evaluate the bootstrap variance estimate
  return(list(var.boot = var.boot, mhat.boot = mhat.boot)) 
  }
var.boot <- bootstrap_variance(x,1000)[[1]] 
var.boot #3.693184e-06
mhat.boot <- bootstrap_variance(x,1000)[[2]] 
mhat.boot

#2. parametric estimates of the quantiles
q1 = qexp(0.1, rate = theta_MLE)
q2 =qexp(0.25, rate = theta_MLE, lower.tail = TRUE, log.p = FALSE)
q3 =qexp(0.5, rate = theta_MLE, lower.tail = TRUE, log.p = FALSE)
q4 =qexp(0.75, rate = theta_MLE, lower.tail = TRUE, log.p = FALSE)
q5 =qexp(0.9, rate = theta_MLE, lower.tail = TRUE, log.p = FALSE)

c(q1,q2,q3,q4,q5)

#3. Using the parametric bootstrap, provide 99% approximate confidence intervals for q(α) - normal CI or bootstrap pivotal CI?
bootstrap_quantiles <- function(X, B, alpha) {
  n <- length(X)
  mhat.boot <- numeric(B) 
  for (j in 1:B){
    X.boot <- rexp(n, theta_MLE) # Sample bootstrap data from Fn
    mhat.boot[j] <- quantile(X.boot, alpha, names = FALSE) # Median bootstrap samples 
  }
  var.boot <- var(mhat.boot) # Evaluate the bootstrap variance estimate
  return(list(var.boot = var.boot, mhat.boot = mhat.boot)) 
}

#normal confidence intervals:
beta = 0.01
var.boot <- bootstrap_quantiles(x,1000, 0.1)[[1]] 
CI=c(qexp(0.1, rate = theta_MLE) - qnorm(1-beta/2)*sqrt(var.boot),qexp(0.1, rate = theta_MLE) + qnorm(1-beta/2)*sqrt(var.boot))
CI #1.366575 2.434573
var.boot <- bootstrap_quantiles(x,1000, 0.25)[[1]] 
CI=c(qexp(0.25, rate = theta_MLE) - qnorm(0.995)*sqrt(var.boot),qexp(0.25, rate = theta_MLE) + qnorm(0.995)*sqrt(var.boot))
CI #4.239425 6.139435
var.boot <- bootstrap_quantiles(x,1000, 0.5)[[1]] 
CI=c(qexp(0.5, rate = theta_MLE) - qnorm(0.995)*sqrt(var.boot),qexp(0.5, rate = theta_MLE) + qnorm(0.995)*sqrt(var.boot))
CI #10.81591 14.19113
var.boot <- bootstrap_quantiles(x,1000, 0.75)[[1]] 
CI=c(qexp(0.75, rate = theta_MLE) - qnorm(0.995)*sqrt(var.boot),qexp(0.75, rate = theta_MLE) + qnorm(0.995)*sqrt(var.boot))
CI #22.10074 27.91335
var.boot <- bootstrap_quantiles(x,1000, 0.9)[[1]] 
CI=c(qexp(0.9, rate = theta_MLE) - qnorm(0.995)*sqrt(var.boot),qexp(0.9, rate = theta_MLE) + qnorm(0.995)*sqrt(var.boot))
CI #36.70103 46.37057

#pivotal confidence intervals:
beta = 0.01
mhat.boot <- bootstrap_quantiles(x,1000, 0.1)[[2]] 
ci.bootpivot = c(2*q1 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*q1 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles(x,1000, 0.25)[[2]] 
ci.bootpivot = c(2*q2 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*q2 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles(x,1000, 0.5)[[2]] 
ci.bootpivot = c(2*q3 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*q3 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles(x,1000, 0.75)[[2]] 
ci.bootpivot = c(2*q4 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*q4 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles(x,1000, 0.9)[[2]] 
ci.bootpivot = c(2*q5 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*q5 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot

#4. Discuss which of the parametric or nonparametric approach is more suitable here
hist(x, freq=FALSE, xlab = "Time between earthquakes (in days)", xlim=c(min(x),max(x)),  ylim=c(0, 0.05), main="")
x0 <- seq(0, 140, length=805)
hx <- dexp(x0, rate = theta_MLE)
par(new=TRUE)
plot(x0, hx, type="l", lty=2, xlab="",
     ylab="", xlim=c(min(x),max(x)),  ylim=c(0, 0.05))
#qqplot
plot(qexp(ppoints(length(x))),sort(x), xlab="Exponential Quantiles", ylab="Observed values")
abline(a=0, b =1/theta_MLE,lty=2)

#.5. prediction on X_(n+1):
bootstrap_h <- function(X, B, q) {
  n <- length(X)
  X.boot <- numeric(B) 
  #X.boot <- matrix(, nrow = B, ncol = n)
  for (j in 1:B){
    #X.boot[j,1:n] <- rexp(n, mhat.boot[j])
    X.boot[j] <- rexp(1, mhat.boot[j])
    #p.boot[j] <- sum(X.boot < q)/n
  }
  p.boot <- sum(X.boot < q)/B
  #p.boot <- sum(X.boot < q)/(n*B)
  mean.p.boot <- mean(p.boot) # Evaluate the bootstrap variance estimate
  return(list(mean.p.boot = mean.p.boot, p.boot = p.boot)) 
}

#h(alpha) for some range of values alpha
mhat.boot <- bootstrap_variance(x,1000)[[2]] 
mhat.boot
q = qexp(0.1, rate = theta_MLE)
h_hat <- bootstrap_h(x, 1000, q)[[2]]
h_hat
q = qexp(0.25, rate = theta_MLE)
h_hat <- bootstrap_h(x, 1000, q)[[1]]
h_hat
q = qexp(0.5, rate = theta_MLE)
h_hat <- bootstrap_h(x, 1000, q)[[1]]
h_hat
q = qexp(0.75, rate = theta_MLE)
h_hat <- bootstrap_h(x, 1000, q)[[1]]
h_hat
q = qexp(0.9, rate = theta_MLE)
h_hat <- bootstrap_h(x, 1000, q)[[1]]
h_hat

#confidence interval for prediction
mean(bootstrap_quantiles(x,1000, 0.005)[[2]])
mean(bootstrap_quantiles(x,1000, 0.995)[[2]])

#alternative approach
bootbootstrap_quantiles <- function(X, B, alpha) {
  n <- length(X)
  m.boot <- numeric(B) 
  for (j in 1:B){
    X.boot <- rexp(n, mhat.boot[j]) # Sample bootstrap data from Fn
    m.boot[j] <- quantile(X.boot, alpha, names = FALSE) # Median bootstrap samples 
  }
  var.boot <- var(m.boot) # Evaluate the bootstrap variance estimate
  return(list(var.boot = var.boot, m.boot = m.boot)) 
}

mean(bootbootstrap_quantiles(x,1000, 0.005)[[2]])
mean(bootbootstrap_quantiles(x,1000, 0.995)[[2]])


x = Service[,2]
hist(x,xlab = "Customer service time (in minutes)", freq = FALSE, xlim= c(min(x), max(x)), ylim=c(0,1.2), main="")
x0 <- seq(min(x), max(x), length=174)
hx <- dgamma(x0, shape=shape, scale=scale)
par(new=TRUE)
plot(x0, hx, type="l", lty=2, xlab="", ylab="", xlim= c(min(x), max(x)), ylim=c(0,1.2))
install.packages('dglm')  
library(dglm)
fit <- dglm(x~1, family=Gamma(link="log"), mustart=mean(x))
summary(fit)

mu <- exp(fit$coefficients)
shape <- exp(-fit$dispersion.fit[[1]])
scale <- mu/shape
c(shape, scale)
rate <- 1/scale

#2. parametric estimates of the quantiles
Q1= qgamma(0.1, shape, scale = scale)
Q2= qgamma(0.25, shape, scale = scale)
Q3= qgamma(0.5, shape, scale = scale)
Q4= qgamma(0.75, shape, scale = scale)
Q5= qgamma(0.9, shape, scale = scale)
c(Q1,Q2,Q3,Q4,Q5)
#3. CIs
bootstrap_quantiles2 <- function(X, B, alpha) {
  n <- length(X)
  mhat.boot <- numeric(B) 
  for (j in 1:B){
    #X.boot <- sample(X,n,replace=TRUE) # Sample bootstrap data from Fn
    #fit <- dglm(X.boot ~ 1, family=Gamma(link="log"), mustart=mean(X.boot))
    X.boot <- rgamma(n, shape = shape, scale = scale)
    mhat.boot[j] <- quantile(X.boot, alpha, names = FALSE)
    #mhat.boot[j] <- qgamma(alpha, shape = exp(-fit$dispersion.fit[[1]]), scale = exp(fit$coefficients)/shape) # Median bootstrap samples 
  }
  var.boot <- var(mhat.boot) # Evaluate the bootstrap variance estimate
  return(list(var.boot = var.boot, mhat.boot = mhat.boot)) 
}

#pivotal confidence intervals:
mhat.boot <- bootstrap_quantiles2(x,1000, 0.1)[[2]] 
ci.bootpivot = c(2*Q1 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*Q1 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles2(x,100, 0.25)[[2]] 
ci.bootpivot = c(2*Q2 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*Q2 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles2(x,100, 0.5)[[2]] 
ci.bootpivot = c(2*Q3 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*Q3 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles2(x,100, 0.75)[[2]] 
ci.bootpivot = c(2*Q4 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*Q4 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot
mhat.boot <- bootstrap_quantiles2(x,100, 0.9)[[2]] 
ci.bootpivot = c(2*Q5 - quantile(mhat.boot, 1-beta/2,names=FALSE) , 2*Q5 - quantile(mhat.boot, beta/2,names=FALSE))
ci.bootpivot

#confidence interval for prediction
mean(bootstrap_quantiles2(x,1000, 0.005)[[2]])
mean(bootstrap_quantiles2(x,1000, 0.995)[[2]])

