#Step 1 - load data and get basic information
apples <- read.table("apple.csv", header = TRUE, sep = ",")
sum(apples$growth == 1) #how many entires have growth present - 32
sum(apples$growth == 0) #how many entires have growth absent - 42
attach(apples)


par(mfrow=c(2,2))
boxplot(ph~growth,xlab='Growth',ylab='pH level')
boxplot(nisin~growth,xlab='Growth',ylab='Nisin level')
boxplot(temperature~growth,xlab='Growth',ylab='temperature')
boxplot(brix~growth,xlab='Growth',ylab='Brix level')

#rows sum to 1
par(mfrow=c(2,2))
barplot(prop.table(table(growth,ph),1),beside=T,col=c(4,2),
        xlab='pH level',ylab='Proportion')
barplot(prop.table(table(growth,nisin),1),beside=T,col=c(4,2),
        xlab='Nisin level',ylab='Proportion')
barplot(prop.table(table(growth,temperature),1),beside=T,col=c(4,2),
        xlab='Temperature',ylab='Proportion')
barplot(prop.table(table(growth,brix),1),beside=T,col=c(4,2),
        xlab='Brix level',ylab='Proportion')

#columns sum to 1
par(mfrow=c(2,2))
barplot(prop.table(table(growth,ph),2),beside=T,col=c(4,2),
        xlab='ph',ylab='Proportion')
barplot(prop.table(table(growth,nisin),2),beside=T,col=c(4,2),
        xlab='nisin',ylab='Proportion')
barplot(prop.table(table(growth,temperature),2),beside=T,col=c(4,2),
        xlab='temperature',ylab='Proportion')
barplot(prop.table(table(growth,brix),2),beside=T,col=c(4,2),
        xlab='brix',ylab='Proportion')


#univariate logistic regression
a1.glm <- glm(growth~ph, family = binomial() )
summary(a1.glm)
a2.glm <- glm(growth~nisin, family = binomial() )
summary(a2.glm)
a3.glm <- glm(growth~temperature, family = binomial() )
summary(a3.glm)
a4.glm <- glm(growth~brix, family = binomial() )
summary(a4.glm)



#model selection
a5.glm <- glm(growth~ph*nisin*temperature*brix, family = binomial ) #M1
summary(a5.glm)
#we get a warning, no p-value is significant 

a6.glm <- glm(growth~(ph+nisin+temperature+brix)^3, family = binomial ) #M2
summary(a6.glm)
#we get a warning
#can we drop any terms from this model?

#removing 3-way interactions:
a7.glm <- glm(growth~(ph+nisin+temperature+brix)^2, family = binomial ) #M3
summary(a7.glm)

#comparing a6 and a7
1-pchisq(66.05 - 39.525, 4)
# p value super small

#removing nisin:temp interactions, as none of nisin:temp, nisin:temp:ph, nisin:temp:brix is significant
a8.glm <- glm(growth~(ph+nisin+brix)^3+(ph+temperature+brix)^3, family = binomial )
summary(a8.glm)
#looks better, only ph, brix, ph:brix not significant at 0.001 level, but ph:brix:temp and ph:brix:nisin are siginificant

#comparing a6 and a8
1-pchisq(45.517 - 39.525, 3)
#p-value 0.1120001
# we can use it as our final model 

#outlier analysis
par(mfrow=c(1,2))

plot((a8.glm)$fitted.values,rstandard(a8.glm),xlab='Fitted Values', 
     ylab='Deviance Residuals',pch=19,col='blue')
identify(a8.glm$fitted.values, rstandard(a8.glm))
abline(h=0)

plot(cooks.distance(a8.glm),col='blue',
     pch=19,ylab="Cook's Distance")
identify( cooks.distance(a8.glm))
abline(h=8/(74-24))
#49 is definitely an outlier

apples2 <- apples[-49,]
a10.glm <- glm(growth~(ph+nisin+brix)^3+(ph+temperature+brix)^3, family = binomial, data = apples2 )
summary(a10.glm)
#all values are so small now! It's not a good model, go back to model selection

#drop all three-way interactions
a11.glm <- glm(growth~(ph+nisin+brix+temperature)^2, family = binomial, data = apples2 )
summary(a11.glm)

#drop all two-way interactions
a12.glm <- glm(growth~(ph+nisin+brix+temperature), family = binomial, data = apples2 )
summary(a12.glm)

#drop temperature?
a13.glm <- glm(growth~(ph+nisin+brix), family = binomial, data = apples2 )
summary(a13.glm)

#test for dropping temperature
1-pchisq(72.539 - 68.532, 1)
# p-value: 0.04531171 - significant, so we won't drop temperature

par(mfrow=c(2,2))
# LEVERAGE
plot(influence(a12.glm)$hat/(5/73),col='blue', pch=19,ylab='Leverage/(p/n)')
abline(h=2)
identify(influence(a12.glm)$hat/(5/73))


# DEVIANCE RESIDUALS
plot((a12.glm)$fitted.values,rstandard(a12.glm),xlab='Fitted Values',
     ylab='Deviance Residuals',pch=19,col='blue')
abline(h=0)

#qq plots not applicable to Bernoulli

# COOKS DISTANCE
plot(cooks.distance(a12.glm),ylim=c(0,0.14),col='blue',
     pch=19,ylab="Cook's Distance")
abline(h=8/63)

# WORKING RESIDUALS
plot(a12.glm$fitted.values, a12.glm$residuals,xlab='Fitted Values',
     ylab='Working Residuals',pch=19,col='blue')
abline(h=0)



#THE END




identify(a8.glm$fitted.values, a8.glm$residuals)


# LEVERAGE
plot(influence(a8.glm)$hat/(12/74),col='blue', pch=19,ylab='Leverage/(p/n)')

#CHANGE 12/74
plot(influence(a12.glm)$hat/(12/74),col='blue', pch=19,ylab='Leverage/(p/n)')

# DEVIANCE RESIDUALS
plot((a8.glm)$fitted.values,rstandard(a8.glm),xlab='Fitted Values',
     ylab='Deviance Residuals',pch=19,col='blue')
identify(a8.glm$fitted.values, rstandard(a8.glm))
abline(h=0)

qqnorm(rstandard(a8.glm),pch=19,main='')
qqline(rstandard(a8.glm))
#new model
plot((a12.glm)$fitted.values,rstandard(a12.glm),xlab='Fitted Values',
     ylab='Deviance Residuals',pch=19,col='blue')
abline(h=0)

#not for bernoulli
qqnorm(rstandard(a12.glm),pch=19,main='')
qqline(rstandard(a12.glm))

# COOKS DISTANCE
plot(cooks.distance(a8.glm),col='blue',
     pch=19,ylab="Cook's Distance")
identify( cooks.distance(a8.glm))
#49 is an outlier

plot(cooks.distance(a10.glm),col='blue',
     pch=19,ylab="Cook's Distance")
# everything is so small


#other plots
a9.glm <- glm(growth~ph + nisin + brix + temperature, family = binomial )
summary(a9.glm)
plot(a9.glm$fitted.values, a9.glm$residuals,xlab='Fitted Values',
     ylab='Working Residuals',pch=19,col='blue')
abline(h=0)
