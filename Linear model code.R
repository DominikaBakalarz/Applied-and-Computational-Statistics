#Linear Model - swim data

#Step 1 - load data and prepare it for the analysis
swim <- read.table("http://www.stats.ox.ac.uk/~laws/LMs/data/swim.txt", header = TRUE)
swim$Shirt <- as.factor(swim$Shirt)
swim$Goggles <- as.factor(swim$Goggles)
swim$Flippers <- as.factor(swim$Flippers)
swim$End <- as.factor(swim$End)

#basic information about or data
str(swim)
head(swim)
summary(swim)

#Step 2 - Exploratory data analysis
attach(swim)
par(mfrow=c(2,2))
boxplot(Time~Shirt,data=swim, xlab="shirt", ylab="time")
boxplot(Time~Goggles,data=swim, 
        xlab="goggles", ylab="time")
boxplot(Time~Flippers,data=swim, 
        xlab="flippers", ylab="time")
boxplot(Time~End,data=swim, 
        xlab="end", ylab="time")
plot(Time~Order,data=swim, 
        xlab="order", ylab="time")

#interaction plots
par(mfrow=c(2,3))
interaction.plot(x.factor=Goggles,trace.factor=Shirt,
                 response=Time,fixed=T)
interaction.plot(x.factor=Flippers,trace.factor=Shirt,
                 response=Time,fixed=T)
interaction.plot(x.factor=End,trace.factor=Shirt,
                 response=Time,fixed=T)
interaction.plot(x.factor=Goggles,trace.factor=Flippers,
                 response=Time,fixed=T)
interaction.plot(x.factor=Goggles,trace.factor=End,
                 response=Time,fixed=T)
interaction.plot(x.factor=End,trace.factor=Flippers,
                 response=Time,fixed=T)

#Step 3 - Model Selection
swim.lm <- lm(Time ~ Order*Shirt*Goggles*Flippers*End, data = swim)
summary(swim.lm) #NA producded - can't use this model
swim2.lm <- lm(Time ~ (Order + Shirt + Goggles + Flippers + End)^3, data = swim)
summary(swim2.lm) #NA producded - can't use this model
swim3.lm <- lm(Time ~ (Order + Shirt + Goggles + Flippers + End)^2, data = swim)
summary(swim3.lm) #starting point for backwards elimination search via ANOVA
#it is called model(*) in the report

#Step 3a - outliers analysis
par(mfrow=c(1,1))
plot(cooks.distance(swim3.lm), 
     main = "swim data", ylab = "Cook's distance")
identify(1:24, cooks.distance(swim3.lm)) 
# (8/(n - 2*p))) is a useless bound in this case as it is negative
#but we can see on the plot that points 9 and 10 are potential outliers
i <- (cooks.distance(swim3.lm) > 0.8)
par(mfrow=c(2,2))
plot(rstandard(swim3.lm) ~ fitted(swim3.lm), main = "Standardised residuals vs Fitted Values",
     xlab = "Fitted values", ylab = "Residuals", lower.panel = NULL, col = 1+i, pch = 1 + 15*i)
qqnorm(resid(swim3.lm), main = "Q-Q plot of (Studentised) residuals", col = 1+i, pch = 1 + 15*i )
qqline(resid(swim3.lm))
plot(hatvalues(swim3.lm), 
     main = "Leverage ", ylab = "Leverage",col = 1+i, pch = 1 + 15*i )
plot(cooks.distance(swim3.lm), 
     main = "Cook's distance", ylab = "Cook's distance", col = 1+i, pch = 1 + 15*i )
swim2 <- swim[-which(i), ] #We remove the outliers
attach(swim2)

#search for potential outliers again:
swim4.lm <- lm(Time ~ (Order + Shirt + Goggles + Flippers + End)^2) #we need a new linear model 
summary(swim4.lm)
par(mfrow=c(2,2))
plot(rstandard(swim4.lm) ~ fitted(swim4.lm), main = "Standardised residuals vs Fitted Values",
     xlab = "Fitted values", ylab = "Residuals", lower.panel = NULL)
qqnorm(resid(swim4.lm), main = "Q-Q plot of (Studentised) residuals" )
qqline(resid(swim4.lm))
plot(hatvalues(swim4.lm), 
     main = "Leverage ", ylab = "Leverage" )
plot(cooks.distance(swim4.lm), 
     main = "Cook's distance", ylab = "Cook's distance")
#it looks better now - we do not remove any other points

#Step 3b - model selection
summary(swim4.lm)
#notice here that no term with End is significant, try removing it:
swim5.lm <- lm(Time ~ (Order + Shirt + Goggles + Flippers)^2)
anova(swim4.lm,swim5.lm) #Direct F-test for null model 4 aganst alternative model 5
#p-value not significant, we can simplify the model
summary(swim5.lm)
swim6.lm <- lm(Time ~ (Shirt + Goggles + Flippers)^2)
summary(swim6.lm)
anova(swim5.lm,swim6.lm) #Direct F-test for null model 5 aganst alternative model 6
# p-value is siginifcant, we can't remove order completely
swim7.lm <- lm(Time ~ (Order + Goggles + Flippers)^2)
summary(swim7.lm)
anova(swim5.lm,swim7.lm) #Direct F-test for null model 5 aganst alternative model 7
# p-value is very siginifcant, we can't remove shirt completely
swim8.lm <- lm(Time ~ Order*Shirt + Goggles*Flippers)
confint(swim8.lm)
summary(swim8.lm)
anova(swim5.lm,swim8.lm) #Direct F-test for null model 5 aganst alternative model 8
#p-value not significant, we can simplify the model
summary(swim8.lm) 
#outliers analysis
par(mfrow=c(2,2))
plot(rstandard(swim8.lm) ~ fitted(swim8.lm), main = "Standardised residuals vs Fitted Values",
     xlab = "Fitted values", ylab = "Residuals", lower.panel = NULL)
qqnorm(resid(swim8.lm), main = "Q-Q plot of (Studentised) residuals" )
qqline(resid(swim8.lm))
plot(hatvalues(swim8.lm),  main = "Leverage ", ylab = "Leverage" )
plot(cooks.distance(swim8.lm), 
     main = "Cook's distance", ylab = "Cook's distance")
#looks good
