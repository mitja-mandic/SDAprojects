# library(ggplot2)
# library(tidyverse)
# library(reshape2)
library(rrcov)
library(robustbase)
library(corrplot)
library(MASS)
library(corrplot)
library(stats)
library(car)
rm(list = ls())

cars_data <- read.table("Project3_data.txt", sep="", header=T)
cars_data$euro_standard <- as.factor(cars_data$euro_standard)
cars_data$transmission_type <- as.factor(cars_data$transmission_type)
cars_data$fuel_type <- as.factor(cars_data$fuel_type)
set.seed(0875362)
data_ind <- sample.int(n=nrow(cars_data), size=500, replace=F)
mydata <- cars_data[data_ind, ]
attach(mydata)

# QUESTION 1

#fit full model
cars.full <- lm(co2~., data = mydata[,-c(1,2,3)])

#insignificant: co_emissions, noise_level, extra_urban_metric, urban_metric, euro_standard4 - euro_standard is factor
cars.full
summary(cars.full)

#multicoll
#remove categorical variables euro_standard, transmission_type, fuel_type and descriptives for correlation plots
mydata.noDescriptives <- within(mydata, rm(manufacturer, model, description, euro_standard, transmission_type, fuel_type))
mydata.noDescriptives <- mydata.noDescriptives[,c(6,1,2,3,4,5,7,8)]
#no response
corrMat <- cor(mydata.noDescriptives[,-1])
round(cor(mydata.noDescriptives)[,"co2"], 2)
det(corrMat)

corrplot(corrMat, tl.col = "black")
sum(cars.full$residuals)
pairs(mydata.noDescriptives)

#vif values for our model. Metrics by far the highest, others relatively low
vif(cars.full)


#standardized residuals. a few outliers, but variance seems homogenous.
cars.full.res <- stdres(cars.full)
hist(cars.full.res)

cutoff <- qnorm(0.995)
plot(cars.full.res, xlab = "Index", ylab = "", pch = 19,
     ylim = c(min(cars.full.res, -cutoff), max(cars.full.res, cutoff)), main = "Standardized residuals")
abline(h = 0, col = "grey", lwd = 1.5, lty = 2)
abline(h = -cutoff, col = "red", lwd = 1.5)
abline(h = cutoff, col = "red", lwd = 1.5)

qqnorm(cars.full.res, ylab = "Standardized residuals", main = "Normal Q-Q Plot of standardized residuals")

#residuals vs fitted values
plot(cars.full$fitted.values, cars.full$residuals, pch=19, xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs fitted values")
abline(h = 0, col = "grey", lwd = 1.5, lty = 2)

mean(cars.full.res%*%t(cars.full.res))

#check the model in general
par(mfrow=c(2,2))
plot(cars.full)
par(mfrow=c(1,1))

######### QUESTION 2 ######### 
######transforming the response

hist(co2)
co2.bc <- boxCox(co2~1)
lambda <- co2.bc$x[which.max(co2.bc[["y"]])]

#####Select transformation here!:

#trans.co2 <- bcPower(co2, lambda)
#trans.co2 <- log(co2)
#trans.co2 <- sqrt(co2)

qqnorm(trans.co2)
hist(trans.co2)

tmydata <- within(mydata, {tco2 <- trans.co2; rm(co2)})


t.cars.full <- lm(tmydata$tco2~., data = tmydata[,-c(1,2,3)])
summary(t.cars.full)
t.cars.full.res <- stdres(t.cars.full)
sum(t.cars.full$residuals) #sum is zero in all cases
plot(t.cars.full.res, xlab = "Index", ylab = "Standardized residual",
     main = "Standardized residuals for the square root of the response", 
     pch = 19,
     ylim = c(min(t.cars.full.res, -cutoff), max(t.cars.full.res, cutoff)))

abline(h = 0, col = "grey", lwd = 1.5, lty = 2)
abline(h = -cutoff, col = "red", lwd = 1.5)
abline(h = cutoff, col = "red", lwd = 1.5)
par(mfrow = c(1,2))
plot(t.cars.full$residuals)
plot(cars.full$residuals)

#calculate MSE
sigma(t.cars.full)

#plot the model, residual plots...
par(mfrow=c(2,2))
plot(t.cars.full)
par(mfrow=c(1,1))


par(mfrow=c(1,2))
qqnorm(t.cars.full.res, main = "t")
qqnorm(cars.full.res)


######### QUESTION 3 ######### 
#variable selection

summary(cars.full)
cars.reduced <- lm(co2 ~ euro_standard + transmission_type + engine_capacity +
                    fuel_type + combined_metric + nox_emissions + extra_urban_metric, 
                    data = mydata)
summary(cars.reduced)
sigma(cars.reduced)

anova(cars.full, cars.reduced)
#can remove: co_emissions, urban_metric, noise_level - that way all the variables are significant

#check whether new model is ok
par(mfrow=c(2,2))
plot(cars.reduced)
par(mfrow=c(1,1))

#question 4

#Can we reomove euro_standard?
anova(cars.reduced) #seems like no
summary(cars.reduced) #euro_standard4 isn't signif diff compared to euro_standard3, 5&6 are different

boxplot(co2~euro_standard, col = c("red", "blue", "yellow", "green"), pch=19, 
        main = "Boxplots of emissions of each engine standard")
#question 5

#confidence interval for beta_1, which is euro_standard4 in our case
confint(cars.reduced,"euro_standard4", level = 0.95)

#question 6

#predict the interval for given variable
x0 <- data.frame(euro_standard = as.factor(4), transmission_type = as.factor("Manual"),
                 engine_capacity = 2196, 
                 fuel_type = as.factor("Petrol"), urban_metric = 9.2, extra_urban_metric = 5.6, combined_metric = 6.9,
                 noise_level = 72, co_emissions = 273.5, nox_emissions = 43)


pred <- predict(cars.reduced, x0, interval = "prediction", level = 0.99)
pred
