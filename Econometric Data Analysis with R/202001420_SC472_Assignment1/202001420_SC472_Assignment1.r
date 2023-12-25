library(wooldridge)
library(stargazer)
library(tidyverse)


#C1
load("/Users/divya/Documents/Semester-6/Econometric Data Analysis with R/R data sets for 5e/401k.RData")
k401k_data <- k401k # Define the data set

# (i) Find the average participation rate and the average match rate in the sample of plans.
participation_rate <- mean(k401k_data$prate) # The average of participation rate
participation_rate
match_rate <- mean(k401k_data$mrate) # The average of match rate
match_rate

# (ii) Now, estimate the simple regression equation
regression1 <- lm(prate ~ mrate, data = k401k_data) # Regression prate on mrate
ggplot(k401k_data, aes(x = mrate, y = prate)) + 
  geom_point(color = 'orange') + 
  geom_smooth(method = 'lm') + 
  labs(x = 'mrate', y = 'prate') # A scatterplot and SRF for mrate & prate
summary(regression1) # The report for results (sample size and R-squared)

# (iii) -> We can answer the question without using R.

# (iv)
83.0755 + 5.8611 * (3.5) # It isn't reasonable because of the characteristics of data 'prate'.

# (v) -> IN the document submitted.
SSE_p <- sum(predict(regression1) - mean(k401k_data$prate)) # SSE for prate

# C2
load("/Users/divya/Documents/Semester-6/Econometric Data Analysis with R/R data sets for 5e/ceosal2.Rdata")
ceosal2_data <- ceosal2 # Define the data set

# (i)
mean(ceosal2_data$salary) # The average of salary
mean(ceosal2_data$ceoten) # The average of tenure

# (ii)
first_tenure <- subset(ceosal2_data, ceoten == '0')
length(first_tenure) # Number of first conten
max(ceosal2_data$ceoten) # Max of conten

# (ii)
ceosal2_data %>% filter(ceoten == '0') %>% length # Number of CEOS in the first years
max(ceosal2_data$ceoten) # The longest tenure

# (iii)
regression2 <- lm(log(salary) ~ ceoten, data = ceosal2_data)
ggplot(ceosal2_data, aes(x = ceoten, y = log(salary))) + 
  geom_point(color = 'red') + 
  geom_smooth(method = 'lm') + 
  labs(x = 'ceoten', y = 'log(salary)') # A scatterplot and SRF for ceoten & log(salary)
summary(regression2) # The report for results(sample size and R-squared)

## The (approximate) predicted percentage increase in salary given one more year as a CEO
regression2$coefficients[2] * 100


# C3
load("/Users/divya/Documents/Semester-6/Econometric Data Analysis with R/R data sets for 5e/sleep75.Rdata")
sleep75_data <- sleep75
# (i)
regression3 <- lm(sleep ~ totwrk, data = sleep75_data)
ggplot(sleep75_data, aes(x = totwrk, y = sleep)) + 
  geom_point(color = 'blue') + 
  geom_smooth(method = 'lm') + 
  labs(x = 'totwrk', y = 'sleep') # A scatterplot and SRF for totwrk & sleep
summary(regression3) # The report for results(sample size and R-squared)

# (ii)
# Two more hour 
regression3$coefficients[2] * 2*60

# One more hour 5 days 
regression3$coefficients[2] * 1*60*5

#C4
load("/Users/divya/Documents/Semester-6/Econometric Data Analysis with R/R data sets for 5e/wage2.Rdata")
w<-wage2 #regressand=wage, regressor=IQ

#(i)
mean(w$wage)
mean(w$IQ)
sd(w$IQ)

#(ii)
reg4<-lm(wage~IQ, data=w)
ggplot(w, aes(x=IQ, y=wage))+geom_point(color='red')+geom_smooth(method = 'lm')+labs(x='IQ',y='wage') #A scatterplot and SRF for IQ & lwage
summary(reg4) #The report for results(sample size and R-squared)

##if IQ increases about 15, then wage changes....
reg4$coefficients[2]*15

## the variation in wage explained by IQ is...

SSE_w<- sum(predict(reg4)-mean(w$wage)) # SSE for wage(You can also check the R-squared)
SSE_w

#(iii)
reg4_1<-lm(lwage~IQ, data=w)
ggplot(w, aes(x=IQ, y=lwage))+geom_point(color='brown')+geom_smooth(method = 'lm')+labs(x='IQ',y='log(wage)') #A scatterplot and SRF for IQ & log(wage)
summary(reg4_1) #The report for results(sample size and R-squared)

##the approximate percentage increase in predicted wage if IQ increases by 15 points
reg4_1$coefficients[2]*100*15 


#C8

#(i)
x<-runif(500, min = 0, max = 10) #Uniform distribution(0,10), n=500
head(x)
mean(x)
sd(x)


#(ii)
u<-rnorm(500,0,36) #Normal distribution(0,36), n=500
mean(u) #Why isn't it zero??
sd(u)

#(iii)
y=1+2*x+u
reg8<-lm(y~x)
dat<-data.frame(x,y)
head(dat)

ggplot(dat, aes(x=x, y=y))+geom_point(color='red')+geom_smooth(method = 'lm')+labs(x='x',y='y') #A scatterplot and SRF for x & y
summary(reg8) #The report for results(sample size and R-squared)

#(iv)
sum(reg8$residuals)
sum(reg8$residuals*x)

#(v)
sum(u)
sum(u*x)

#(vi)->Repeat (i)~(iii) with new random variables
x <- runif(500, min = 0, max = 10) # Uniform distribution(0,10), n=500
head(x) 
mean(x)
sd(x)

u <- rnorm(500, 0, 36) # Normal distribution(0,36), n=500
mean(u) 
sd(u)

y <- 1 + 2 * x + u  
reg8_1 <- lm(y ~ x)
dat <- data.frame(x, y)
head(dat)
ggplot(dat, aes(x = x, y = y)) + geom_point(color = 'red') + geom_smooth(method = 'lm') + labs(x = 'x', y = 'y') # A scatterplot and SRF for x & y
summary(reg8_1) # The report for results(sample size and R-squared)

