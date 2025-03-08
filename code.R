setwd("C:/Users/annaz/Documents/HSE-2/Experimetrics/Homeworks")
data <- read.csv("fs_sim_ZHURBA.csv")



##Firstly, let's assume the model without gender and age
## Just simple model with 2 types of people in data:
LogLK <- function(par, data){
  d <- data$d
  bottom_1 <- data$bottom_1
  bottom_2 <- data$bottom_2
  top_1 <- data$top_1
  top_2 <- data$top_2
  
  alpha1 <- par[1]
  beta1 <- par[2]
  alpha2 <- par[3]
  beta2 <- par[4]
  lambda <- par[5]
  theta1 <- par[6]
  
  #Probability of being type 1
  p1 <- inv.logit(theta1)
  p2 <- 1-p1
  
  #Determine whether advantageous or disadvantageous inequality
  a_top <- ifelse(top_2>top_1,1,0)
  a_bottom <- ifelse(bottom_2>bottom_1,1,0)
  
  #Utility of first type
  u_bottom1 <- bottom_1 - a_bottom*alpha1*(bottom_2-bottom_1) - (1-a_bottom)*beta1*(bottom_1-bottom_2)
  u_top1 <- top_1 - a_top*alpha1*(top_2-top_1) - (1-a_top)*beta1*(top_1-top_2)
  diff_ia1 <- u_bottom1 - u_top1
  
  #Utility of second type
  u_bottom2 <- bottom_1 - a_bottom*alpha2*(bottom_2-bottom_1) - (1-a_bottom)*beta2*(bottom_1-bottom_2)
  u_top2 <- top_1 - a_top*alpha2*(top_2-top_1) - (1-a_top)*beta2*(top_1-top_2)
  diff_ia2 <- u_bottom2 - u_top2
  
  #Likelihood conditional on type
  t1 <- d*(inv.logit(lambda*diff_ia1))+(1-d)*(inv.logit(-lambda*diff_ia1))
  t2 <- d*(inv.logit(lambda*diff_ia2))+(1-d)*(inv.logit(-lambda*diff_ia2))
  
  likelihoods <- p1*t1+(1-p1)*t2
  ll <- sum(log(likelihoods))
  return(-ll)
}

#Maximize log-likelihood
mleK <- optim(c(0.1,0,0.25,0.3,1,0),LogLK,data=data, hessian=TRUE, method="BFGS")
mleK$par
mleK$value

p1 <- inv.logit(mleK$par[6])
p2 <- 1-p1
p1
p2
p1+p2

stdK <- sqrt(diag(solve(mleK$hessian)))
stdK

## Check statistic significance of our estimates
#Two-tailed test for significance of estimates (H0: estimate=0)

ts_1 <- mean(mleK2$par[1]/stdK[1])
2*min((pnorm(ts_1)), (pnorm(ts_1, lower.tail = FALSE))) 

ts_2 <- mean(mleK$par[2]/stdK[2])
2*min((pnorm(ts_2)), (pnorm(ts_2, lower.tail = FALSE))) 

ts_3 <- mean(mleK$par[3]/stdK[3])
2*min((pnorm(ts_3)), (pnorm(ts_3, lower.tail = FALSE))) 

ts_4 <- mean(mleK$par[4]/stdK[4])
2*min((pnorm(ts_4)), (pnorm(ts_4, lower.tail = FALSE))) 

ts_5 <- mean(mleK$par[5]/stdK[5])
2*min((pnorm(ts_5)), (pnorm(ts_5, lower.tail = FALSE))) 

ts_6 <- mean(mleK$par[6]/stdK[6])
2*min((pnorm(ts_6)), (pnorm(ts_6, lower.tail = FALSE))) 


#BOOTSTRAPPED STANDARD ERRORS WITH CLUSTERING
set.seed(1)
nboots <- 100
id <- sort(unique(data$subject))
bootpars <- NULL

for(i in 1:nboots){
  bootid <- sample(id,length(id),replace=TRUE)
  selectrows <- unlist(lapply(seq(1:length(bootid)), function(n) which(data$subject == bootid[n])))
  bootdata <- data[selectrows,]
  bootpars <- rbind(bootpars,optim(c(0.1,0,0.25,0.3,1,0),LogLK, data=bootdata, hessian=FALSE)$par)
}

alpha1_err = sqrt(sum((bootpars[,1]-mleK$par[1])^2)/nrow(bootpars))
beta1_err = sqrt(sum((bootpars[,2]-mleK$par[2])^2)/nrow(bootpars))
alpha2_err = sqrt(sum((bootpars[,3]-mleK$par[3])^2)/nrow(bootpars))
beta2_err = sqrt(sum((bootpars[,4]-mleK$par[4])^2)/nrow(bootpars))
lambda_err = sqrt(sum((bootpars[,5]-mleK$par[5])^2)/nrow(bootpars))
theta1_err = sqrt(sum((bootpars[,6]-mleK$par[6])^2)/nrow(bootpars))

alpha1_err
beta1_err
alpha2_err
beta2_err
lambda_err
theta1_err


## Check theta1
library(msm)
estmean <- (mleK$par[6])
estvar <- (stdK[6])^2  
formula <- ~ exp(x1) / (1 + exp(x1))
## for theta1
std_err <- deltamethod(formula, estmean, estvar)
print(std_err)
#Check theta1 is significantly greater than 0
1-pnorm(inv.logit(estmean)/std_err)


## Check statistic significance of our estimates with bootstrap errors
#Two-tailed test for significance of estimates (H0: estimate=0)
ts_1 <- mean(mleK2$par[1]/alpha1_err)
2*min((pnorm(ts_1)), (pnorm(ts_1, lower.tail = FALSE))) 

ts_2 <- mean(mleK$par[2]/beta1_err)
2*min((pnorm(ts_2)), (pnorm(ts_2, lower.tail = FALSE))) 

ts_3 <- mean(mleK$par[3]/alpha2_err)
2*min((pnorm(ts_3)), (pnorm(ts_3, lower.tail = FALSE))) 

ts_4 <- mean(mleK$par[4]/beta2_err)
2*min((pnorm(ts_4)), (pnorm(ts_4, lower.tail = FALSE))) 

ts_5 <- mean(mleK$par[5]/lambda_err)
2*min((pnorm(ts_5)), (pnorm(ts_5, lower.tail = FALSE))) 

ts_6 <- mean(mleK$par[6]/theta1_err)
2*min((pnorm(ts_6)), (pnorm(ts_6, lower.tail = FALSE))) 







## Now let's estimate models with gender and age
###ESTIMATE ALPHA AND BETA FOR SINGLE TYPE###
library(boot)
LogLK1 <- function(par, data){
  d <- data$d
  bottom_1 <- data$bottom_1
  bottom_2 <- data$bottom_2
  top_1 <- data$top_1
  top_2 <- data$top_2
  age <- data$age
  male <- data$male
  male_p <- data$male_p
  
  alpha <- par[1]
  beta <- par[2]
  lambda <- par[3]
  b0 <- par[4]
  b1 <- par[5]
  b2 <- par[6]
  
  #Determine whether advantageous or disadvantageous inequality
  a_top <- ifelse(top_2>top_1,1,0)
  a_bottom <- ifelse(bottom_2>bottom_1,1,0)
  gender <- ifelse(male == male_p,1,0)
  age <- ifelse(age>=35,1,0)
  
  #Compute utility difference between top and bottom
  u_bottom <- (bottom_1 - a_bottom*alpha*(bottom_2-bottom_1) - ((1-a_bottom)*beta*(bottom_1-bottom_2)))^(b0+b1*gender+b2*age)
  u_top <- (top_1 - (b0+b1*gender+b2*age)*a_top*alpha*(top_2-top_1) - (b0+b1*gender+b2*age)*(1-a_top)*beta*(top_1-top_2))^(b0+b1*gender+b2*age)
  diff_ia <- u_bottom - u_top
  
  likelihoods <- d*(inv.logit(lambda*diff_ia))+(1-d)*(inv.logit(-lambda*diff_ia))
  ll <- sum(log(likelihoods))
  return(-ll)
}
#Maximize log-likelihood
mleK1 <- optim(c(0,0,1,0.5,0.5,0.5),LogLK1,data=data, hessian=TRUE)
mleK1
stdK1 <- sqrt(diag(solve(mleK111$hessian)))
stdK1



###ESTIMATE ALPHA AND BETA IF TWO TYPES###
library(boot)
LogLK2 <- function(par, data){
  d <- data$d
  bottom_1 <- data$bottom_1
  bottom_2 <- data$bottom_2
  top_1 <- data$top_1
  top_2 <- data$top_2
  age <- data$age
  male <- data$male
  male_p <- data$male_p
  
  alpha1 <- par[1]
  beta1 <- par[2]
  alpha2 <- par[3]
  beta2 <- par[4]
  
  lambda <- par[5]
  theta1 <- par[6]
  b0 <- par[7]
  b1<-par[8]
  b2<-par[9]
  
  a0 <- par[10]
  a1<-par[11]
  a2<-par[12]
  
  #Probability of being type 1
  p1 <- inv.logit(theta1)
  #Probability of being type 2
  p2 <- 1-p1
  
  #Determine whether advantageous or disadvantageous inequality
  a_top <- ifelse(top_2>top_1,1,0)
  a_bottom <- ifelse(bottom_2>bottom_1,1,0)
  gender <- ifelse(male == male_p,1,0)
  age <- ifelse(age>30,1,0)
  
  #Compute utility for type 1
  u_bottom1 <- (bottom_1 - a_bottom*alpha1*(bottom_2-bottom_1) - (1-a_bottom)*beta1*(bottom_1-bottom_2))^(a0+a1*gender+a2*age)
  u_top1 <- (top_1 - a_top*alpha1*(top_2-top_1) - (1-a_top)*beta1*(top_1-top_2))^(a0+a1*gender+a2*age)
  diff_ia1 <- u_bottom1 - u_top1
  
  #Compute utility for type 2
  u_bottom2 <- (bottom_1 - a_bottom*alpha2*(bottom_2-bottom_1) - (1-a_bottom)*beta2*(bottom_1-bottom_2))^(b0+b1*gender+b2*age)
  u_top2 <- (top_1 - a_top*alpha2*(top_2-top_1) - (1-a_top)*beta2*(top_1-top_2))^(b0+b1*gender+b2*age)
  diff_ia2 <- u_bottom2 - u_top2
  
  
  #Likelihood conditional on type
  t1 <- d*(inv.logit(lambda*diff_ia1))+(1-d)*(inv.logit(-lambda*diff_ia1))
  t2 <- d*(inv.logit(lambda*diff_ia2))+(1-d)*(inv.logit(-lambda*diff_ia2))
  
  likelihoods <- p1*t1+(1-p1)*t2
  ll <- sum(log(likelihoods))
  return(-ll)
}
#Maximize log-likelihood
mleK2 <- optim(c(0.1,0,0.25,0.3,1,0,0.4,0.25,0.8,0,0,0),LogLK2,data=data, hessian=TRUE)
#, method="BFGS"
mleK2

## Look at probabilities (shares) of each type
p1 <- inv.logit(mleK2$par[6])
p2 <- 1-p1
p1
p2
p1+p2

## Simple standart errors:
stdK2 <- sqrt(diag(solve(mleK2$hessian)))
stdK2

## Check theta1
library(msm)
estmean <- (mleK2$par[6])
estvar <- (stdK2[6])^2  
formula <- ~ exp(x1) / (1 + exp(x1))
## for theta1
std_err <- deltamethod(formula, estmean, estvar)
print(std_err)
#Check theta1 is significantly greater than 0
1-pnorm(inv.logit(estmean)/std_err)

## We can see that stamdart errors can't be camputed, thus let's use bootstrap errors
#BOOTSTRAPPED STANDARD ERRORS WITH CLUSTERING (take into account the dependence of observations)
set.seed(1)
nboots <- 100
ids <- sort(unique(data$subject))
bootpars1 <- NULL

for(i in 1:nboots){
  bootid <- sample(ids,length(ids),replace=TRUE)
  selectrows <- unlist(lapply(seq(1:length(bootid)), function(n) which(data$subject == bootid[n])))
  bootdata <- data[selectrows,]
  bootpars1 <- rbind(bootpars1,optim(c(0.1,0,0.25,0.3,1,0,0.4,0.25,0.8,0,0,0),LogLK2, data=bootdata, hessian=FALSE)$par)
}


# Estimates errors with bootstrap:
alpha1_err = sqrt(sum((bootpars1[,1]-mleK2$par[1])^2)/nrow(bootpars1))
beta1_err = sqrt(sum((bootpars1[,2]-mleK2$par[2])^2)/nrow(bootpars1))
alpha2_err = sqrt(sum((bootpars1[,3]-mleK2$par[3])^2)/nrow(bootpars1))
beta2_err = sqrt(sum((bootpars1[,4]-mleK2$par[4])^2)/nrow(bootpars1))
lambda_err = sqrt(sum((bootpars1[,5]-mleK2$par[5])^2)/nrow(bootpars1))
theta1_err = sqrt(sum((bootpars1[,6]-mleK2$par[6])^2)/nrow(bootpars1))

b0_err = sqrt(sum((bootpars1[,7]-mleK2$par[7])^2)/nrow(bootpars1))
b1_err = sqrt(sum((bootpars1[,8]-mleK2$par[8])^2)/nrow(bootpars1))
b2_err = sqrt(sum((bootpars1[,9]-mleK2$par[9])^2)/nrow(bootpars1))
a0_err = sqrt(sum((bootpars1[,10]-mleK2$par[10])^2)/nrow(bootpars1))
a1_err = sqrt(sum((bootpars1[,11]-mleK2$par[11])^2)/nrow(bootpars1))
a2_err = sqrt(sum((bootpars1[,12]-mleK2$par[12])^2)/nrow(bootpars1))

alpha1_err
beta1_err
alpha2_err
beta2_err
lambda_err
theta1_err
b0_err
b1_err
b2_err
a0_err
a1_err
a2_err


## Check statistic significance of our estimates
#Two-tailed test for significance of estimates (H0: estimate=0)
ts_1 <- mean(mleK2$par[1]/alpha1_err)
2*min((pnorm(ts_1)), (pnorm(ts_1, lower.tail = FALSE))) 

ts_2 <- mean(mleK$par[2]/beta1_err)
2*min((pnorm(ts_2)), (pnorm(ts_2, lower.tail = FALSE))) 

ts_3 <- mean(mleK$par[3]/alpha2_err)
2*min((pnorm(ts_3)), (pnorm(ts_3, lower.tail = FALSE))) 

ts_4 <- mean(mleK$par[4]/beta2_err)
2*min((pnorm(ts_4)), (pnorm(ts_4, lower.tail = FALSE))) 

ts_5 <- mean(mleK$par[5]/lambda_err)
2*min((pnorm(ts_5)), (pnorm(ts_5, lower.tail = FALSE))) 

ts_6 <- mean(mleK$par[6]/theta1_err)
2*min((pnorm(ts_6)), (pnorm(ts_6, lower.tail = FALSE))) 


ts_7 <- mean(mleK2$par[7]/b0_err)
2*min((pnorm(ts_1)), (pnorm(ts_7, lower.tail = FALSE))) 

ts_8 <- mean(mleK2$par[8]/b1_err)
2*min((pnorm(ts_8)), (pnorm(ts_8, lower.tail = FALSE))) 

ts_9 <- mean(mleK2$par[9]/b2_err)
2*min((pnorm(ts_9)), (pnorm(ts_9, lower.tail = FALSE))) 

ts_10 <- mean(mleK2$par[10]/a0_err)
2*min((pnorm(ts_10)), (pnorm(ts_10, lower.tail = FALSE))) 

ts_11 <- mean(mleK2$par[11]/a1_err)
2*min((pnorm(ts_11)), (pnorm(ts_11, lower.tail = FALSE))) 

ts_12 <- mean(mleK2$par[12]/a2_err)
2*min((pnorm(ts_12)), (pnorm(ts_12, lower.tail = FALSE))) 
##Just two estimates are significant: beta2 (p-value = 0.00003717166) theta1 (p-value = 0.000001288363)


## Now let's compare our models with Vuong and Clarke tests
#### For model without age and gender ####
d <- data$d
bottom_1 <- data$bottom_1
bottom_2 <- data$bottom_2
top_1 <- data$top_1
top_2 <- data$top_2

alpha1 <- mleK$par[1]
beta1 <- mleK$par[2]
alpha2 <- mleK$par[3]
beta2 <- mleK$par[4]
lambda <- mleK$par[5]
theta1 <- mleK$par[6]

#Probability of being type 1
p1 <- inv.logit(theta1)
p2 <- 1-p1

#Determine whether advantageous or disadvantageous inequality
a_top <- ifelse(top_2>top_1,1,0)
a_bottom <- ifelse(bottom_2>bottom_1,1,0)

#Utility of first type
u_bottom1 <- bottom_1 - a_bottom*alpha1*(bottom_2-bottom_1) - (1-a_bottom)*beta1*(bottom_1-bottom_2)
u_top1 <- top_1 - a_top*alpha1*(top_2-top_1) - (1-a_top)*beta1*(top_1-top_2)
diff_ia1 <- u_bottom1 - u_top1

#Utility of second type
u_bottom2 <- bottom_1 - a_bottom*alpha2*(bottom_2-bottom_1) - (1-a_bottom)*beta2*(bottom_1-bottom_2)
u_top2 <- top_1 - a_top*alpha2*(top_2-top_1) - (1-a_top)*beta2*(top_1-top_2)
diff_ia2 <- u_bottom2 - u_top2

#Likelihood conditional on type
t1 <- d*(inv.logit(lambda*diff_ia1))+(1-d)*(inv.logit(-lambda*diff_ia1))
t2 <- d*(inv.logit(lambda*diff_ia2))+(1-d)*(inv.logit(-lambda*diff_ia2))

likelihoods <- p1*t1+(1-p1)*t2
mean(likelihoods)

#### For model with gender and age ####
d <- data$d
bottom_1 <- data$bottom_1
bottom_2 <- data$bottom_2
top_1 <- data$top_1
top_2 <- data$top_2
age <- data$age
male <- data$male
male_p <- data$male_p

alpha1 <- mleK2$par[1]
beta1 <- mleK2$par[2]
alpha2 <- mleK2$par[3]
beta2 <- mleK2$par[4]

lambda <- mleK2$par[5]
theta1 <- mleK2$par[6]
b0 <- mleK2$par[7]
b1<-mleK2$par[8]
b2<-mleK2$par[9]

a0 <- mleK2$par[10]
a1<-mleK2$par[11]
a2<-mleK2$par[12]

#Probability of being type 1
p1 <- inv.logit(theta1)
#Probability of being type 2
p2 <- 1-p1

#Determine whether advantageous or disadvantageous inequality
a_top <- ifelse(top_2>top_1,1,0)
a_bottom <- ifelse(bottom_2>bottom_1,1,0)
gender <- ifelse(male == male_p,1,0)
age <- ifelse(age>30,1,0)

#Compute utility for type 1
u_bottom1 <- (bottom_1 - a_bottom*alpha1*(bottom_2-bottom_1) - (1-a_bottom)*beta1*(bottom_1-bottom_2))^(a0+a1*gender+a2*age)
u_top1 <- (top_1 - a_top*alpha1*(top_2-top_1) - (1-a_top)*beta1*(top_1-top_2))^(a0+a1*gender+a2*age)
diff_ia1 <- u_bottom1 - u_top1

#Compute utility for type 2
u_bottom2 <- (bottom_1 - a_bottom*alpha2*(bottom_2-bottom_1) - (1-a_bottom)*beta2*(bottom_1-bottom_2))^(b0+b1*gender+b2*age)
u_top2 <- (top_1 - a_top*alpha2*(top_2-top_1) - (1-a_top)*beta2*(top_1-top_2))^(b0+b1*gender+b2*age)
diff_ia2 <- u_bottom2 - u_top2


#Likelihood conditional on type
t1 <- d*(inv.logit(lambda*diff_ia1))+(1-d)*(inv.logit(-lambda*diff_ia1))
t2 <- d*(inv.logit(lambda*diff_ia2))+(1-d)*(inv.logit(-lambda*diff_ia2))

likelihoods_g_a <- p1*t1+(1-p1)*t2


## Vuong statistic 
#Calculator numerator of Vuong statistic
mean(likelihoods_g_a)
mean(likelihoods)
v_num <- sum(log(likelihoods/likelihoods_g_a))/sqrt(length(likelihoods))
v_den <- sqrt(sum((log(likelihoods/likelihoods_g_a)^2-(sum(log(likelihoods/likelihoods_g_a))/length(likelihoods))^2))/length(likelihoods))
2*(1-pnorm(v_num/v_den))
## Vuong test showed that One model capture data better than another. Need to compare the value of likelihoods


# Clarke Test - nonparametric test, which compare the medians of average individual likelihoods.
# H0: true median difference is not equal to 0 ( Both models explain the data equally well).
library(DescTools)
median(likelihoods_g_a)
median(likelihoods)
SignTest(likelihoods_g_a,likelihoods)
## Clarke Test showed that one model capture data better than another 

## Let's compare the Log Likelihoods:
mleK2$value
mleK$value

## Model with gender and age is better!






