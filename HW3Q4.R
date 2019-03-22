#########################################################
#############   BTRY/STSCI 4520   #######################
#############      Homework 3     #######################
############# Due: March 21, 2018 #######################
#########################################################


# Instructions: save this file in the format HW3Q4.R. 
# Complete each question using code below the question number.
# You need only upload this file to CMS. 

# Note, we assume your working directory contains any files
# that accompany this one. 

# Further note: 10% will be deducted if your file produces
# an error when run. If your code produces an error and you 
# cannot find it, comment out that portion of the code and we
# will give partial credit for it. 

# Do not use either the function set.seed() or rm(list=ls())
# in your code. 



##################################
# Question 4: Bootstrap Cautions #
##################################

# While the bootstrap is often regarded as better than pretending
# that your data is normal, it is not perfect. The theory that
# says it works is asymptotic -- it works as the number of
# samples becomes large. 

# Here we will examine an example of this that Tim 
# Hesterberg (senior statistician at Google) presented when 
# he visited Cornell in October 2018. 

# The problem is to bootstrap exponentially distributed data. Where
# how you bootstrap can make a difference. For this question
# we will generate exponentially distributed random variables using
# rexp(n)  where the mean of these should be 1. 

# Throughout the following use nboot = nsim = 1000.

# a) Write a function to construct bootstrap 95% confidencence intervals
# for the mean of a data vector X using normal and percentile standard 
# errors. 

q4afunction = function(X,nboot){
  n = length(X)
  mean.boot = rep(0,nboot)
  for(i in 1:nboot){
    mean.boot[i] = mean(X[sample(n,replace=TRUE)])
  }
  se.boot = sd(mean.boot)
  b0.025 = as.numeric(quantile(mean.boot,0.025))
  b0.975 = as.numeric(quantile(mean.boot,0.975))
  normal.CI = mean(X) + c(-1,1)*qnorm(0.975)*se.boot
  percentile.CI = 2*mean(X)-c(b0.975, b0.025)
	return( list(normal.CI = normal.CI, percentile.CI = percentile.CI))
}



# b) Write a function to conduct a simulation in which you generate
# data set X using rexp(n) each of nsim times and report the percentage
# of times that the mean of 1 is contained in the confidence intervals
# calculated from your function in part a. 

q4bfunction = function(n,nsim,nboot){
  normal.boot = matrix(rep(0,2*nsim),nsim,2)
  per.boot = matrix(rep(0,2*nsim),nsim,2)
  norm.check = rep(0,nsim)
  per.check = rep(0,nsim)
  for(i in 1:nsim){
    X = rexp(n)
    re = q4afunction(X,nboot)
    normal.boot[i,] = re$normal.CI
    per.boot[i,] = re$percentile.CI
    norm.check[i] = (normal.boot[i,1]<=1 &normal.boot[i,2]>=1)
    per.check[i] = (per.boot[i,1]<=1 &per.boot[i,2]>=1)
  }
  return(list( normal.coverage =mean(norm.check), percentile.coverage = mean(per.check)) )
}


# Run this function using n = 10, 500, 100. Do the bootstrap intervals
# cover the true value of 1 95% of the time?
nsim = 1000
nboot = 1000
b1 = q4bfunction(10,nsim,nboot)
b2 = q4bfunction(100,nsim,nboot)
b3 = q4bfunction(500,nsim,nboot)
## I used nsim=1000, and nboot=1000. With the increase of n, the percentage increased, it is more likely to cover the true value of 1 95% of time
## when n=10 only 85% and 83.6% of time covers the true value using normal and percentile confidence interval, which means we were overconfident about how well we were estimating the mean. 
## but when n=500 94.5% and 94.4% of time covers the true value using normal and percentile confidence interval



# c) A suggested alternative is to bootstrap t-statistics. That is,
# if our data is X, and we bootstrap this to get a data set Z, we 
# then form the statistic
#
#     T  =  (mean(Z) - mean(X))/sd(Z)
#
# we then obtain T(0.025) and T(0.975)  (the 2.5% and 97.5% quantiles)
# from the bootstrapped T-statistics and use the following as our
# confidence interval:
#
#    [mean(X) - T(0.957)*sd(X),  mean(X) + T(0.025)*sd(X)]  
#
# Write a function to caculate t-bootstrap confidence intervals for
# a data set X

q4cfunction = function(X,nboot){
  n = length(X)
  mean.boot = rep(0,nboot)
  t.boot = rep(0,nboot)
  for(i in 1:nboot){
    X.boot = X[sample(n,replace=TRUE)]
    mean.boot[i] = mean(X.boot)
    t.boot[i] = (mean.boot[i]-mean(X))/sd(X.boot)
  }
se = sd(mean.boot)
t.upper = as.numeric(quantile(t.boot,0.975))
t.lower = as.numeric(quantile(t.boot,0.025))
interval = c(mean(X)-t.upper*sd(X),mean(X)-t.lower*sd(X))
return(interval)
}

# d) Write a function to repeat the simulation from part b. Does the 
# coverage improve?

q4dfunction = function(n,nsim,nboot){
  int.boot = matrix(rep(0,2*nsim),nsim,2)
  check =  rep(0,nsim)
  for(i in 1:nsim){
    X = rexp(n)
    int.boot[i,] = q4cfunction(X,nboot)
    check[i] = (int.boot[i,1]<=1 &int.boot[i,2]>=1)
  }
return(mean(check))
} 

## Answer:
## Converge does improve. We are more confident about how well we were estimating the mean. 
## when n=10 93.9% of time covers the true value using bootstrap t-statistics, it is closer to the value 95%
## when comparing with the result of normal and percentile confidence interval
## and with the increase of n, it is more likely close to the value95% 
## wg, n=500, the result is 95.6%, which is really close to 95% confidence level.
d1 = q4dfunction(10,1000,1000)
d2 = q4dfunction(100,1000,1000)
d3 = q4dfunction(500,1000,1000)
