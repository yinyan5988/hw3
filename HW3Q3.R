#########################################################
#############   BTRY/STSCI 4520   #######################
#############      Homework 3     #######################
############# Due: March 21, 2018 #######################
#########################################################


# Instructions: save this file in the format HW3Q3.R. 
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



##############################
# Question 3: Knockoff Tests #
##############################

# We've seen that you generally can't conduct permutation 
# tests for individual covariates in multiple regression. 
# In one case, however, this is possible. The data in 

peanuts = read.table('peanuts.txt',sep=',',head=TRUE)

# contains data on the yield of peanuts (in pounds) from
# each plant. It's height and whether it had been treated
# with a control,  fast-release fertilizer or slow-release
# fertilizer. 

# Typically, we'd be interested in the effect of treatment,
# but here we will test whether there is a relationship 
# between height and yield, while controlling for 
# fertilizer. 

# The key is that we can avoid breaking the relationship
# between fertilizer and height by permuting height 
# WITHIN a level of fertilizer. This will break any 
# relationship between height and yield while keeping 
# the relationship between the covariates. 

# A) Write a function to conduct a permutation test for 
# the coefficient of height while controlling for 
# fertilizer using the framework above. Report 
# the observed t-statistic, the p-value and the 
# critical value for your test.  We will use
# the estimated coefficient of Height as a statistic.



PartialPermutation = function(data,nperm)
{
  coef.vec = rep(0,nperm)
  obs.coef = summary(lm(data[,1]~data[,2]+data[,3]))$coefficients[2,1]
  for(i in 1:nperm){
    c = sample(10)
    s = sample(11:20)
    f = sample(21:30)
    ord = c(c,s,f)
    re = summary(lm(data[,1]~data[ord,2]+data[,3]))
    coef.vec[i] = re$coefficients[2,1]
  }
  crit.value = as.numeric(quantile(coef.vec,0.975))
  pvalue = mean(abs(coef.vec)>obs.coef)
  return( list( obs.coef = obs.coef,
                pvalue = pvalue,
                crit.value = crit.value) )
}

# Test this for yourself using nperm = 1000. 
re = PartialPermutation(peanuts,nperm=1000)
## p value is 0 and critical value is about 0.0218

# This is a special case of the idea of knock-offs
# suggested in 2014 in 
#
#  https://arxiv.org/abs/1404.5609
#
# which were designed for False Discovery Rates. The
# basic idea is that you need to create a new version 
# of x_i that has the same relationship to the other
# covariates, but no relationship to y.  In knockoffs,
# you generate a new  x*_i|x_{-i}   (meaning x_i given 
# all the other covariates). This x*_i has no more 
# information about y than x_{-i} does and so we can
# compare how much it appears to predict to how much
# the real x_i does.  




# B) We can apply the same logic to create a boostrap 
# that respects the balance of the experiment assigned
# to fertilizer. 
#
# Write a function to produce a bootstrap confidence
# interval for the slope of yeild on height while 
# accounting for fertilizer. You should bootstrap 
# values WITHIN each level of fertilizer.
#
# Your function should take in a data set and output
# the estimated slope and both normal-theory and 
# percentile bootstrap confidence intervals. 

PartialBoot = function(data,nboot)
{
  obs.coef = summary(lm(data[,1]~data[,2]+data[,3]))$coefficients[2,1]
  c.boot = rep(0,nboot)
  for(i in 1:nboot){
    c = sample(10,replace = TRUE)
    s = sample(11:20,replace = TRUE)
    f = sample(21:30,replace = TRUE)
    ord = c(c,s,f)
    re = summary(lm(data[ord,1]~data[ord,2]+data[,3]))
    c.boot[i] = re$coefficients[2,1]
  }
  p.upper <-as.numeric(quantile(c.boot,1-0.05/2))
  p.lower <- as.numeric(quantile(c.boot,0.05/2))
  c.se = sd(c.boot)
   return( list( obs.coef = obs.coef ,
                 normal.interval = obs.coef + c(-1,1)*qnorm(0.975)*c.se,
				 percentile.interval = c(2*obs.coef-p.upper, 2*obs.coef-p.lower)))
}


# Is there a substantial difference between the 
# bootstrap intervals? 

## Answer: The two confidence intervals are really similar has no substantial difference. 

