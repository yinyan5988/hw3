#########################################################
#############   BTRY/STSCI 4520   #######################
#############      Homework 3     #######################
############# Due: March 21, 2018 #######################
#########################################################


# Instructions: save this file in the format HW3Q1.R. 
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



#####################################################################
# Question 1:  Critical Values, Random Effects and Multiple Testing #
#####################################################################


# Here we will examine a more careful approach to the multiple
# testing correction Giles presented in lecture. 

# Recall that Giles' collaborator is conducting experiments to try 
# to see differences in flourescense between malaria-infected cells 
# treated with Triacin-C and untreated cells.  To carry this out,
# she uses 4 batches of cells which she divides into 9 different
# treatment groups (3 stages and 3 times, so there is more structure,
# but we'll keep it simple and look at 9). In each treatment she divides
# the cells in two and treats half and leaves the other half untreated
# and then measures the difference in log-fluorescense between the two. 

# We are interested in the best way of detecting whether any of the 9 
# different groups exhibit differences between treatments. What makes
# this challenging is that the batches are common across all groups
# so that we can't treat them as independent. 

# A standard model for this would be to specify the following model 
# for Y_ij, the difference between treatment and control in the i'th
# batch and j'th group

#    Y_ij = d_j + b_i + e_ij

# where 
#
#   - d_j is the difference between treatment and control in the 
# j'th group (e.g. how much difference do we see between treatment
# and control in trophs at 15 minutes). We expect this to be the 
# same every time we simulate the experiment. For our null hypothesis
# we will assume that the d_j = 0. 
#
#   - b_i is the effect for whatever might be different in the i'th
# batch; maybe the microscope has a bit of grime on it that reduces
# the fluorescence readings, or the treatment is slightly less well
# spread out, or something. These effects will be different every 
# time we simulate the eximeriment and we assume these can be generated
# from  N(0,sigma^2_b)
#
#  - e_ij is the idiosyncratic errors (maybe the placement of the 
# unique to each observation. They will also be different between
# every simulation and can be generated from N(0,sigma^2_e)

# Here we will try to examine the best (ie, most powerful) way
# to control the family-wise error rate. 


# A) Write a function that simulates the group-wise averages from 
# one such experiment. It should take arguments 
#
#   dvec -- a vector of differences for each group
#   nbatch -- the number of batches used
#   sigma2b -- between batch variance
#   sigma2e -- between experiment variance
#
# and return a vector of length (dvec) giving the average difference
# in each group. 
#
# You can use the fact that if X_1,...,X_n are N(0,s), then bar(X)
# is N(0,s/n) -- our solutions will do so. 

Q1afunction = function(dvec=rep(0,9),nbatch=4,sigma2b,sigma2e)
{
 Y = rep(0,length(dvec))
 b = rnorm(1,0,sd = sqrt(sigma2b/nbatch))
 for(i in 1:length(dvec)){
   Y[i] = dvec[i]+rnorm(1,0,sd = sqrt(sigma2e/nbatch))+b
 }
 return(Y)
}


# B) Use this to simulate a critical value to use to control the
# family wise error rate. To do this, 
#
#   1. For each of nsim simulations, generate data with dvec being 0. 
#   2. In each simulation, for each group form the normalized 
# statistic by dividing by  sqrt((sigma2b + sigma2e)/nbatch) -- the 
# standard deviation of the average.
#   3. Calculate the maximum (over groups) of these statistics
#   4. Find the 95% quantile (over simulations) of these maxima. 
#
# Code this up in the following function that returns the critical
# value


Q1bfunction = function(nsim=1000,nbatch=4,sigma2b,sigma2e,ngroup=9)
{
  re = matrix(rep(0,nsim*ngroup),nrow=nsim)
  for(i in 1:nsim){
    re[i,] = Q1afunction(dvec=rep(0,ngroup),nbatch=nbatch,sigma2b,sigma2e)/sqrt((sigma2b + sigma2e)/nbatch) 
  }
  re.max = apply(re,MARGIN=1,max)
  re.max.crit = quantile(re.max,.95)
  return(as.numeric(re.max.crit))
}

# Carry this out with nsim = 1000 simulations using sigma2b = 0.5 and 
# sigma2e = 1 and nbatch=4 (approximately the result in Giles' data) 
# report the critical value that you obtain below. 

max.crit = Q1bfunction(1000,4,0.5,1)
max.crit
## the critical value is about 2.47 when I run this function
# How does this compare to a one-sided normal critical value after
# employing a Bonferroni correction for making 9 tests? Report this
# value in

alpha_bonf = 0.05/9
bonf.crit =  qnorm(1-alpha_bonf)
bonf.crit
##Critical value getting from Bonferroni corresction is very close to using maximum of statistics which is 2.54


# C) Was using the same batch across all experiments a good idea?
# Let's look at power. Let's suppose that d_j = 2 for each effect. 
# Write a function that evaluates the probability that each group
# is correctly rejected using both the critical values from 
# the simulation in Part B (a vector of 9 probabilities in 
# max.percent, and the Bonferroni correction (also a vector of
# 9 probabilities given in bonf percent)

Q1cfunction = function(dvec = rep(2,9),nsim,nbatch,sigma2b,sigma2e)
{
  re = matrix(rep(0,nsim*length(dvec)),nrow=nsim)
  max.crit = Q1bfunction(nsim,nbatch,sigma2b,sigma2e,ngroup = length(dvec))
  bonf.crit =  qnorm(1-(0.05/length(dvec)))  
  for(i in 1:nsim){
    re[i,] = Q1afunction(dvec=dvec,nbatch=nbatch,sigma2b,sigma2e)/sqrt((sigma2b + sigma2e)/nbatch) 
  }
  max.percent = apply(re,2,function(x) mean(x>max.crit))
  bonf.percent = apply(re,2,function(x) mean(x>bonf.crit))
  return(list( max.percent =max.percent ,  bonf.percent = bonf.percent ) )
}

# Run this function using our settings above. 
pro.rej = Q1cfunction(dvec = rep(2,9),nsim=1000,nbatch=4,0.5,1)
pro.rej

# Run it again using sigma2b = 0 and sigma2e = 1.5
pro.rej1 = Q1cfunction(dvec = rep(2,9),nsim=1000,nbatch=4,0,1.5)
pro.rej1


# Comment on which experimental protocol gives you a better chance
# of detecting real effects.  Is Bonferroni less conservative 
# in one protocol than the other?


##Answer: When using sigma2b = 0.5 and sigma2e = 1 gives a better chance of detecting real effects. That protocol
## has smaller variance in error than the seond one. Although the variance of random effect+variance of error is same
## The result of Bonferroni correction in two protocol it very similar, but the second one is less convervative since
## it reject more test than the first one