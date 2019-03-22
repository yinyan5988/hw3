#########################################################
#############   BTRY/STSCI 4520   #######################
#############      Homework 3     #######################
############# Due: March 21, 2018 #######################
#########################################################


# Instructions: save this file in the format HW3Q2.R. 
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



################################
# Question 2: Multiple Testing #
################################

# This question is associated with recent controversies
# over poor statistical practice, particularly in 
# psychology and nutrition. See a piece from February 25, 2018
# in buzzfeed at:
#
# https://www.buzzfeed.com/stephaniemlee/brian-wansink-cornell-p-hacking?utm_term=.nxoJ87k0pG#.ilJpwEGJWX
#
# and a somewhat less aggressive view of the report at
# 
# http://andrewgelman.com/2018/02/27/no-researchers-not-typically-set-prove-specific-hypothesis-study-begins/
# 
# The article describes e-mails that discuss routinely scanning
# many possible hypotheses looking for an effect. 
#
# Here we'll examine one particular case study of an attempt to 
# find significance of a relationship between eating and TV watching 
# in which 400 mediation analyses were trialled. 
#
# In mediation analysis, one seeks to find an effect that is partly 
# obscured by a different covariate. We'll be a bit simpler, and 
# look at trying to find a relationship between y and x while
# controling for another possible covariate z, where there are 400 
# possible z's to look at. 

# We'll consider the following setup. 10 subjects are assigned to 
# watch TV while eating pizza and 10 subjects assigned to read. We
# record the total caloric intake of each subject along with 
# 400 survey variables about their demographics etc. 

# Here y is caloric intake, x is a binary indicator of whether
# the subject watched TV and we have 400 z_i's.  We'll simulate data
# where y has no relationship to x, but might be related to z_i. That
# is we set

x = rep( c(0,1), 1, each = 10)

# and generate y as 20 standard normals. We'll generate each z 
# column by
# 
# z = rnorm(1,sd=0.5)*y + rnorm(20)
#
# The first random number here is a 'coefficient' for a relationship 
# with y, except we will end up treating y as a response. 
 
# a) Write a function to generate data as above, run a linear 
# regression y ~ x + z_i  and report the p-value of the 
# coefficient for x for i = 1,...,400. 



Q2a.function = function(){
y = rnorm(20,0,1)
beta = rnorm(400,0,0.5)
e = matrix(rnorm(400*20),20,400)
p.value = rep(0,400)
for(i in 1:400){
  z=beta[i]* y + e[,i]
  p.value[i] = summary(lm(y~x+z))$coefficients[2,4]
}
return(p.value)
}


# b) Conduct 100 simulations to determine 
# 
#  i)  How often at least one significant effect is found for x
#  when conducting all tests at the 0.05 level. 
# 
#  ii) How often at least one significant effect is found for x
#  when making a Bonferroni correction.
#
#  iii) The average size of the list of confounders that make 
#  x significant if you apply a false discovery rate procedure 
#  and control the FDR at 0.1.  To do so, obtain the p-values for 
#  x and use this to produce q-values; you can then select the 
#  confounders with q-value less than Q. 
#
# Write your procedure in the following function

Q2simfunc = function(nsim)
{
  p.value = matrix(rep(0,nsim*400),nsim,400)
  q.value = matrix(rep(0,nsim*400),nsim,400)
  p.sort = matrix(rep(0,nsim*400),nsim,400)
for(i in 1:nsim){
  p.value[i,] = Q2a.function()
  p.sort[i,] = sort(p.value[i,])
}
sig.effect = mean(apply(p.value,1,function(x) any(x<0.05)))
bon.effect = mean(apply(p.value,1,function(x) any(x<0.05/400)))
m = 1:400
list.len = rep(0,nsim)
for(i in 1:nsim){
  q.value[i,] = p.sort[i,]*400/m
}
list.len = apply(q.value,1,function(x) sum(x<0.1))
fdr.size = mean(list.len)
  return( list( one.significant = sig.effect,
                bonferroni.significant = bon.effect,
				fdr.size = fdr.size ) ) 
}
 


Q2simfunc = function(nsim)
{
  p.value = matrix(rep(0,nsim*400),nsim,400)
  q.value = matrix(rep(0,nsim*400),nsim,400)
  p.sort = matrix(rep(0,nsim*400),nsim,400)
  for(i in 1:nsim){
    p.value[i,] = Q2a.function()
    p.sort[i,] = sort(p.value[i,])
  }
  sig.effect = mean(apply(p.value,1,function(x) any(x<0.05)))
  bon.effect = mean(apply(p.value,1,function(x) any(x<0.05/400)))
  list.len = rep(0,nsim)
  for(i in 1:nsim){
    for(j in 1:400){
      q.value[i,j] = 400*p.sort[i,j]/sum(p.sort[i,]<p.sort[i,j])
    }
  }
  list.len = apply(q.value,1,function(x) sum(x<0.1))
  fdr.size = mean(list.len)
  return( list( one.significant = sig.effect,
                bonferroni.significant = bon.effect,
                fdr.size = fdr.size ) ) 
}
 
# c) What do you believe is the correct procedure to carry 
# out in this case?

##with 400 z doing regression, we have more chance to find something significant. 
##even if x doesn't play a role in generating y. 
##using 0.05 as significant level is not proper, it would increase fdr.
## using bonferroni correction is a too conservative way for multi test. In this question
## it performs good because null pypothesis is the right thing. In the simulation, only 18 times at least one significant effect is found
## But for fdr, the average confounders with q-value less than Q is 12.434, it is small comparing with 400 tests
##It is better to both control q value and also p value. 


