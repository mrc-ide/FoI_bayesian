#================================================================================================================================#
#=========================== Simple cataltyic model fit to multiple dataset (joint fitting)  (Dixon et al. ) ====================#
#================================================================================================================================# 

#=====================================================#
#                 LIBRARIES                           #
library(dplyr)
library(ggplot2)
library(MASS) # for covariance matrix
library(data.table)


#======================#
#   Load data          #

#SI_data <- read.csv("")  # read in csv 
#View(SI_data)

# subset Ab-EITB datasets #
data_joint <- SI_data[SI_data$Dataset_name == ("Garcia et al. 2003") | 
                        SI_data$Dataset_name == ("Jayashi et al. 2012") | 
                        SI_data$Dataset_name == ("Lescano et al. 2007") | 
                        SI_data$Dataset_name == ("Taico et al. 2003") | 
                        SI_data$Dataset_name == ("Sarti et al. 2000"),]

# note code used below for example fitting to 5 datasets using Ab-EITB (will need modification for other datasets)


#======================#
# Generate p'(a)       #

predicted_prev_func <- function(data, par){
  sp <- par[1]
  se <- par[2]
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda_par <- par[3:length(par)]
  # Repeat each site specific lambda for each age group in each dataset
  lambda_all <- lambda_par[data_joint$ID]
  # Prediction
  tp <-  1 - exp(-lambda_all * data_joint$age_month)  # 'True' (modelled) prevalence - p(a) - as a function of the catalytic model 
  op <- (1-sp) + (se+sp-1)*tp                         # Observed prevalence - p'(a) given by Diggle et al 2011 Epi Res Int
  op
}


#======================#
#  Binomial Likelihood #

# produce a single loglikelihood to evaluate from multiple site-specific lambda values#
loglike_simple <- function(data, par){
  predicted_seroprevalence = predicted_prev_func(data, par)    # produce a vector of predicted seroprev (by age x dataset)
  sum(dbinom(data_joint$x, data_joint$n, predicted_seroprevalence, log=T)) # produces individual binomail LL for each site specific lambda for age in each dataset, then SUMMED (x= positive pigs, n= total pigs)
}


#======================#
#        prior         #

prior <- function(par) {
  # Sp and Se are the first to parameters inthe vector of parameters
  sp <- par[1]
  se <- par[2]
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda_par <- par[3:length(par)]
  
  lambda_prior = sum(dunif(lambda_par, min = 0.0001, max = 1, log = T)) # uniform prior distributions
  se_prior = dbeta(se,9.5, 1.198198, log = T)                           # beta prior distributions
  sp_prior = dbeta(sp,38.6, 41.31718, log = T)                          # beta prior distributions
  
  return(sum(c(lambda_prior, se_prior, sp_prior)))
}

#======================#
#   Posterior          #

posterior_function_simple <- function(data, par){
  loglike_simple(data, par) + prior(par)
}

#======================#
#   Proposal           #

proposal_function_simple <- function(par, cov) {
  
  ## draw propopsals all from a multivariate normal 
  repeat {
    proposed <- mvrnorm(1, par, cov)
    if(all(proposed>0 & all(proposed[1:length(par)]<1))){break} 
  }
  
  return(proposed)
  
}  

#======================#
#    MCMC              #

# Run MCMC Function 
MCMC_simple_model <- function(inits,  number_of_iterations, cov) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol=length(inits))   # ncol is number of parameters
  MCMC_output[1,] <- inits
  Acceptances <- vector(length = number_of_iterations + 1)
  LogLikelihood_storage <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    proposed_parameter_value <- proposal_function_simple(MCMC_output[i,], cov)           # new proposed paramater value(s) with a given s.d. (step size)
    
    current_likelihood <- loglike_simple(data=data, MCMC_output[i,])                     # likelihood 
    
    current_posterior <- posterior_function_simple(data=data, MCMC_output[i,])           # current posterior from MCMC
    
    proposed_posterior <- posterior_function_simple(data=data, proposed_parameter_value) # proposed posterior with new proposed par value
    
    likelihood_ratio = exp(proposed_posterior - current_posterior);
    
    if(i %% (number_of_iterations/20) == 0){
      message(round(i/number_of_iterations*100), ' % completed')
    }
    
    if(runif(1) < likelihood_ratio) {
      # likelihood ratio comparison step (exponentiated because on log scale) 
      MCMC_output[i + 1,] <- proposed_parameter_value
      Acceptances[i] <- 1
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- proposed_posterior
      
    } else{
      
      MCMC_output[i + 1,] <- MCMC_output[i,]
      Acceptances[i] <- 0
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- current_posterior
      
    }
    
  } 
  
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Likelihood_Output"]] <- LogLikelihood_storage
  list[["Posterior_Output"]] <- Logposterior_storage
  return(list)
  
}

inits1 <- c(0.4, 0.9, 0.1, 0.001, 0.1, 0.001, 0.1)  # initial starting values for (sp, se, lambda_site_1: lambda_site_n)
inits2 <- c(0.8, 0.6, 0.001, 0.1, 0.001, 0.1, 0.001)


sd <- 0.0035                                        # standard deviation of proposal distribution 
cov <- diag(sd^2, 2+length(unique(data_joint$ID)))   # covariance matrix from multivariate proposal distribution

niter <- 1000000 # number of iterations
burnin <- 100000 # burnin

#=======================#
#      replicate        #

# set.seed() # Uncomment and add seed for determinism in MCMC
#sessionInfo()

#=======================#
#   Run MCMC            #

simple_out_chain1 <- MCMC_simple_model(inits1, niter, cov)  # initiate the MCMC
simple_out_chain2 <- MCMC_simple_model(inits2, niter, cov)  # initiate the MCMC

#==============================#
#     Output                   #

# acceptance ratio (target ~ 0.25)
sum(simple_out_chain1$Acceptances)/niter
sum(simple_out_chain2$Acceptances)/niter

# Chains plots for a) sens and spec and b) dataset-specific lambdas (need to edit) #
par(mfrow=(c(1,length(inits1[1:2]))))
for (i in 1:length(inits1[1:2])) {
  if (i==1) {
    ylab="specificty"
  } else {
    ylab="sensitivity"
  } 
  plot(simple_out_chain1$MCMC_Output[,i], type = "l", ylab=ylab, xlab="iter")
  lines(simple_out_chain2$MCMC_Output[,i], col="red")
}

par(mfrow=(c(1,length(unique(data_joint$ID)))))
for (i in 3:(length(unique(data_joint$ID))+2)) {
  if (i==3) {
    ylab="lambda (site 1)"
  } else if (i==4) {
    ylab="lambda (site 2)"
  } else if (i==5) {
    ylab="lambda (site 3)"
  } else if (i==6) {
    ylab="lambda (site 4)"
  } else {
    ylab="lambda (site 5)"
  }
  plot(simple_out_chain1$MCMC_Output[,i], type = "l", ylab=ylab, xlab="iter")
  lines(simple_out_chain2$MCMC_Output[,i], col="red")
}

# histogram plots for a) sens and spec and b) dataset-specific lambdas (need to edit) #

par(mfrow=(c(1,length(inits1[1:2]))))
for (i in 1:length(inits1[1:2])) {
  if (i==1) {
    ylab="specificty"
  } else {
    ylab="sensitivity"
  } 
  hist(c(simple_out_chain1$MCMC_Output[burnin:niter,i],simple_out_chain2$MCMC_Output[burnin:niter,i]), 
       xlab = ylab, main="", breaks=100)
}

par(mfrow=(c(1,length(unique(data_joint$ID)))))
for (i in 3:(length(unique(data_joint$ID))+2)) {
  if (i==3) {
    ylab="lambda (site 1)"
  } else if (i==4) {
    ylab="lambda (site 2)"
  } else if (i==5) {
    ylab="lambda (site 3)"
  } else if (i==6) {
    ylab="lambda (site 4)"
  } else {
    ylab="lambda (site 5)"
  }
  hist(c(simple_out_chain1$MCMC_Output[burnin:niter,i],simple_out_chain2$MCMC_Output[burnin:niter,i]), 
       xlab = ylab, main="", breaks=50)
}

#==================#
# autocorrelation  # 

# chain 1
par(mfrow=c(2,4))
acf(tail(simple_out_chain1$MCMC_Output[,1], 9500)) # autocorrelation function (y-axis) ~ lag (x-axis) plot for each parameter/chain
acf(tail(simple_out_chain1$MCMC_Output[,2], 9500)) # want ACF to drop with increasing k (or lag)
acf(tail(simple_out_chain1$MCMC_Output[,3], 9500)) # see tutorial for info: http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html 
acf(tail(simple_out_chain1$MCMC_Output[,4], 9500))
acf(tail(simple_out_chain1$MCMC_Output[,5], 9500))
acf(tail(simple_out_chain1$MCMC_Output[,6], 9500))
acf(tail(simple_out_chain1$MCMC_Output[,7], 9500))

# chain 2
par(mfrow=c(2,4))
acf(tail(simple_out_chain2$MCMC_Output[,1], 9500))
acf(tail(simple_out_chain2$MCMC_Output[,2], 9500))
acf(tail(simple_out_chain2$MCMC_Output[,3], 9500))
acf(tail(simple_out_chain2$MCMC_Output[,4], 9500))
acf(tail(simple_out_chain2$MCMC_Output[,5], 9500))
acf(tail(simple_out_chain2$MCMC_Output[,6], 9500))
acf(tail(simple_out_chain2$MCMC_Output[,7], 9500))

# calculate  autocorrelation within chains # 
autocor.sp.simple.catalytic.t1 <- cor(simple_out_chain1$MCMC_Output[,1][-1],simple_out_chain1$MCMC_Output[,1][-length(simple_out_chain1$MCMC_Output[,1])])
autocor.se.simple.catalytic.t1 <- cor(simple_out_chain1$MCMC_Output[,2][-1],simple_out_chain1$MCMC_Output[,2][-length(simple_out_chain1$MCMC_Output[,2])])
autocor.lambda1.simple.catalytic.t1 <- cor(simple_out_chain1$MCMC_Output[,3][-1],simple_out_chain1$MCMC_Output[,3][-length(simple_out_chain1$MCMC_Output[,3])])
autocor.lambda2.simple.catalytic.t1 <- cor(simple_out_chain1$MCMC_Output[,4][-1],simple_out_chain1$MCMC_Output[,4][-length(simple_out_chain1$MCMC_Output[,4])])
autocor.lambda3.simple.catalytic.t1 <- cor(simple_out_chain1$MCMC_Output[,5][-1],simple_out_chain1$MCMC_Output[,5][-length(simple_out_chain1$MCMC_Output[,5])])
autocor.lambda4.simple.catalytic.t1 <- cor(simple_out_chain1$MCMC_Output[,6][-1],simple_out_chain1$MCMC_Output[,6][-length(simple_out_chain1$MCMC_Output[,6])])
autocor.lambda5.simple.catalytic.t1 <- cor(simple_out_chain1$MCMC_Output[,7][-1],simple_out_chain1$MCMC_Output[,7][-length(simple_out_chain1$MCMC_Output[,7])])

# plot autocorrelation # 
par(mfrow=c(2,4))
plot(simple_out_chain1$MCMC_Output[,1][-1],simple_out_chain1$MCMC_Output[,1][-length(simple_out_chain1$MCMC_Output[,1])],main=paste("sp from simple model autocorrelation =",signif(autocor.sp.simple.catalytic.t1,3)),xlab="",ylab="")
plot(simple_out_chain1$MCMC_Output[,2][-1],simple_out_chain1$MCMC_Output[,2][-length(simple_out_chain1$MCMC_Output[,2])],main=paste("se from simple model autocorrelation =",signif(autocor.se.simple.catalytic.t1,3)),xlab="",ylab="")
plot(simple_out_chain1$MCMC_Output[,3][-1],simple_out_chain1$MCMC_Output[,3][-length(simple_out_chain1$MCMC_Output[,3])],main=paste("lambda (dataset 1) from simple model autocorrelation =",signif(autocor.lambda1.simple.catalytic.t1,3)),xlab="",ylab="")
plot(simple_out_chain1$MCMC_Output[,4][-1],simple_out_chain1$MCMC_Output[,4][-length(simple_out_chain1$MCMC_Output[,4])],main=paste("lambda (dataset 2) from simple model autocorrelation =",signif(autocor.lambda2.simple.catalytic.t1,3)),xlab="",ylab="")
plot(simple_out_chain1$MCMC_Output[,5][-1],simple_out_chain1$MCMC_Output[,5][-length(simple_out_chain1$MCMC_Output[,5])],main=paste("lambda (dataset 3) from simple model autocorrelation =",signif(autocor.lambda3.simple.catalytic.t1,3)),xlab="",ylab="")
plot(simple_out_chain1$MCMC_Output[,6][-1],simple_out_chain1$MCMC_Output[,6][-length(simple_out_chain1$MCMC_Output[,6])],main=paste("lambda (dataset 4) from simple model autocorrelation =",signif(autocor.lambda4.simple.catalytic.t1,3)),xlab="",ylab="")
plot(simple_out_chain1$MCMC_Output[,7][-1],simple_out_chain1$MCMC_Output[,7][-length(simple_out_chain1$MCMC_Output[,7])],main=paste("lambda (dataset 5) from simple model autocorrelation =",signif(autocor.lambda4.simple.catalytic.t1,3)),xlab="",ylab="")

#=====================================================================#
# correlation  (large number of combination so choose a selection)    # 

## sp ~ se ##
par(mfrow=c(1,1))
cor.t1_simple.spse <- cor(simple_out_chain1$MCMC_Output[, 1],simple_out_chain1$MCMC_Output[, 2])
plot(simple_out_chain1$MCMC_Output[, 1],simple_out_chain1$MCMC_Output[, 2],main=paste("Correlation for sp ~ se (chain 1, simple model)=",signif(cor.t1_simple.spse,3)),xlab="sp",ylab="se")
cor.t1_simple.spse

## sp ~ lambda site 1 ##
cor.t1_simple.splambda1 <- cor(simple_out_chain1$MCMC_Output[, 1],simple_out_chain1$MCMC_Output[, 3])
plot(simple_out_chain1$MCMC_Output[, 1],simple_out_chain1$MCMC_Output[, 3],main=paste("Correlation for sp ~ lambda site 1 (chain 1, simple model)=",signif(cor.t1_simple.splambda1,3)),xlab="sp",ylab="lambda site 1")
cor.t1_simple.splambda1

## se ~ lambda site 1 ##
cor.t1_simple.selambda1 <- cor(simple_out_chain1$MCMC_Output[, 2],simple_out_chain1$MCMC_Output[, 3])
plot(simple_out_chain1$MCMC_Output[, 2],simple_out_chain1$MCMC_Output[, 3],main=paste("Correlation for se ~ lambda site 1 (chain 1, simple model)=",signif(cor.t1_simple.splambda1,3)),xlab="se",ylab="lambda site 1")
cor.t1_simple.splambda1

# repeat/modify above code to look at correlation between other parameter combinations #

#======================#
# Processing of chains #

chains1_output <- simple_out_chain1$MCMC_Output
chains2_output <- simple_out_chain2$MCMC_Output

# functions to call #
plot_chains<-function(run1, run2){
  par(mfrow=c(ncol(run1),1))
  
  for(i in 1:ncol(run1)){
    plot(run1[,i], t='l', col='deeppink',
         ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
         xlab='', ylab=paste('Parameter', i, sep=' '))
    lines(run2[,i], col='dodgerblue')
  }
  
}

Burn<-function(chains, burnin){
  chains[-(1:burnin),]
}

Downsample<-function(chains, sample){
  chains[seq(1, nrow(chains), sample),]
}

Process_chains<-function(run1, run2, burnin, sample){
  C1<-Burn(run1, burnin)
  C2<-Burn(run2, burnin)
  
  
  C1<-Downsample(C1, sample)
  C2<-Downsample(C2, sample)
  
  
  return(list(C1, C2))
}

# Processing of chains function (# modify burnin and sub-sampling (reduce memory requirement of chains & autocorrelation))
PC_simple <-Process_chains(chains1_output, chains2_output, burnin=100000, sample=250)
plot_chains(PC_simple[[1]], PC_simple[[2]])

# replot ACF ~ lag to check autocorrelation after thinning #
# chain 1
par(mfrow=c(2,4))
acf(tail(PC_simple[[1]][,1], 9500)) # autocorrelation function (y-axis) ~ lag (x-axis) plot for each parameter/chain
acf(tail(PC_simple[[1]][,2], 9500)) # want ACF to drop with increasing k (or lag)
acf(tail(PC_simple[[1]][,3], 9500)) # see tutorial for info: http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html 
acf(tail(PC_simple[[1]][,4], 9500))
acf(tail(PC_simple[[1]][,5], 9500))
acf(tail(PC_simple[[1]][,6], 9500))
acf(tail(PC_simple[[1]][,7], 9500))

# chain 2
par(mfrow=c(2,4))
acf(tail(PC_simple[[2]][,1], 9500)) # autocorrelation function (y-axis) ~ lag (x-axis) plot for each parameter/chain
acf(tail(PC_simple[[2]][,2], 9500)) # want ACF to drop with increasing k (or lag)
acf(tail(PC_simple[[2]][,3], 9500)) # see tutorial for info: http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html 
acf(tail(PC_simple[[2]][,4], 9500))
acf(tail(PC_simple[[2]][,5], 9500))
acf(tail(PC_simple[[2]][,6], 9500))
acf(tail(PC_simple[[2]][,7], 9500))

autocor.sp.simple.t1 <- cor(PC_simple[[1]][,1][-1],PC_simple[[1]][,1][-length(PC_simple[[1]][,1])])
autocor.se.simple.t1 <- cor(PC_simple[[1]][,2][-1],PC_simple[[1]][,2][-length(PC_simple[[1]][,2])])
autocor.lambda1.simple.t1 <- cor(PC_simple[[1]][,3][-1],PC_simple[[1]][,3][-length(PC_simple[[1]][,3])])
autocor.lambda2.simple.t1 <- cor(PC_simple[[1]][,4][-1],PC_simple[[1]][,4][-length(PC_simple[[1]][,4])])
autocor.lambda3.simple.t1 <- cor(PC_simple[[1]][,5][-1],PC_simple[[1]][,5][-length(PC_simple[[1]][,5])])
autocor.lambda4.simple.t1 <- cor(PC_simple[[1]][,6][-1],PC_simple[[1]][,6][-length(PC_simple[[1]][,6])])
autocor.lambda5.simple.t1 <- cor(PC_simple[[1]][,7][-1],PC_simple[[1]][,7][-length(PC_simple[[1]][,7])])

signif(autocor.sp.simple.t1,3)
signif(autocor.se.simple.t1,3)
signif(autocor.lambda1.simple.t1,3)
signif(autocor.lambda2.simple.t1,3)
signif(autocor.lambda3.simple.t1,3)
signif(autocor.lambda4.simple.t1,3)
signif(autocor.lambda5.simple.t1,3)

#====================================================#
#   Post processing posterior plotting               #

# best parameter point estimates from distributions and credible intervals ##
par(mfrow=c(1,4))
sp_simple<-quantile(c(PC_simple[[1]][,1], PC_simple[[2]][,1]), c(0.025,0.5,0.975))
se_simple<-quantile(c(PC_simple[[1]][,2], PC_simple[[2]][,2]), c(0.025,0.5,0.975))
lambda1_simple<-quantile(c(PC_simple[[1]][,3], PC_simple[[2]][,3]), c(0.025,0.5,0.975))
lambda2_simple<-quantile(c(PC_simple[[1]][,4], PC_simple[[2]][,4]), c(0.025,0.5,0.975))
lambda3_simple<-quantile(c(PC_simple[[1]][,5], PC_simple[[2]][,5]), c(0.025,0.5,0.975))
lambda4_simple<-quantile(c(PC_simple[[1]][,6], PC_simple[[2]][,6]), c(0.025,0.5,0.975))
lambda5_simple<-quantile(c(PC_simple[[1]][,7], PC_simple[[2]][,7]), c(0.025,0.5,0.975))

sp_simple
se_simple
lambda1_simple
lambda2_simple
lambda3_simple
lambda4_simple
lambda5_simple

# only medians
sp.median<-quantile(c(PC_simple[[1]][,1], PC_simple[[2]][,1]), c(0.5))
se.median<-quantile(c(PC_simple[[1]][,2], PC_simple[[2]][,2]), c(0.5))
lambda1.median<-quantile(c(PC_simple[[1]][,3], PC_simple[[2]][,3]), c(0.5))
lambda2.median<-quantile(c(PC_simple[[1]][,4], PC_simple[[2]][,4]), c(0.5))
lambda3.median<-quantile(c(PC_simple[[1]][,5], PC_simple[[2]][,5]), c(0.5))
lambda4.median<-quantile(c(PC_simple[[1]][,6], PC_simple[[2]][,6]), c(0.5))
lambda5.median<-quantile(c(PC_simple[[1]][,7], PC_simple[[2]][,7]), c(0.5))

# plot posterior distributions #
par(mfrow=c(1,2))
hist(c(PC_simple[[1]][,1], PC_simple[[2]][,1]), breaks=30, xlab='specifcity')  # Parameter 1 - spec

hist(c(PC_simple[[1]][,2], PC_simple[[2]][,2]), breaks=30, xlab='sensitivity')      # Parameter 2 - sens

par(mfrow=c(1,5))
hist(c(PC_simple[[1]][,3], PC_simple[[2]][,3]), breaks=30, xlab='lambda (site 1)')      # Parameter 2 - lambda (site 1)

hist(c(PC_simple[[1]][,4], PC_simple[[2]][,4]), breaks=30, xlab='lambda (site 2)')      # Parameter 2 - lambda (site 2)

hist(c(PC_simple[[1]][,5], PC_simple[[2]][,5]), breaks=30, xlab='lambda (site 3)')      # Parameter 2 - lambda (site 3)

hist(c(PC_simple[[1]][,6], PC_simple[[2]][,6]), breaks=30, xlab='lambda (site 4)')      # Parameter 2 - lambda (site 4)

hist(c(PC_simple[[1]][,7], PC_simple[[2]][,7]), breaks=30, xlab='lambda (site 5)')      # Parameter 2 - lambda (site 5)


#===================================================#
#     Predicted prevalence curves                   #

# extract individual datasets from main data #
n <- length(unique(data_joint$ID))
eval(parse(text = paste0("data", seq(1:n), " <- ", split(data_joint, data_joint$ID))))
eval(parse(text = paste0("data", seq(1:n), " <-  as.data.frame(data", seq(1:5), ")"))) # change seq to 1: number of individual datasets

# specify new predicted prev function (for p'(a)) to work in dataframe below # 

predicted_prev_func_singledataset <- function(age, par){
  # Sp and Se are the first to parameters in the vector of parameters
  sp <- par[1]
  se <- par[2]
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda <- par[3]
  # Prediction
  tp <-  1 - exp(-lambda * age)  #'True' (modelled) prevalence - p(a) - as a function of the catalytic model 
  op <- (1-sp) + (se+sp-1)*tp          # Observed prevalence - p'(a) given by Diggle et al 2011 Epi Res Int
  op
}

# Predicted prevalence curve using posterior median point estimates (for each dataset) #
age_month <- seq(from=0, to=40, by=0.005)  # ages to produce predicted prevalence over

fitted_curve_df <- as.data.frame(age_month) 

predicted_simple1 <- full_join(fitted_curve_df, data1) 
predicted_simple2 <- full_join(fitted_curve_df, data2) 
predicted_simple3 <- full_join(fitted_curve_df, data3) 
predicted_simple4 <- full_join(fitted_curve_df, data4) 
predicted_simple5 <- full_join(fitted_curve_df, data5) 

predicted_simple1$predicted <- sapply(1:nrow(predicted_simple1), function(i) predicted_prev_func_singledataset(age=predicted_simple1$age_month[i], c(sp.median, se.median, lambda1.median)))
predicted_simple2$predicted <- sapply(1:nrow(predicted_simple2), function(i) predicted_prev_func_singledataset(age=predicted_simple2$age_month[i], c(sp.median, se.median, lambda2.median)))
predicted_simple3$predicted <- sapply(1:nrow(predicted_simple3), function(i) predicted_prev_func_singledataset(age=predicted_simple3$age_month[i], c(sp.median, se.median, lambda3.median)))
predicted_simple4$predicted <- sapply(1:nrow(predicted_simple4), function(i) predicted_prev_func_singledataset(age=predicted_simple4$age_month[i], c(sp.median, se.median, lambda4.median)))
predicted_simple5$predicted <- sapply(1:nrow(predicted_simple5), function(i) predicted_prev_func_singledataset(age=predicted_simple5$age_month[i], c(sp.median, se.median, lambda5.median)))

predicted_simple1$dataset_name <- rep(as.factor("data1"))
predicted_simple2$dataset_name <- rep(as.factor("data2"))
predicted_simple3$dataset_name <- rep(as.factor("data3"))
predicted_simple4$dataset_name <- rep(as.factor("data4"))
predicted_simple5$dataset_name <- rep(as.factor("data5"))

predicted_simple <- rbind(predicted_simple1, predicted_simple2, predicted_simple3, predicted_simple4, predicted_simple5) # master dataset

#================================================================================================================#
# create uncertainty around point estimates (using credible interval) of model run by subsampling from posterior # 

# dataset 1 #
subsampled_model <- matrix(NA, nrow = length(PC_simple[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_simple[[1]][,1])){
  
  single_model_output <- predicted_prev_func_singledataset(seq(0, 40, 0.005),c(PC_simple[[1]][i,1],PC_simple[[1]][i,2],PC_simple[[1]][i,3]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
simple_CrI1 <- cbind(Lower_sub, Upper_sub)
simple_CrI1 <- as.data.frame(simple_CrI1)
simple_CrI1$age_month <- fitted_curve_df$age_month
simple_CrI1$dataset_name <- rep(as.factor("data1"))


# dataset 2 #
subsampled_model <- matrix(NA, nrow = length(PC_simple[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_simple[[1]][,1])){
  
  single_model_output <- predicted_prev_func_singledataset(seq(0, 40, 0.005),c(PC_simple[[1]][i,1],PC_simple[[1]][i,2],PC_simple[[1]][i,4]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
simple_CrI2 <- cbind(Lower_sub, Upper_sub)
simple_CrI2 <- as.data.frame(simple_CrI2)
simple_CrI2$age_month <- fitted_curve_df$age_month
simple_CrI2$dataset_name <- rep(as.factor("data2"))


# dataset 3 #
subsampled_model <- matrix(NA, nrow = length(PC_simple[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_simple[[1]][,1])){
  
  single_model_output <- predicted_prev_func_singledataset(seq(0, 40, 0.005),c(PC_simple[[1]][i,1],PC_simple[[1]][i,2],PC_simple[[1]][i,5]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
simple_CrI3 <- cbind(Lower_sub, Upper_sub)
simple_CrI3 <- as.data.frame(simple_CrI3)
simple_CrI3$age_month <- fitted_curve_df$age_month
simple_CrI3$dataset_name <- rep(as.factor("data3"))


# dataset 4 #
subsampled_model <- matrix(NA, nrow = length(PC_simple[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_simple[[1]][,1])){
  
  single_model_output <- predicted_prev_func_singledataset(seq(0, 40, 0.005),c(PC_simple[[1]][i,1],PC_simple[[1]][i,2],PC_simple[[1]][i,6]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
simple_CrI4 <- cbind(Lower_sub, Upper_sub)
simple_CrI4 <- as.data.frame(simple_CrI4)
simple_CrI4$age_month <- fitted_curve_df$age_month
simple_CrI4$dataset_name <- rep(as.factor("data4"))

# dataset 5 #
subsampled_model <- matrix(NA, nrow = length(PC_simple[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_simple[[1]][,1])){
  
  single_model_output <- predicted_prev_func_singledataset(seq(0, 40, 0.005),c(PC_simple[[1]][i,1],PC_simple[[1]][i,2],PC_simple[[1]][i,7]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
simple_CrI5 <- cbind(Lower_sub, Upper_sub)
simple_CrI5 <- as.data.frame(simple_CrI5)
simple_CrI5$age_month <- fitted_curve_df$age_month
simple_CrI5$dataset_name <- rep(as.factor("data5"))


simple_CrI <- rbind(simple_CrI1, simple_CrI2, simple_CrI3, simple_CrI4, simple_CrI5) # master
names(simple_CrI) <- c("lower_credible", "upper_credible", "age_month", "dataset_name")

# plot #

ggplot() +    
  geom_line(data=predicted_simple,aes(x=age_month, y=predicted))+
  geom_ribbon(data=simple_CrI,aes(x=age_month,ymin=lower_credible, ymax=upper_credible),alpha=0.1)+
  geom_point(data=predicted_simple, aes(x=age_month, y=Observed_prevalence))+
  geom_errorbar(data=predicted_simple,aes(x=age_month, y=Observed_prevalence, ymin=lower, ymax=upper))+
  facet_wrap(~dataset_name, scales = "free")+
  ylim(0,1)+
  labs(x="Age (months) of pig host", y="(Sero)prevalence (0-1)", fill = "Model fitted to data", colour= "Model fitted to data", size="Sample size at 
       each age group")+
  theme_bw() +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=11, face= "bold"),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size = 18, hjust=0.5),
        axis.title.x = element_text(size = 18, face= "bold"),
        axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
        legend.position = c(0.83, 0.15),
        legend.title=element_text(size=20), 
        legend.text=element_text(size=16))