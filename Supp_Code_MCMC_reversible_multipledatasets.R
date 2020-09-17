#================================================================================================================================#
#=========================== Reversible cataltyic model fit to multiple dataset (joint fitting)  (Dixon et al. ) ================#
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
  
  # define lambda final parameter position in parameter vector input
  lambda_end <- (length(unique(data_joint$ID))+2)
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda_par <- par[3:lambda_end]
  # Repeat each site specific lambda for each age group in each dataset
  lambda_all <- lambda_par[data_joint$ID]
  
  
  # define rho starting parameter position in parameter vector input
  rho_start <- (2+length(unique(data_joint$ID))+1)
  # define rho final parameter position in parameter vector input
  rho_end <- rho_start + (length(unique(data_joint$ID))-1)
  # Site-specific rhos are parameters 
  rho_par <- par[(rho_start:rho_end)]
  # Repeat each site specific lambda for each age group in each dataset
  rho_all <- rho_par[data_joint$ID]
  
  # Prediction
  tp <-  (lambda_all/(lambda_all + rho_all)) *(1 - exp(-(lambda_all + rho_all) * (data_joint$age_month))) # 'True' (modelled) prevalence - p(a) - as a function of the catalytic model 
  op <- (1-sp) + (se+sp-1)*tp                                                                             # Observed prevalence - p'(a) given by Diggle et al 2011 Epi Res Int
  op
}


#======================#
#  Binomial Likelihood #

# produce a single loglikelihood to evaluate from multiple site-specific lambda values#
loglike_reversible <- function(data, par){
  predicted_seroprevalence = predicted_prev_func(data, par)    # produce a vector of predicted seroprev (by age x dataset)
  sum(dbinom(data_joint$x, data_joint$n, predicted_seroprevalence, log=T)) # produces individual binomail LL for each site specific lambda for age in each dataset, then SUMMED (x= positive pigs, n= total pigs)
}


#======================#
#        prior         #

prior <- function(par) {
  # Sp and Se are the first to parameters inthe vector of parameters
  sp <- par[1]
  se <- par[2]
  # define lambda final parameter position in parameter vector input
  lambda_end <- (length(unique(data_joint$ID))+2)
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda_par <- par[3:lambda_end]
  
  # define rho starting parameter position in parameter vector input
  rho_start <- (2+length(unique(data_joint$ID))+1)
  # define rho final parameter position in parameter vector input
  rho_end <- rho_start + (length(unique(data_joint$ID))-1)
  # Site specific rhos
  rho_par <- par[(rho_start:rho_end)]
  
  # define prior distributions 
  lambda_prior = sum(dunif(lambda_par, min = 0.0001, max = 1, log = T)) # uniform prior distributions
  rho_prior = sum(dunif(rho_par, min = 0.00001, max = 1, log = T))      # Uniform prior distributions
  se_prior = dbeta(se,9.5, 1.198198, log = T)                           # beta prior distributions
  sp_prior = dbeta(sp,38.6, 41.31718, log = T)                          # beta prior distributions
  
  return(sum(c(lambda_prior, rho_prior, se_prior, sp_prior)))
}



#======================#
#   Posterior          #

posterior_function_reversible <- function(data, par){
  loglike_reversible(data, par) + prior(par)
}

#======================#
#   Proposal           #

proposal_function_reversible <- function(par, cov) {
  
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
MCMC_reversible_model <- function(inits,  number_of_iterations, cov) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol=length(inits))   # ncol is number of parameters
  MCMC_output[1,] <- inits
  Acceptances <- vector(length = number_of_iterations + 1)
  LogLikelihood_storage <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    proposed_parameter_value <- proposal_function_reversible(MCMC_output[i,], cov)           # new proposed paramater value(s) with a given s.d. (step size)
    
    current_likelihood <- loglike_reversible(data=data, MCMC_output[i,])                     # likelihood 
    
    current_posterior <- posterior_function_reversible(data=data, MCMC_output[i,])           # current posterior from MCMC
    
    proposed_posterior <- posterior_function_reversible(data=data, proposed_parameter_value) # proposed posterior with new proposed par value
    
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

inits1 <- c(0.4, 0.9, 0.1, 0.001, 0.1, 0.001, 0.1,0.1, 0.001, 0.1, 0.001, 0.1)  # initial starting values for (sp, se, lambda_site_1: lambda_site_n)
inits2 <- c(0.8, 0.6, 0.001, 0.1, 0.001, 0.1, 0.001, 0.01, 0.1, 0.01, 0.1, 0.001)


sd <- 0.005                                              # standard deviation of proposal distribution 
cov <- diag(sd^2, 2+(2*length(unique(data_joint$ID))))   # covariance matrix from multivariate proposal distribution

niter <- 1000000 # number of iterations
burnin <- 100000 # burnin

#=======================#
#      replicate        #

# set.seed() # Uncomment and add seed for determinism in MCMC
#sessionInfo()

#=======================#
#   Run MCMC            #

rev_out_chain1 <- MCMC_reversible_model(inits1, niter, cov)  # initiate the MCMC
rev_out_chain2 <- MCMC_reversible_model(inits2, niter, cov)  # initiate the MCMC

#==============================#
#     Output                   #

# acceptance ratio (target ~ 0.25)
sum(rev_out_chain1$Acceptances)/niter
sum(rev_out_chain2$Acceptances)/niter

# sens and spec chains #
par(mfrow=(c(1,length(inits1[1:2]))))
for (i in 1:length(inits1[1:2])) {
  if (i==1) {
    ylab="specificty"
  } else {
    ylab="sensitivity"
  } 
  plot(rev_out_chain1$MCMC_Output[,i], type = "l", ylab=ylab, xlab="iter")
  lines(rev_out_chain2$MCMC_Output[,i], col="red")
}

# lambda chains # 
lambda_end <- (length(unique(data_joint$ID))+2)# define lambda final parameter position in parameter vector input
par(mfrow=(c(1,length(unique(data_joint$ID)))))
for (i in 3:lambda_end) {
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
  plot(rev_out_chain1$MCMC_Output[,i], type = "l", ylab=ylab, xlab="iter")
  lines(rev_out_chain2$MCMC_Output[,i], col="red")
}

# rho chains #
rho_start <- (2+length(unique(data_joint$ID))+1)        # define rho starting parameter position in parameter vector input
rho_end <- rho_start + (length(unique(data_joint$ID))-1)# define rho final parameter position in parameter vector input

par(mfrow=(c(1,length(unique(data_joint$ID)))))
for (i in rho_start:rho_end) {
  if (i==8) {
    ylab="rho (site 1)"
  } else if (i==9) {
    ylab="rho (site 2)"
  } else if (i==10) {
    ylab="rho (site 3)"
  } else if (i==11) {
    ylab="rho (site 4)"
  } else {
    ylab="rho (site 5)"
  }
  plot(rev_out_chain1$MCMC_Output[,i], type = "l", ylab=ylab, xlab="iter")
  lines(rev_out_chain2$MCMC_Output[,i], col="red")
}

# histogram plots #

# sens and spec posteriors # 
par(mfrow=(c(1,length(inits1[1:2]))))
for (i in 1:length(inits1[1:2])) {
  if (i==1) {
    ylab="specificty"
  } else {
    ylab="sensitivity"
  } 
  hist(c(rev_out_chain1$MCMC_Output[burnin:niter,i],rev_out_chain2$MCMC_Output[burnin:niter,i]), 
       xlab = ylab, main="", breaks=100)
}

# lambda posteriors # 
par(mfrow=(c(1,length(unique(data_joint$ID)))))
for (i in 3:lambda_end) {
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
  hist(c(rev_out_chain1$MCMC_Output[burnin:niter,i],rev_out_chain2$MCMC_Output[burnin:niter,i]), 
       xlab = ylab, main="", breaks=50)
}

# rho posteriors # 
par(mfrow=(c(1,length(unique(data_joint$ID)))))
for (i in rho_start:rho_end) {
  if (i==8) {
    ylab="rho (site 1)"
  } else if (i==9) {
    ylab="rho (site 2)"
  } else if (i==10) {
    ylab="rho (site 3)"
  } else if (i==11) {
    ylab="rho (site 4)"
  } else {
    ylab="rho (site 5)"
  }
  hist(c(rev_out_chain1$MCMC_Output[burnin:niter,i],rev_out_chain2$MCMC_Output[burnin:niter,i]), 
       xlab = ylab, main="", breaks=50)
}
#==================#
# autocorrelation  # 

# chain 1
par(mfrow=c(3,4))
acf(tail(rev_out_chain1$MCMC_Output[,1], 9500)) # autocorrelation function (y-axis) ~ lag (x-axis) plot for each parameter/chain
acf(tail(rev_out_chain1$MCMC_Output[,2], 9500)) # want ACF to drop with increasing k (or lag)
acf(tail(rev_out_chain1$MCMC_Output[,3], 9500)) # see tutorial for info: http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html 
acf(tail(rev_out_chain1$MCMC_Output[,4], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,5], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,6], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,7], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,8], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,9], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,10], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,11], 9500))
acf(tail(rev_out_chain1$MCMC_Output[,12], 9500))

# chain 2
par(mfrow=c(3,4))
acf(tail(rev_out_chain2$MCMC_Output[,1], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,2], 9500)) 
acf(tail(rev_out_chain2$MCMC_Output[,3], 9500))  
acf(tail(rev_out_chain2$MCMC_Output[,4], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,5], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,6], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,7], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,8], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,9], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,10], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,11], 9500))
acf(tail(rev_out_chain2$MCMC_Output[,12], 9500))

# calculate  autocorrelation within chains # 
autocor.sp.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,1][-1],rev_out_chain1$MCMC_Output[,1][-length(rev_out_chain1$MCMC_Output[,1])])
autocor.se.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,2][-1],rev_out_chain1$MCMC_Output[,2][-length(rev_out_chain1$MCMC_Output[,2])])
autocor.lambda1.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,3][-1],rev_out_chain1$MCMC_Output[,3][-length(rev_out_chain1$MCMC_Output[,3])])
autocor.lambda2.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,4][-1],rev_out_chain1$MCMC_Output[,4][-length(rev_out_chain1$MCMC_Output[,4])])
autocor.lambda3.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,5][-1],rev_out_chain1$MCMC_Output[,5][-length(rev_out_chain1$MCMC_Output[,5])])
autocor.lambda4.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,6][-1],rev_out_chain1$MCMC_Output[,6][-length(rev_out_chain1$MCMC_Output[,6])])
autocor.lambda5.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,7][-1],rev_out_chain1$MCMC_Output[,7][-length(rev_out_chain1$MCMC_Output[,7])])
autocor.rho1.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,8][-1],rev_out_chain1$MCMC_Output[,8][-length(rev_out_chain1$MCMC_Output[,8])])
autocor.rho2.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,9][-1],rev_out_chain1$MCMC_Output[,9][-length(rev_out_chain1$MCMC_Output[,9])])
autocor.rho3.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,10][-1],rev_out_chain1$MCMC_Output[,10][-length(rev_out_chain1$MCMC_Output[,10])])
autocor.rho4.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,11][-1],rev_out_chain1$MCMC_Output[,11][-length(rev_out_chain1$MCMC_Output[,11])])
autocor.rho5.rev.catalytic.t1 <- cor(rev_out_chain1$MCMC_Output[,12][-1],rev_out_chain1$MCMC_Output[,12][-length(rev_out_chain1$MCMC_Output[,12])])


# plot autocorrelation # 
par(mfrow=c(1,2))
plot(rev_out_chain1$MCMC_Output[,1][-1],rev_out_chain1$MCMC_Output[,1][-length(rev_out_chain1$MCMC_Output[,1])],main=paste("sp from reversible model autocorrelation =",signif(autocor.sp.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,2][-1],rev_out_chain1$MCMC_Output[,2][-length(rev_out_chain1$MCMC_Output[,2])],main=paste("se from reversible model autocorrelation =",signif(autocor.se.rev.catalytic.t1,3)),xlab="",ylab="")

par(mfrow=c(2,3))
plot(rev_out_chain1$MCMC_Output[,3][-1],rev_out_chain1$MCMC_Output[,3][-length(rev_out_chain1$MCMC_Output[,3])],main=paste("lambda (dataset 1) from reversible model autocorrelation =",signif(autocor.lambda1.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,4][-1],rev_out_chain1$MCMC_Output[,4][-length(rev_out_chain1$MCMC_Output[,4])],main=paste("lambda (dataset 2) from reversible model autocorrelation =",signif(autocor.lambda2.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,5][-1],rev_out_chain1$MCMC_Output[,5][-length(rev_out_chain1$MCMC_Output[,5])],main=paste("lambda (dataset 3) from reversible model autocorrelation =",signif(autocor.lambda3.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,6][-1],rev_out_chain1$MCMC_Output[,6][-length(rev_out_chain1$MCMC_Output[,6])],main=paste("lambda (dataset 4) from reversible model autocorrelation =",signif(autocor.lambda4.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,7][-1],rev_out_chain1$MCMC_Output[,7][-length(rev_out_chain1$MCMC_Output[,7])],main=paste("lambda (dataset 5) from reversible model autocorrelation =",signif(autocor.lambda4.rev.catalytic.t1,3)),xlab="",ylab="")

par(mfrow=c(2,3))
plot(rev_out_chain1$MCMC_Output[,8][-1],rev_out_chain1$MCMC_Output[,8][-length(rev_out_chain1$MCMC_Output[,8])],main=paste("rho (dataset 1) from reversible model autocorrelation =",signif(autocor.rho1.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,9][-1],rev_out_chain1$MCMC_Output[,9][-length(rev_out_chain1$MCMC_Output[,9])],main=paste("rho (dataset 2) from reversible model autocorrelation =",signif(autocor.rho2.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,10][-1],rev_out_chain1$MCMC_Output[,10][-length(rev_out_chain1$MCMC_Output[,10])],main=paste("rho (dataset 3) from reversible model autocorrelation =",signif(autocor.rho3.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,11][-1],rev_out_chain1$MCMC_Output[,11][-length(rev_out_chain1$MCMC_Output[,11])],main=paste("rho (dataset 4) from reversible model autocorrelation =",signif(autocor.rho4.rev.catalytic.t1,3)),xlab="",ylab="")
plot(rev_out_chain1$MCMC_Output[,12][-1],rev_out_chain1$MCMC_Output[,12][-length(rev_out_chain1$MCMC_Output[,12])],main=paste("rho (dataset 5) from reversible model autocorrelation =",signif(autocor.rho4.rev.catalytic.t1,3)),xlab="",ylab="")

#=====================================================================#
# correlation  (large number of combination so choose a selection)    # 

## sp ~ se ##
par(mfrow=c(1,1))
cor.t1_rev.spse <- cor(rev_out_chain1$MCMC_Output[, 1],rev_out_chain1$MCMC_Output[, 2])
plot(rev_out_chain1$MCMC_Output[, 1],rev_out_chain1$MCMC_Output[, 2],main=paste("Correlation for sp ~ se (chain 1, reversible model)=",signif(cor.t1_rev.spse,3)),xlab="sp",ylab="se")
cor.t1_simple.spse

## sp ~ lambda site 1 ##
cor.t1_rev.splambda1 <- cor(rev_out_chain1$MCMC_Output[, 1],rev_out_chain1$MCMC_Output[, 3])
plot(rev_out_chain1$MCMC_Output[, 1],rev_out_chain1$MCMC_Output[, 3],main=paste("Correlation for sp ~ lambda site 1 (chain 1, reversible model)=",signif(cor.t1_rev.splambda1,3)),xlab="sp",ylab="lambda site 1")
cor.t1_simple.splambda1

## se ~ lambda site 1 ##
cor.t1_rev.selambda1 <- cor(rev_out_chain1$MCMC_Output[, 2],rev_out_chain1$MCMC_Output[, 3])
plot(rev_out_chain1$MCMC_Output[, 2],rev_out_chain1$MCMC_Output[, 3],main=paste("Correlation for se ~ lambda site 1 (chain 1, reversible model)=",signif(cor.t1_rev.splambda1,3)),xlab="se",ylab="lambda site 1")
cor.t1_rev.splambda1

## sp ~ rho site 1 ##
cor.t1_rev.sprho1 <- cor(rev_out_chain1$MCMC_Output[, 1],rev_out_chain1$MCMC_Output[, 8])
plot(rev_out_chain1$MCMC_Output[, 1],rev_out_chain1$MCMC_Output[, 8],main=paste("Correlation for sp ~ rho site 1 (chain 1, reversible model)=",signif(cor.t1_rev.sprho1,3)),xlab="sp",ylab="rho site 1")
cor.t1_rev.sprho1

## sp ~ rho site 1 ##
cor.t1_rev.serho1 <- cor(rev_out_chain1$MCMC_Output[, 2],rev_out_chain1$MCMC_Output[, 8])
plot(rev_out_chain1$MCMC_Output[, 2],rev_out_chain1$MCMC_Output[, 8],main=paste("Correlation for se ~ rho site 1 (chain 1, reversible model)=",signif(cor.t1_rev.serho1,3)),xlab="se",ylab="rho site 1")
cor.t1_rev.serho1

## lambda ~ rho site 1 ##
cor.t1_rev.lamrho1 <- cor(rev_out_chain1$MCMC_Output[, 3],rev_out_chain1$MCMC_Output[, 8])
plot(rev_out_chain1$MCMC_Output[, 3],rev_out_chain1$MCMC_Output[, 8],main=paste("Correlation for lambda ~ rho site 1 (chain 1, reversible model)=",signif(cor.t1_rev.lamrho1,3)),xlab="lambda site 1",ylab="rho site 1")
cor.t1_rev.lamrho1

# repeat/modify above code to look at correlation between other parameter combinations #

#======================#
# Processing of chains #

chains1_output <- rev_out_chain1$MCMC_Output
chains2_output <- rev_out_chain2$MCMC_Output

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
PC_rev <-Process_chains(chains1_output, chains2_output, burnin=100000, sample=500)
plot_chains(PC_rev[[1]], PC_rev[[2]])

# replot ACF ~ lag to check autocorrelation after thinning #
# chain 1
par(mfrow=c(3,4))
acf(tail(PC_rev[[1]][,1], 9500)) # autocorrelation function (y-axis) ~ lag (x-axis) plot for each parameter/chain
acf(tail(PC_rev[[1]][,2], 9500)) # want ACF to drop with increasing k (or lag)
acf(tail(PC_rev[[1]][,3], 9500)) # see tutorial for info: http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html 
acf(tail(PC_rev[[1]][,4], 9500))
acf(tail(PC_rev[[1]][,5], 9500))
acf(tail(PC_rev[[1]][,6], 9500))
acf(tail(PC_rev[[1]][,7], 9500))
acf(tail(PC_rev[[1]][,8], 9500))
acf(tail(PC_rev[[1]][,9], 9500))
acf(tail(PC_rev[[1]][,10], 9500))
acf(tail(PC_rev[[1]][,11], 9500))
acf(tail(PC_rev[[1]][,12], 9500))

# chain 2
par(mfrow=c(3,4))
acf(tail(PC_rev[[2]][,1], 9500)) # autocorrelation function (y-axis) ~ lag (x-axis) plot for each parameter/chain
acf(tail(PC_rev[[2]][,2], 9500)) # want ACF to drop with increasing k (or lag)
acf(tail(PC_rev[[2]][,3], 9500)) # see tutorial for info: http://sbfnk.github.io/mfiidd/mcmc_diagnostics.html 
acf(tail(PC_rev[[2]][,4], 9500))
acf(tail(PC_rev[[2]][,5], 9500))
acf(tail(PC_rev[[2]][,6], 9500))
acf(tail(PC_rev[[2]][,7], 9500))
acf(tail(PC_rev[[2]][,8], 9500))
acf(tail(PC_rev[[2]][,9], 9500))
acf(tail(PC_rev[[2]][,10], 9500))
acf(tail(PC_rev[[2]][,11], 9500))
acf(tail(PC_rev[[2]][,12], 9500))

autocor.sp.rev.t1 <- cor(PC_rev[[1]][,1][-1],PC_rev[[1]][,1][-length(PC_rev[[1]][,1])])
autocor.se.rev.t1 <- cor(PC_rev[[1]][,2][-1],PC_rev[[1]][,2][-length(PC_rev[[1]][,2])])
autocor.lambda1.rev.t1 <- cor(PC_rev[[1]][,3][-1],PC_rev[[1]][,3][-length(PC_rev[[1]][,3])])
autocor.lambda2.rev.t1 <- cor(PC_rev[[1]][,4][-1],PC_rev[[1]][,4][-length(PC_rev[[1]][,4])])
autocor.lambda3.rev.t1 <- cor(PC_rev[[1]][,5][-1],PC_rev[[1]][,5][-length(PC_rev[[1]][,5])])
autocor.lambda4.rev.t1 <- cor(PC_rev[[1]][,6][-1],PC_rev[[1]][,6][-length(PC_rev[[1]][,6])])
autocor.lambda5.rev.t1 <- cor(PC_rev[[1]][,7][-1],PC_rev[[1]][,7][-length(PC_rev[[1]][,7])])
autocor.rho1.rev.t1 <- cor(PC_rev[[1]][,8][-1],PC_rev[[1]][,8][-length(PC_rev[[1]][,8])])
autocor.rho2.rev.t1 <- cor(PC_rev[[1]][,9][-1],PC_rev[[1]][,9][-length(PC_rev[[1]][,9])])
autocor.rho3.rev.t1 <- cor(PC_rev[[1]][,10][-1],PC_rev[[1]][,10][-length(PC_rev[[1]][,10])])
autocor.rho4.rev.t1 <- cor(PC_rev[[1]][,11][-1],PC_rev[[1]][,11][-length(PC_rev[[1]][,11])])
autocor.rho5.rev.t1 <- cor(PC_rev[[1]][,12][-1],PC_rev[[1]][,12][-length(PC_rev[[1]][,12])])


signif(autocor.sp.rev.t1,3)
signif(autocor.se.rev.t1,3)
signif(autocor.lambda1.rev.t1,3)
signif(autocor.lambda2.rev.t1,3)
signif(autocor.lambda3.rev.t1,3)
signif(autocor.lambda4.rev.t1,3)
signif(autocor.lambda5.rev.t1,3)
signif(autocor.rho1.rev.t1,3)
signif(autocor.rho2.rev.t1,3)
signif(autocor.rho3.rev.t1,3)
signif(autocor.rho4.rev.t1,3)
signif(autocor.rho5.rev.t1,3)

#====================================================#
#   Post processing posterior plotting               #

# best parameter point estimates from distributions and credible intervals ##
par(mfrow=c(1,4))
sp_rev<-quantile(c(PC_rev[[1]][,1], PC_rev[[2]][,1]), c(0.025,0.5,0.975))
se_rev<-quantile(c(PC_rev[[1]][,2], PC_rev[[2]][,2]), c(0.025,0.5,0.975))
lambda1_rev<-quantile(c(PC_rev[[1]][,3], PC_rev[[2]][,3]), c(0.025,0.5,0.975))
lambda2_rev<-quantile(c(PC_rev[[1]][,4], PC_rev[[2]][,4]), c(0.025,0.5,0.975))
lambda3_rev<-quantile(c(PC_rev[[1]][,5], PC_rev[[2]][,5]), c(0.025,0.5,0.975))
lambda4_rev<-quantile(c(PC_rev[[1]][,6], PC_rev[[2]][,6]), c(0.025,0.5,0.975))
lambda5_rev<-quantile(c(PC_rev[[1]][,7], PC_rev[[2]][,7]), c(0.025,0.5,0.975))
rho1_rev<-quantile(c(PC_rev[[1]][,8], PC_rev[[2]][,8]), c(0.025,0.5,0.975))
rho2_rev<-quantile(c(PC_rev[[1]][,9], PC_rev[[2]][,9]), c(0.025,0.5,0.975))
rho3_rev<-quantile(c(PC_rev[[1]][,10], PC_rev[[2]][,10]), c(0.025,0.5,0.975))
rho4_rev<-quantile(c(PC_rev[[1]][,11], PC_rev[[2]][,11]), c(0.025,0.5,0.975))
rho5_rev<-quantile(c(PC_rev[[1]][,12], PC_rev[[2]][,12]), c(0.025,0.5,0.975))

sp_rev
se_rev
lambda1_rev
lambda2_rev
lambda3_rev
lambda4_rev
lambda5_rev
rho1_rev
rho2_rev
rho3_rev
rho4_rev
rho5_rev

# only medians
sp.median<-quantile(c(PC_rev[[1]][,1], PC_rev[[2]][,1]), c(0.5))
se.median<-quantile(c(PC_rev[[1]][,2], PC_rev[[2]][,2]), c(0.5))
lambda1.median<-quantile(c(PC_rev[[1]][,3], PC_rev[[2]][,3]), c(0.5))
lambda2.median<-quantile(c(PC_rev[[1]][,4], PC_rev[[2]][,4]), c(0.5))
lambda3.median<-quantile(c(PC_rev[[1]][,5], PC_rev[[2]][,5]), c(0.5))
lambda4.median<-quantile(c(PC_rev[[1]][,6], PC_rev[[2]][,6]), c(0.5))
lambda5.median<-quantile(c(PC_rev[[1]][,7], PC_rev[[2]][,7]), c(0.5))
rho1.median<-quantile(c(PC_rev[[1]][,8], PC_rev[[2]][,8]), c(0.5))
rho2.median<-quantile(c(PC_rev[[1]][,9], PC_rev[[2]][,9]), c(0.5))
rho3.median<-quantile(c(PC_rev[[1]][,10], PC_rev[[2]][,10]), c(0.5))
rho4.median<-quantile(c(PC_rev[[1]][,11], PC_rev[[2]][,11]), c(0.5))
rho5.median<-quantile(c(PC_rev[[1]][,12], PC_rev[[2]][,12]), c(0.5))

# plot posterior distributions #
par(mfrow=c(1,2))
hist(c(PC_rev[[1]][,1], PC_rev[[2]][,1]), breaks=30, xlab='specifcity')           # Parameter 1 - spec

hist(c(PC_rev[[1]][,2], PC_rev[[2]][,2]), breaks=30, xlab='sensitivity')          # Parameter 2 - sens

par(mfrow=c(2,3))
hist(c(PC_rev[[1]][,3], PC_rev[[2]][,3]), breaks=30, xlab='lambda (site 1)')      # Parameter 3 - lambda (site 1)

hist(c(PC_rev[[1]][,4], PC_rev[[2]][,4]), breaks=30, xlab='lambda (site 2)')      # Parameter 4 - lambda (site 2)

hist(c(PC_rev[[1]][,5], PC_rev[[2]][,5]), breaks=30, xlab='lambda (site 3)')      # Parameter 5 - lambda (site 3)

hist(c(PC_rev[[1]][,6], PC_rev[[2]][,6]), breaks=30, xlab='lambda (site 4)')      # Parameter 6 - lambda (site 4)

hist(c(PC_rev[[1]][,7], PC_rev[[2]][,7]), breaks=30, xlab='lambda (site 5)')      # Parameter 7 - lambda (site 5)

par(mfrow=c(2,3))
hist(c(PC_rev[[1]][,8], PC_rev[[2]][,8]), breaks=30, xlab='rho (site 1)')               # Parameter 8 - rho (site 1)

hist(c(PC_rev[[1]][,9], PC_rev[[2]][,9]), breaks=30, xlab='rho (site 2)')               # Parameter 9 - rho (site 2)

hist(c(PC_rev[[1]][,10], PC_rev[[2]][,10]), breaks=30, xlab='rho (site 3)')             # Parameter 10 - rho (site 3)

hist(c(PC_rev[[1]][,11], PC_rev[[2]][,11]), breaks=30, xlab='rho (site 4)')             # Parameter 11 - rho (site 4)

hist(c(PC_rev[[1]][,12], PC_rev[[2]][,12]), breaks=30, xlab='rho (site 5)')             # Parameter 12 - rho (site 5)

#===================================================#
#     Predicted prevalence curves                   #

# extract individual datasets from main data #
n <- length(unique(data_joint$ID))
eval(parse(text = paste0("data", seq(1:n), " <- ", split(data_joint, data_joint$ID))))
eval(parse(text = paste0("data", seq(1:n), " <-  as.data.frame(data", seq(1:5), ")"))) # change seq to 1: number of individual datasets

# specify new predicted prev function (for p'(a)) to work in dataframe below # 

predicted_prev_reversible_func_singledataset <- function(age, par){
  # Sp and Se are the first to parameters in the vector of parameters
  sp <- par[1]
  se <- par[2]
  # Site-specific lambda 
  lambda <- par[3]
  # Site-specific rho 
  rho <- par[4]
  
  # Prediction
  tp <-  (lambda/(lambda + rho)) *(1 - exp(-(lambda + rho) * age))  #'True' (modelled) prevalence - p(a) - as a function of the catalytic model 
  op <- (1-sp) + (se+sp-1)*tp                                       # Observed prevalence - p'(a) given by Diggle et al 2011 Epi Res Int
  op
}

# Predicted prevalence curve using posterior median point estimates (for each dataset) #
age_month <- seq(from=0, to=40, by=0.005)  # ages to produce predicted prevalence over

fitted_curve_df <- as.data.frame(age_month) 

predicted_rev1 <- full_join(fitted_curve_df, data1) 
predicted_rev2 <- full_join(fitted_curve_df, data2) 
predicted_rev3 <- full_join(fitted_curve_df, data3) 
predicted_rev4 <- full_join(fitted_curve_df, data4) 
predicted_rev5 <- full_join(fitted_curve_df, data5) 

predicted_rev1$predicted <- sapply(1:nrow(predicted_rev1), function(i) predicted_prev_reversible_func_singledataset(age=predicted_rev1$age_month[i], c(sp.median, se.median, lambda1.median, rho1.median)))
predicted_rev2$predicted <- sapply(1:nrow(predicted_rev2), function(i) predicted_prev_reversible_func_singledataset(age=predicted_rev2$age_month[i], c(sp.median, se.median, lambda2.median, rho2.median)))
predicted_rev3$predicted <- sapply(1:nrow(predicted_rev3), function(i) predicted_prev_reversible_func_singledataset(age=predicted_rev3$age_month[i], c(sp.median, se.median, lambda3.median, rho3.median)))
predicted_rev4$predicted <- sapply(1:nrow(predicted_rev4), function(i) predicted_prev_reversible_func_singledataset(age=predicted_rev4$age_month[i], c(sp.median, se.median, lambda4.median, rho4.median)))
predicted_rev5$predicted <- sapply(1:nrow(predicted_rev5), function(i) predicted_prev_reversible_func_singledataset(age=predicted_rev5$age_month[i], c(sp.median, se.median, lambda5.median, rho5.median)))

predicted_rev1$dataset_name <- rep(as.factor("data1"))
predicted_rev2$dataset_name <- rep(as.factor("data2"))
predicted_rev3$dataset_name <- rep(as.factor("data3"))
predicted_rev4$dataset_name <- rep(as.factor("data4"))
predicted_rev5$dataset_name <- rep(as.factor("data5"))

predicted_reversible <- rbind(predicted_rev1, predicted_rev2, predicted_rev3, predicted_rev4, predicted_rev5) # master dataset

#================================================================================================================#
# create uncertainty around point estimates (using credible interval) of model run by subsampling from posterior # 

# dataset 1 #
subsampled_model <- matrix(NA, nrow = length(PC_rev[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_rev[[1]][,1])){
  
  single_model_output <- predicted_prev_reversible_func_singledataset(seq(0, 40, 0.005),c(PC_rev[[1]][i,1],PC_rev[[1]][i,2],PC_rev[[1]][i,3],PC_rev[[1]][i,8]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
reversible_CrI1 <- cbind(Lower_sub, Upper_sub)
reversible_CrI1 <- as.data.frame(reversible_CrI1)
reversible_CrI1$age_month <- fitted_curve_df$age_month
reversible_CrI1$dataset_name <- rep(as.factor("data1"))


# dataset 2 #
subsampled_model <- matrix(NA, nrow = length(PC_rev[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_rev[[1]][,1])){
  
  single_model_output <- predicted_prev_reversible_func_singledataset(seq(0, 40, 0.005),c(PC_rev[[1]][i,1],PC_rev[[1]][i,2],PC_rev[[1]][i,4],PC_rev[[1]][i,9]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
reversible_CrI2 <- cbind(Lower_sub, Upper_sub)
reversible_CrI2 <- as.data.frame(reversible_CrI2)
reversible_CrI2$age_month <- fitted_curve_df$age_month
reversible_CrI2$dataset_name <- rep(as.factor("data2"))


# dataset 3 #
subsampled_model <- matrix(NA, nrow = length(PC_rev[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_rev[[1]][,1])){
  
  single_model_output <- predicted_prev_reversible_func_singledataset(seq(0, 40, 0.005),c(PC_rev[[1]][i,1],PC_rev[[1]][i,2],PC_rev[[1]][i,5],PC_rev[[1]][i,10]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
reversible_CrI3 <- cbind(Lower_sub, Upper_sub)
reversible_CrI3 <- as.data.frame(reversible_CrI3)
reversible_CrI3$age_month <- fitted_curve_df$age_month
reversible_CrI3$dataset_name <- rep(as.factor("data3"))


# dataset 4 #
subsampled_model <- matrix(NA, nrow = length(PC_rev[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_rev[[1]][,1])){
  
  single_model_output <- predicted_prev_reversible_func_singledataset(seq(0, 40, 0.005),c(PC_rev[[1]][i,1],PC_rev[[1]][i,2],PC_rev[[1]][i,6],PC_rev[[1]][i,11]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
reversible_CrI4 <- cbind(Lower_sub, Upper_sub)
reversible_CrI4 <- as.data.frame(reversible_CrI4)
reversible_CrI4$age_month <- fitted_curve_df$age_month
reversible_CrI4$dataset_name <- rep(as.factor("data4"))

# dataset 5 #
subsampled_model <- matrix(NA, nrow = length(PC_rev[[1]][,1]), ncol = length(seq(0, 40, 0.005))) # create matrix to store uncertainty output

for (i in 1:length(PC_rev[[1]][,1])){
  
  single_model_output <- predicted_prev_reversible_func_singledataset(seq(0, 40, 0.005),c(PC_rev[[1]][i,1],PC_rev[[1]][i,2],PC_rev[[1]][i,7],PC_rev[[1]][i,12]))
  subsampled_model[i, ] <- single_model_output
  
}

lower_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.025)
upper_credible_interval <- apply(subsampled_model, MARGIN = 2, quantile, prob = 0.975)
Lower_sub <- as.data.frame(lower_credible_interval)
Upper_sub <- as.data.frame(upper_credible_interval)
reversible_CrI5 <- cbind(Lower_sub, Upper_sub)
reversible_CrI5 <- as.data.frame(reversible_CrI5)
reversible_CrI5$age_month <- fitted_curve_df$age_month
reversible_CrI5$dataset_name <- rep(as.factor("data5"))


reversible_CrI <- rbind(reversible_CrI1, reversible_CrI2, reversible_CrI3, reversible_CrI4, reversible_CrI5) # master
names(reversible_CrI) <- c("lower_credible", "upper_credible", "age_month", "dataset_name")

# plot #

ggplot() +    
  geom_line(data=predicted_reversible,aes(x=age_month, y=predicted))+
  geom_ribbon(data=reversible_CrI,aes(x=age_month,ymin=lower_credible, ymax=upper_credible),alpha=0.1)+
  geom_point(data=predicted_reversible, aes(x=age_month, y=Observed_prevalence))+
  geom_errorbar(data=predicted_reversible,aes(x=age_month, y=Observed_prevalence, ymin=lower, ymax=upper))+
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