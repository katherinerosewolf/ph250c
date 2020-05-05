##### PB HLTH 250C: Advanced Epidemiologic Methods
##### R script for Bayesian Bias Analysis Lab Homework
##### 23 April 2020

library(R2jags)
library(epiR)
library(survival)
library(SurvRegCensCov)
library(epiR)
library(tictoc) # To time things

load("senssamp.rdata")
colnames(senssamp) <- tolower(colnames(senssamp))

# Correction to the dataset from last week's assignment:
senssamp$mortime2[senssamp$mortime2==0] <- .01 # Survival times need to be >1 

########## Standard analysis of less than definitive therapy on breast cancer mortality:
##### Frequentist Analysis (Optional, for comparison):
obs.model <- survreg(Surv(mortime2, bccause) ~ defnther + excat1 + agecat1 + agecat2,
                   data=senssamp, dist="weibull")
ConvertWeibull(obs.model)$HR

##### Standard Bayesian analysis (Weibull model):
jags.weibull <- function(){
     # SAMPLING DISTRIBUTION
     for (i in 1:N) {
          log(lambda[i]) <- b[1] + b[2]*defnther[i] + b[3]*excat1[i] + b[4]*agecat1[i] +
               b[5]*agecat2[i];
          
          censored[i] ~ dinterval(t[i], c[i]);
          t[i] ~ dweib(shape, lambda[i]);
     }
     
     # Prior on betas (log-HR):
     b[1:N.x] ~ dmnorm(mu.b[1:N.x], tau.b[1:N.x, 1:N.x]) # multivariate normal prior
     
     # Prior on shape parameter for Weibull:
     shape ~ dlnorm(0,8); # Log-normal
     
     # Calculate HRs:
     for (l in 1:N.x) {
          HR[l] <- exp(b[l]);
     }
}

# Parameters and priors
N.x <- 5                      # Number of predictors in outcome model
N <- nrow(senssamp)                  # Total number of observations
mu.b <- rep(0, N.x)           # Prior mean for regression parameters
tau.b <- diag(10^-2, N.x)     # Prior precision on regression parameters

# Outcome data:
# Define death time and censoring time:
senssamp$t <- senssamp$c <- senssamp$mortime2
senssamp$t[senssamp$bccause == 0] <- NA # If censored (death indicator = 0) then death time is missing
senssamp$c[senssamp$bccause == 1] <- 6 # If death (death = 1) then censoring time > death time
senssamp$censored <- 1 - senssamp$bccause # Indicator of censoring (=1)

# Data and initial values for JAGS
weibull.model.data <- list(N=N, N.x=N.x, 
                              mu.b=mu.b, tau.b=tau.b, 
                              t=senssamp$t, c=senssamp$c, censored=senssamp$censored, 
                              defnther=senssamp$defnther, excat1=senssamp$excat1,
                              agecat1=senssamp$agecat1, agecat2=senssamp$agecat2)

set.seed(123)
tic()
standard.weibull <- jags.parallel(model=jags.weibull, data=weibull.model.data, 
                         parameters.to.save = c("b","HR","shape"),
                         n.iter=100000, n.thin=10, n.chains=3,
                         jags.seed=123)
toc()
print(standard.weibull)
# plot(as.mcmc(standard.weibull), ask=T) # Uncomment for plots
# autocorr.plot(as.mcmc(standard.weibull))

############################################################
########## Bayesian bias analysis of unmeasured confounding
############################################################

##### Bias Analysis 1
jags.weibull.conf <- function(){
        for (i in 1:N) {
                # SAMPLING DISTRIBUTION
                
                # ***** MODIFY THE FOLLOWING EXPRESSION FOR THE HAZARD TO
                # ***** INCLUDE THE EFFECT OF THE UNMEASURED CONFOUNDER
                # ***** ASSUMING THE LOG-HAZARD RATIO IS CALLED b.u
                # ***** AND THE COVARIATE IS CALLED U:
                log(lambda[i]) <- b[1] + b[2]*defnther[i] + b[3]*excat1[i] +
                        b[4]*agecat1[i] + b[5]*agecat2[i];
                
                censored[i] ~ dinterval(t[i], c[i]);
                t[i] ~ dweib(shape, lambda[i]);
                
                # Distributon for the unmeasured confounder:
                pi.u[i] <- p.0*(1-defnther[i]) +   # Definitive therapy
                        p.1*defnther[i]         # Less than definitive therapy
                
                U[i] ~ dbin(pi.u[i],1) # Sample the unmeasured confounder
        }
        
        # Priors on betas:
        b[1:N.x] ~ dmnorm(mu.b[1:N.x], tau.b[1:N.x, 1:N.x]) # multivariate normal prior
        b.u ~ dnorm(mu.b.u, tau.b.u) # the prior for the bias parameter b.u
        
        # Shape parameter for Weibull:
        shape ~ dlnorm(0,8);
        
        # Priors for bias parameters:
        # ***** SPECIFY PRIORS FOR THE BIAS PARAMETERS
        # ***** THAT DEFINE THE DISTRIBUTION OF U WITH RESPECT
        # ***** TO THE EXPOSURE GROUPS.
        # ***** DEFINE THESE IN TERMS OF HYPERPARAMETERS:
        # ***** a.0, b.0 (for p.U0)
        # ***** a.1, b.1 (for p.U1)
        p.0 ~  # Definitive therapy, did not die
                p.1 ~  # Less than definitive therapy, did not die
                
                # Calculate HRs:
                for (l in 1:N.x) {
                        HR[l] <- exp(b[l]);
                }
        HR.u <- exp(b.u);
}

# ***** ACCORDING TO THE ABOVE DESCRIPTION OF THE BIAS PARAMETERS:
mu.b.u <- # ***** SPECIFY THE VALUE FOR THE MEAN OF THE PRIOR FOR b.u 
tau.b.u <- # ***** SPECIFY THE VALUE FOR THE PRECISION FOR THE PRIOR FOR b.u
        
# Uses epi.betabuster to obtain parameters for Beta distribution
# that satisfies the stated inputs (mode and lower quantile):
# ***** COMPLETE THE FOLLOWING TWO FUNCTION CALLS:
beta.p0 <- epi.betabuster(mode= , conf= , 
                          greaterthan = TRUE, x= )
beta.p1 <- epi.betabuster(mode= , conf= , 
                          greaterthan = TRUE, x= )

# Pulls shape parameters out of above:
a.0 <- beta.p0$shape1; b.0 <- beta.p0$shape2
a.1 <- beta.p1$shape1; b.1 <- beta.p1$shape2

# Confirm that above parameterizations yield distribution on bias
# parameters that are consistent with beliefs:
round(qlnorm(c(0.025, 0.975), mu.b.u, 1/sqrt(tau.b.u))) # Quantiles for prior HR.u
round(qbeta(c(.025, .975), a.0, b.0), 2) # Quantiles from prior for p.0
round(qbeta(c(.025, .975), a.1, b.1), 2) # Quantiles from prior for p.1

# Add these hyperparameters to the previous data list:  
weibull.model.data2 <- weibull.model.data
weibull.model.data2[["mu.b.u"]] <- mu.b.u 
weibull.model.data2[["tau.b.u"]] <- tau.b.u
weibull.model.data2[["a.0"]] <- a.0
weibull.model.data2[["b.0"]] <- b.0
weibull.model.data2[["a.1"]] <- a.1
weibull.model.data2[["b.1"]] <- b.1

set.seed(123)
tic()
conf.weibull <- jags.parallel(model=jags.weibull.conf, data=weibull.model.data2, 
                              parameters.to.save = c("b","b.u","HR","shape","HR.u",
                                                     "p.0","p.1"),
                              n.iter=100000, n.thin=5, n.chains=3,
                              jags.seed=123)
toc()
print(conf.weibull)
plot(as.mcmc(conf.weibull), ask=T)

##### Bias Analysis 2
# ***** ACCORDING TO THE ABOVE DESCRIPTION OF THE BIAS PARAMETERS:
mu.b.u.2 <- # ***** SPECIFY THE VALUE FOR THE MEAN OF THE PRIOR FOR b.u 
tau.b.u.2 <- # ***** SPECIFY THE VALUE FOR THE PRECISION FOR THE PRIOR FOR b.u

# Checking
round(qlnorm(c(0.025, 0.975), mu.b.u.2, 1/sqrt(tau.b.u.2)))

weibull.model.data3 <- weibull.model.data2
weibull.model.data3[["mu.b.u"]] <- mu.b.u.2 # refreshing hyperparameters
weibull.model.data3[["tau.b.u"]] <- tau.b.u.2

set.seed(123)
tic()
conf.weibull.2 <- jags.parallel(model=jags.weibull.conf, data=weibull.model.data3, 
                              parameters.to.save = c("b","b.u","HR","shape","HR.u",
                                                     "p.0","p.1"),
                              n.iter=100000, n.thin=10, n.chains=3,
                              jags.seed=123)
toc()
print(conf.weibull.2)
# plot(as.mcmc(conf.weibull.2), ask=T)

### Summarize:
# Standard model summary
standard.model <- summary(as.mcmc(standard.weibull))$quantiles["HR[2]", c(3,1,5)]
round(standard.model,2)

# Bias analysis summary
to.keep <- c("HR[2]", "HR.u", "p.0","p.1")

BA1 <- summary(as.mcmc(conf.weibull))$quantiles[to.keep, c(3,1,5)]
BA2 <- summary(as.mcmc(conf.weibull.2))$quantiles[to.keep, c(3,1,5)]

BA.all <- cbind(BA1, BA2)
colnames(BA.all) <- rep(c("Median","2.5%","97.5%"),3)
round(BA.all,2)
