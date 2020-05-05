################################################
### PUBHLTH 250C: Bias Analysis Example Code
### Patrick T. Bradshaw
################################################

########## Code for Probabilistic Bias Analysis

expit <- function(x) exp(x)/(1+exp(x))


##### Example based on Gustafson and McCandless (2010), coffee and pancreatic cancer:
# Generate data from case-control study:
# Repeat relevant covariate pattern according to number in cells.
data.coffee <- rbind( matrix(rep(c(1,1), 347),ncol=2, byrow=T), # Exposed cases
                      matrix(rep(c(1,0), 20), ncol=2, byrow=T), # Unexposed cases
                      matrix(rep(c(0,1), 555), ncol=2, byrow=T), # Exposed controls
                      matrix(rep(c(0,0), 88), ncol=2, byrow=T)) # Unexposed controls
colnames(data.coffee) <- c("Y","X") # Give names
data.coffee <- data.frame(data.coffee) # Make data frame
data.coffee$id <- seq(1,nrow(data.coffee)) # create id variable

table(data.coffee$Y, data.coffee$X) # Check counts

N.obs <- nrow(data.coffee)

# Odds ratios from standard analysis:
exp(coef(glm(Y~X, data=data.coffee, family=binomial)))[2]

##################################################
##### Selection bias
##################################################
set.seed(123)

N.samp <- 10000

# Generate samples from bias parameters (these don't depend on iteration, 
# so can be done outside of loop).
theta.0 <- rnorm(N.samp, -5, 1)
theta.x <- rnorm(N.samp, -0.3, .1)
theta.y <- rnorm(N.samp, 8, 1.5)

# Initialize variable to store results:
OR.xy <- beta.xy <- rep(NA, N.samp)

for (i in 1:N.samp){
  # Sample U based on bias parameters for exposure effect:
  p.s <- expit(theta.0[i] + theta.x[i]*data.coffee$X + theta.y[i]*data.coffee$Y)
  w <- 1/p.s
  
  fit <- glm(Y ~ X, weights=w,
             data=data.coffee, 
             family=binomial)
  
  b.xy <- coef(fit)[2]; 
  se.b.xy <- sqrt(diag(vcov(fit)))[2];
  
  # Resample using estimated parameters (to add sampling variability)
  # and transform to OR.
  beta.xy[i] <- rnorm(1, b.xy, se.b.xy)
  OR.xy[i] <- exp(beta.xy[i]) # MCSA
}
# I keep finding an extreme OR (probably a bad set of weights)--this sifts that out:
quantile(OR.xy[OR.xy<1000], c(.5, 0.025, 0.975))
mean(beta.xy[beta.xy<1000])
sd(beta.xy[beta.xy<1000])

# pdf("MCSA_selection.pdf")
plot(density(OR.xy[OR.xy<1000]), main="Density of Bias-corrected OR (coffee and pancreatic cancer)",
     xlab="OR",ylab="Frequency")
# dev.off()







################################################
### PB HLTH 250C: Bayesian Bias Analysis Example Code
### Patrick T. Bradshaw
################################################

library(episensr)
library(R2jags)
library(coda)

##################################################
##### Unmeasured confounding
##################################################
##### Example based on Greenland (ME3), resins and lung cancer:
# Generate data from case-control study:
# Repeat relevant covariate pattern according to number in cells.
data.resins <- rbind( matrix(rep(c(1,1), 45),ncol=2, byrow=T), # Exposed cases
                      matrix(rep(c(1,0), 94), ncol=2, byrow=T), # Unexposed cases
                      matrix(rep(c(0,1), 257), ncol=2, byrow=T), # Exposed controls
                      matrix(rep(c(0,0), 945), ncol=2, byrow=T)) # Unexposed controls
colnames(data.resins) <- c("Y","X") # Give names
data.resins <- data.frame(data.resins) # Make data frame

data.resins$id <- seq(1,nrow(data.resins))

table(data.resins$Y, data.resins$X) # Check counts

N.obs <- nrow(data.resins)

### JAGS
library(R2jags)
confounding.bayes <- function() {
  # SAMPLING DISTRIBUTION
  for (i in 1:N) {
    logit(p[i]) <- b[1] + b[2]*X[i] + b[3]*U[i];
    Y[i] ~ dbern(p[i]);
    
    # Confounder model :
    p.U[i] <- p1*equals(X[i],1) + p0*equals(X[i],0);
    U[i] ~ dbern(p.U[i]);
  }
  
  # PRIORS ON BETAS
  b[1:2] ~ dmnorm(mu.b[1:2], tau.b[1:2,1:2]) # multivariate normal prior
  b[3] ~ dnorm(mu.b3, tau.b3)
  
  # Priors on prevalence:
  p1 ~ dunif(.6, .8);
  p0 ~ dunif(.4, .6);
  
  # Calculate ORs:
  for (l in 1:3) {
    OR[l] <- exp(b[l]);
  }
}

N <- nrow(data.resins);

# Parameters on the priors:
mu.b <- rep(0,2);          # Prior mean of beta
tau.b <- diag(.001,2);     # Prior precision

mu.b3 <- log(5);  # relationship between U and Y oods ratio
tau.b3 <- ((log(6) - log(4))/(2*1.96))^-2

# Lists for JAGS
data.conf <- list(N=nrow(data.resins), Y=data.resins$Y, X=data.resins$X,
                  mu.b=mu.b, tau.b=tau.b, mu.b3=mu.b3, tau.b3=tau.b3)
parameters.conf <- c("b", "OR","p1","p0") # Parameters to keep track of

set.seed(123)
# The following takes quite a bit of time:
sim.conf <- jags(data=data.conf,
                 parameters.to.save=parameters.conf,
                 n.iter=10000, model.file=confounding.bayes, 
                 n.thin=5, jags.seed=110410, n.chains = 2)

print(sim.conf,digits=4)

mcmc.conf <- as.mcmc(sim.conf)
plot(mcmc.conf, ask=F)

mcmc.conf.quantiles <- round(summary(mcmc.conf)$quantiles[, c(3,1,5)],2)
mcmc.conf.quantiles

# This matches the Bayesian parameterizaton above:
set.seed(123)
probsens.conf(case=data.resins$Y, exposed = data.resins$X,
              reps=10000,
              prev.exp = list("uniform", c(.6, .8)),
              prev.nexp = list("uniform",c(.4, .6)),
              risk = list("log-normal", c(mu.b3, sqrt(1/tau.b3))))