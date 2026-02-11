# Load libraries ----
library(R2jags)
library(tidyverse)

# Data ----
# The year is 2016 and these are the gusher data. 

# About 6 months ago my partner started buying Minion-themed Fruit Gushers 
# for the occasional snack for our oldest child. In keeping with the 
# Minions theme, the fruit snacks in each package come in two colors: 
# blue and yellow. After having personally dispatched no fewer than a 
# dozen of these from the goody cabinet in my garage, I am convinced that 
# there is a lower probability of getting a blue gusher than a yellow gusher. 
# Everybody knows that blue is better than yellow when it comes to sugary 
# snacks, with the possible exception of Warheads, so I feel like I have been 
# ripped off. Needless to say, this is the most important issue in my life 
# during an election year, and I think we deserve some answers before we 
# commit another $7.99 to General Mills and, by extension, the superfund that 
# is undoubtedly supporting one of the presidential candidates.

# So, the question I want answered during this lab is this: is the 
# probability of obtaining a blue Minion-themed Fruit Gusher (i.e. success) 
# significantly different than one would expect (0.5) if there is no bias at 
# the gusher factory in the distribution of blue and yellow gushers in these 
# packages. To answer this question, I bought a case of Minions Gushers 
# wholesale. We are going to house said case of snacks, catch a wicked 
# sugar rush, and determine if the house has stacked the deck against blue 
# Gushers in an effort to promote their candidate of choice 
# (whose constituency clearly prefers yellow over blue, boooooo!). We actually
# did do this in class.
#
# - The blue gusher bandit

# Read in the data
gushers = read.csv('data/gushers.csv')

# Have a look
head(gushers, 10)

# We have a dataframe with 36 observations of two variables. These data
# represent binomial trials. The number of successes ('blue') is provided for 
# total number of gushers in each package ('total').

# Specify a model in BUGS language ----
# Binomial logistic regresion model where theta is probability
# of success (proportion of blue gushers), and n is number of trials 
# (number of gushers in the package)
# The sink function is used to write a text string to an external model 
# file - another example of how you can do this.
sink("models/gushers_model.txt")          
cat("model {                       
      # Likelihood
        for(i in 1:n.packs){      
          y[i] ~ dbinom( theta, n[i] )   
        }                         

        theta ~ dbeta(a, b)       
        a <- 10                    
        b <- 10
    }",
    fill = TRUE)
sink()

# Now specify a data set for use in JAGS ----
g_data = list(
  y = gushers$blue,
  n.packs = nrow(gushers),
  n = gushers$total
)

# Parameters monitored ----
parameters <- c("theta")

# Provide initial values ----
inits <- function(){
  list(
    theta = runif(1, 0, 1)
    )
  }

# MCMC settings ----
ni <- 33000     # Number of draws from posterior (for each chain)
nt <- 3         # Thinning rate
nb <- 3000      # Number of draws to discard as burn-in
nc <- 3         # Number of chains

# Call jags and run the model ----
gush_model <- jags(g_data, inits=inits, parameters, "models/gushers_model.txt",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
                   working.directory = getwd())

# Print a summary of the model
print(gush_model)

# Have a look at the diagnostics ----
# Just like everything else, the model is an object with named
# elements:
names(gush_model)
names(gush_model$BUGSoutput)
names(gush_model$BUGSoutput$sims.list)

# We access the posterior distribution from the mcmc simulations
# like this
post_gush = gush_model$BUGSoutput$sims.list$theta

# We can also access the individual Markov chains
chain1 = post_gush[1:(length(post_gush)/nc)]
chain2 = post_gush[((length(post_gush)/nc)+1):((length(post_gush)/nc)*2)]
chain3 = post_gush[((length(post_gush)/nc)*2+1):((length(post_gush)/nc)*3)]

# Plot the chains
plot(chain1, col="blue", type='l', ylab=expression(theta),
     xlab='iteration')
lines(chain2, col="red")
lines(chain3, col='green')

# Estimates ----
# Now that we know things converged nicely, let's have a look at the results
# Plot a histogram of the posterior
hist(post_gush, main = '', xlab = expression(theta), ylab = '', yaxt='n',
     col='gray87')
abline(v=mean(post_gush), col='blue', lwd=2)
abline(v=quantile(post_gush, 0.025), col='red', lty=2, lwd=2)
abline(v=quantile(post_gush, 0.975), col='red', lty=2, lwd=2)

# Plot a density curve of the posterior
plot(density(post_gush), main='', xlab = expression(theta), col='blue')

# Dang, I was wrong again :/