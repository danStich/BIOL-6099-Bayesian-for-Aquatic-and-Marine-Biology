# Front-end needs -----
# Load necessary packages
  library(R2jags)

# Make some convenience functions ----
# Create function to invert logit link function
  inv.logit = function(x){
    exp(x)/(1+exp(x))
  }

# Make a function to get lower 95% credible limit with short name
  low = function(x){
    quantile(x, probs=c(0.025))
  }

# Make a function to get upper 95% credible limit with short name
  up = function(x){
    quantile(x, probs=c(0.975))
  }

# Data -----
# Read in the data file
# This file contains data that were collected to settle an argument over
# whether one side of our "Scramble" board game set was luckier than the
# other after allegations of "cheating" were thrown around. The values
# column are the number of pieces that popped out of the board at the end
# of the timer, ind is the color of the pieces, and total is the total 
# starting number of pieces on each side at the end of the timer (we got 
# really good at it, and then really bored of it). Google the game if
# you're lost.
  scramble <- read.csv('data/projectScramble.csv')
  scram <- data.frame(stack(scramble, c('green', 'pink')), total=18)

# Model -----
# Binomial logistic regression model.
# This shows an example of how you can just save your model to a 
# string instead of a separate file. It also shows how we can use
# "derived" quantities to make inference (difference between colors).
  modelString = "
    model{
      # Likelihood
        for(i in 1:N){
          Y[i] ~ dbinom(p[color[i]], n[i])
        }

      # Priors
      for(j in 1:ncolors){
        p[j] ~ dbeta(1, 1)
      }

      # Difference between colors
      diff <- p[1] - p[2]
      
    }"

# Package the data for JAGS
  scramble_data = list(
    Y = scram$values,
    n = scram$total,
    N = nrow(scram),
    color=as.numeric(as.factor(scram$ind)),
    ncolors=length(unique(scram$ind))
  )

# Parameters monitored
  params = c('p', 'diff')

# Initial values for parameters
  inits <- function(){
    list(
      p = rbeta(2, 1, 1)
    )
  }

# MCMC settings
  ni <- 50000      # Number of draws from posterior (for each chain)
  nt <- 5          # Thinning rate
  nb <- 25000      # Number of draws to discard as burn-in
  nc <- 3          # Number of chains

# Call jags and run the model
  scramble_mod <- jags(data = scramble_data, inits=inits, params,
                 textConnection(modelString),
                 n.chains = nc, n.thin = nt,
                 n.iter = ni, n.burnin = nb,
                 working.directory = getwd())
  
  print(scramble_mod)

# Results ----  
  difference <- scramble_mod$BUGSoutput$sims.list$diff

# Graph the difference between green and pink sides of the board. If
# zero is within the 95% uncertainty interval (red lines) then the
# game set is not biased (even if certain people *were* cheating)
  par(mar=c(4,4,4,1))
  hist(difference, col='gray87', xlim=c(-.2, .2), axes=F,
       ylim = c(0, 3000), main='Project scramble',
       xlab='Difference (green-pink)',
       border = 'gray90')
  axis(1, pos=0)
  axis(2, las=2, pos=-.2)
  abline(v=mean(difference), col='blue', lty=1, lwd=2)
  abline(v=quantile(difference, c(0.025, 0.975)),
         col='red', lty=2, lwd=2)
