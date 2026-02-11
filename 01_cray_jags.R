# Package load ----
library(R2jags)
# library(JagsUI)
library(tidyverse)
library(lubridate)

# Data manipulation ----
# . Read data ----
cray <-  read.csv("data/cray.csv", stringsAsFactors = FALSE)
cpues <- read.csv('data/cpue.csv', stringsAsFactors = FALSE)
cpues$waterbody_site <- paste0(cpues$waterbody, " - ", cpues$site)
  
# Merge the data sets 
cray_data <- merge(cray, cpues, by=c('date', 'waterbody', 'site'))
nrow(cray_data)
  
# Log transform length and mass
cray_data$loglength <- log10(cray$length)
cray_data$logmass <- log10(cray$mass)
  
# Get rid of NA values
cray_data <- cray_data %>% 
  filter(!is.na(loglength) & !is.na(logmass))
  
# Transform date and get year
cray_data$date <- as.Date(as.character(cray_data$date), format = "%m/%d/%Y")
cray_data$year <- year(cray_data$date)
cray_data$days <- yday(cray_data$date)

# Add date to grouping variable
cray_data$waterbody_site_year <- paste0(cray_data$waterbody_site, " - ",
                                        cray_data$year)

# Bundle data
mod_data <- list(
  y = cray_data$logmass,
  x = cray_data$loglength,    
  pop = as.numeric(as.factor(cray_data$waterbody_site)), 
  N = nrow(cray_data),
  n_pop = length(unique(cray_data$waterbody_site))
  )



# Model calibration ----
# Parameters to estimate
params <- c("alpha", "beta", "sigma")

# Fit the model
fit <- jags(data = mod_data,
            parameters.to.save = params,
            model.file = 'models/cray_jags_hierarchical',
            n.chains = 3,
            n.iter = 5000,
            n.burnin = 2500,
            n.thin = 5)

# Result ----
print(fit)

# Extract posteriors
posts <- fit$BUGSoutput$sims.list

alpha_fit <- apply(posts$alpha, 2, median)
beta_fit <- apply(posts$beta, 2, median)

# Making predictions ----
# Make a new sequence of log10 lengths
new_lengths <- seq(min(cray_data$loglength), 
                   max(cray_data$loglength),
                   0.01)

# Create an empty matrix to hold preds for each population
preds <- matrix(nrow = length(new_lengths), ncol = length(alpha_fit))

# For each new length value (rows) in each population (columns)
# preds is the mean outcome of the linear predictor (y_hat from model)
for(i in 1:nrow(preds)){
  for(t in 1:ncol(preds)){
    preds[i, t] <- alpha_fit[t] + beta_fit[t] * new_lengths[i]
  }
}

# Re-organize predictions into a data frame and tidy it up
means <- data.frame(preds)
means$loglength <- new_lengths

plotter <- pivot_longer(means, cols = 1:ncol(preds),
                        names_to = "pop", 
                        values_to = "logmass",
                        names_prefix = "X") %>% 
  mutate(pop = as.numeric(pop))

plotter$waterbody_site <- sort(unique(cray_data$waterbody_site))[plotter$pop]
plotter$length <- 10^plotter$loglength
plotter$mass <- 10^plotter$logmass


# Plot predictions against raw data ----
ggplot(plotter, aes(x = length, y = mass, color = waterbody_site)) +
  geom_line() +
  geom_point(data = cray_data, aes(x = length, y = mass)) +
  facet_wrap(~waterbody_site)





