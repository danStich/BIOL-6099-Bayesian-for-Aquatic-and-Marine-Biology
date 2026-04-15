# Package load ----
library(tidyverse)
library(reshape)
library(R2jags)
library(RColorBrewer)

# Data ----
whs <- read.csv("data/whs.csv")

# Get sites with more than two samples 
# per year of monitoring
counts <- whs %>% 
  group_by(water, year) %>% 
  summarize(n_sample = n()) %>% 
  filter(n_sample >= 3, water != "Unnamed Water")

# Retain only sites with n >= 3 in sample data
whs <- whs %>% 
  filter(water %in% counts$water)

# Create a variable for replicates within year
whs <- whs %>%
  group_by(water, year) %>%
  mutate(rep = row_number()) %>%
  ungroup()

# Data formatting ----
# Observed presence/absence
# Create a wide-form capture history
whs_occ_mat <- reshape::cast(data = whs,
                             formula = water ~ year,
                             value = "Detected",
                             drop = FALSE,
                             fill = NA, 
                             add.missing = TRUE,
                             fun.aggregate = max)

# Number of surveys for m-array likelihood
k_mat <- whs %>% 
  group_by(water, year) %>% 
  summarize(reps = max(rep))

k_mat <-reshape::cast(data = k_mat,
                      formula = water ~ year,
                      value = "reps",
                      drop = FALSE,
                      fill = NA, 
                      add.missing = TRUE,
                      fun.aggregate = mean)

k_mat <- k_mat[, 2:ncol(k_mat)]
k_mat[is.na(k_mat)] <- 0



# Response
y <- as.matrix(whs_occ_mat[, 2:ncol(whs_occ_mat)])

# Initial values for occupancy state
z_init <- y
z_init[z_init == 0] <- 1

  
# Package the data for jags
jags_data <- list(
  y = y,
  n_site = length(unique(whs$water)),
  n_year = length(unique(whs$year)),
  K = k_mat
)

# Function for initial values
inits <- function(){
  list(
    z = z_init
  )
}

# Parameters to save
params <- c("psi", "eps", "gamma", "p")

# Run the model ----
jags_fit <- jags(data = jags_data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "models/05_dynamic_occupancy_jags",
                 n.chains = 3,
                 n.iter = 5000,
                 n.burnin = 2500,
                 n.thin = 1)

print(jags_fit, digits = 3)


# For Class ----
posts <- jags_fit$BUGSoutput$sims.list

names(posts)

psi_ests <- reshape::melt(posts$psi)
names(psi_ests) <- c("iteration", "site", "year", "estimate")

psi_summary <- psi_ests %>% 
  group_by(site, year) %>% 
  summarize(fit = mean(estimate),
            lwr = quantile(estimate, 0.025),
            upr = quantile(estimate, 0.975))


# Get human-readable sites and years
psi_summary$water <- sort(unique(whs$water))[psi_summary$site]
psi_summary$year <- sort(unique(whs$year))[psi_summary$year]


# Get coordinates for site-specific predictions
sites <- whs[!duplicated(whs$water), c("water", "nytme", "nytmn")]

psi_spatial_preds <- merge(psi_summary,
                           sites,
                           by = c("water"))


# Color ramp
myPalette <- colorRampPalette(rev(brewer.pal(5, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(5), limits=c(0, 1))

# Plot
psi_spatial_preds %>% 
  arrange(fit) %>% 
  filter(year %in% c(1987, 1992, 1997, 2002, 2007, 
                     2012, 2017, 2022, 2025)) %>%
  ggplot(aes(x = nytme, y = nytmn, z = fit, color = fit)) +
    geom_point(size = 2, alpha = 0.5) +
    coord_sf() +
    sc + 
    theme_bw() +
    theme(legend.position = "top",
          legend.direction = "horizontal") +
    facet_wrap(~year)

