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
whs_occ <- whs %>% 
  filter(water %in% counts$water)

# Create a variable for replicates within year
whs <- whs %>%
  group_by(water, year) %>%
  mutate(rep = row_number()) %>%
  ungroup()

whs$site_year <- as.numeric( as.factor(
  paste0(whs$water, whs$year)
))


# Data formatting ----
# Observed presence/absence
y <- whs_occ$Detected

# Initial values for occupancy state
z_init <- y
z_init[z_init == 0] <- 1

# Package the data for jags
jags_data <- list(
  y = y,
  n_site = length(unique(whs_occ$water)),
  site = as.numeric(as.factor(whs_occ$water)),
  n_year = length(unique(whs_occ$year)),
  year = as.numeric(as.factor(whs_occ$year)),
  n_obs = nrow(whs_occ),
  site_year = whs$site_year
)

# Function for initial values
inits <- function(){
  list(
    z = z_init
  )
}

# Parameters to save
params <- c("psi", "p")

# Run the model ----
jags_fit <- jags(data = jags_data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "models/05_single_species_occ_temporal_jags",
                 n.chains = 3,
                 n.iter = 1000,
                 n.burnin = 500,
                 n.thin = 10)

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
psi_summary$water <- sort(unique(whs_occ$water))[psi_summary$site]
psi_summary$year <- sort(unique(whs_occ$year))[psi_summary$year]


# Get coordinates for site-specific predictions
sites <- whs_occ[!duplicated(whs_occ$water), c("water", "nytme", "nytmn")]

psi_spatial_preds <- merge(psi_summary,
                           sites,
                           by = c("water"))


# Color ramp
myPalette <- colorRampPalette(rev(brewer.pal(5, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(5), limits=c(0.10, 1))

ggplot(psi_spatial_preds,
       aes(x = nytme, y = nytmn, z = fit, color = fit)) +
  geom_point(size = 4, alpha = 0.25) +
  coord_sf() +
  sc + 
  theme_bw() +
  theme(legend.position = "top",
        legend.direction = "horizontal") +
  facet_wrap(~year)

