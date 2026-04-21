# Package load ----
library(tidyverse)
library(reshape)
library(R2jags)
library(RColorBrewer)
library(fields)
library(mgcv)

# Data ----
# White suckahs in the upper Susquehanna River 1987-2025
whs <- read.csv("data/whs.csv")

# Get sites with more than two samples 
# per year of monitoring and
counts <- whs %>% 
  filter(water != "Unnamed Water") %>% 
  group_by(water, year) %>% 
  summarize(n_sample = n()) %>% 
  filter(n_sample >= 3)

# Retain only sites with n >= 3 in sample data based
# on which waterbodies are included in the counts df above
whs <- whs %>% 
  filter(water %in% counts$water)

# Create a variable for replicates within year
# so we know how many times species could have been observed.
# We'll use this to collapse the Bernoulli likelihood
# to a true binomial so we can account for sparse design
# wherein not all sites were sampled equally (or at all) in
# all years
whs <- whs %>%
  group_by(water, year) %>%
  mutate(rep = row_number()) %>%
  ungroup()

# Data formatting ----
# . Observed presence/absence ----
# Create a wide-form capture history with sites
# in the rows and years in the columns
whs_occ_mat <- reshape::cast(data = whs,
                             formula = water ~ year,
                             value = "Detected",
                             drop = FALSE,
                             fill = NA, 
                             add.missing = TRUE,
                             fun.aggregate = max)

# . Number of surveys for likelihood ----
# Get number of surveys (reps) per year for
# each of the waterbodies 
k_mat <- whs %>% 
  group_by(water, year) %>% 
  summarize(reps = max(rep))

# Cast the previous into wide format so the structure
# matches the actual detection history with water in row
# and year in columns to match the capture matrix
k_mat <-reshape::cast(data = k_mat,
                      formula = water ~ year,
                      value = "reps",
                      drop = FALSE,
                      fill = NA, 
                      add.missing = TRUE,
                      fun.aggregate = mean)

# Just keep the data columns
k_mat <- k_mat[, 2:ncol(k_mat)]

# Mark NA as zero (reps) if the site was not surveyed in
# a given year
k_mat[is.na(k_mat)] <- 0

# . Gaussian radial basis function (spatial) ----
# Creates a 2-d distance-based grid that will allow
# parameters to vary spatially across a landscape

# Get coordinate pairs for each waterbody
sites <- whs %>% 
  select(water, nytme, nytmn) %>% 
  distinct() %>% 
  arrange(water)

# Cast coordinates as a matrix for cover design
coords <- as.matrix(sites[, c("nytme", "nytmn")])

# Knots - keeping it small to match dimensionality and extent
# of observations in the watershed. Fewer knots would be a 
# smoother surface and more knots would be a bumpier surface
# with more bends
knots <- fields::cover.design(coords, nd = 20)$design

# Scale and center coordinates for numerical stability
coords_sc <- scale(coords)
knots_sc  <- scale(
  knots,
  center = attr(coords_sc, "scaled:center"),
  scale  = attr(coords_sc, "scaled:scale")
)

# Create a distance matrix based on the scaled
# coordinates and the georeferenced knots (locations of
# potential bends in surface)
r <- rdist(coords_sc, knots_sc)

# Choose sigma based on knot spacing
d_knots <- rdist(knots_sc)
sigma <- mean(d_knots[lower.tri(d_knots)])

# Gaussian radial basis function (RBF)
X_smooth <- exp(-r^2 / (2 * sigma^2))

# Center the basis function for numerical stability
X_smooth <- scale(X_smooth, center = TRUE, scale = FALSE)

# Number of columns in basis
K_smooth = ncol(X_smooth)


# . Thin-plate spline basis function (temporal) ----
# Choose basis dimension and update iteratively as needed.
# This is the number of knots: same as above except now 
# we are working with a line instead of a surface. Start
# small
ky_smooth <- 5

# Build smooth object and remove intercept to reduce
# redundancy/colinearity in posteriors
sm <- smoothCon(
  s(year, bs = "tp", k = ky_smooth),
  data = data.frame(year = whs$year),
  absorb.cons = TRUE)[[1]]

# Get the basis from output
Xy_smooth <- sm$X

# Number of columns
Ky_smooth <- ncol(Xy_smooth)

# Center the basis for numerical stability like we did for the RBF
Xy_smooth <- scale(Xy_smooth, center = TRUE, scale = FALSE)


# . JAGS formatting ----
# Response (observed presence or absence)
y <- as.matrix(whs_occ_mat[, 2:ncol(whs_occ_mat)])

# .. Initial values ----
# NOTE: only needed if you retain the latent state z
# I marginalized it out of the likelihood in this 
# example to make the model run fast(er)
z_init <- y
z_init[z_init == 0] <- 1

  
# .. Package the data for jags ----
jags_data <- list(
  y = y,                              # Response
  n_site = length(unique(whs$water)), # Number of sites
  n_year = length(unique(whs$year)),  # Number of years
  K = k_mat,                          # Number of reps per site & year
  X_smooth = X_smooth,                # Spatial basis function
  K_smooth = K_smooth,                # Number of knots in spatial basis
  Xy_smooth = Xy_smooth,              # Temporal basis function
  Ky_smooth = Ky_smooth               # Number of knots in temp basis
)

# Function for initial values
inits <- function(){
  list(
    z = z_init
  )
}

# Parameters to save
# Just saving essentials for faster offloading
# these model results get biggg
params <- c("psi", "eps", "gamma", "p")

# Run the model ----
jags_fit <- jags(data = jags_data,
                 inits = NULL,
                 parameters.to.save = params,
                 model.file = "models/06_dynamic_occupancy_spatial_jags",
                 n.chains = 3,
                 n.iter = 500,
                 n.burnin = 250,
                 n.thin = 1)


# Results ----
# . Summary print ----
print(jags_fit, digits = 3)


# . Posterior of expected psi ----
# Extract MCMC list containing posteriors
posts <- jags_fit$BUGSoutput$sims.list


# .. Posterior summary ----
# Flatten the MCMC list into a dataframe for tidy workflows
psi_ests <- reshape::melt(posts$psi)
names(psi_ests) <- c("iteration", "site", "year", "estimate")

# Posterior summary for psi using our usual workflow
psi_summary <- psi_ests %>% 
  group_by(site, year) %>% 
  summarize(fit = mean(estimate),
            lwr = quantile(estimate, 0.025),
            upr = quantile(estimate, 0.975))


# .. Referencing posteriors to sites & years ----
# Get human-readable sites and years
psi_summary$water <- sort(unique(whs$water))[psi_summary$site]
psi_summary$year <- sort(unique(whs$year))[psi_summary$year]


psi_spatial_preds <- merge(psi_summary,
                           sites,
                           by = c("water"))


# .. Plot of spatio-temporal psi predictions ----
# Color ramp ofc
myPalette <- colorRampPalette(rev(brewer.pal(5, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(10), 
                             limits=c(min(psi_spatial_preds$fit),
                                      max(psi_spatial_preds$fit)))

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
          legend.direction = "horizontal",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~year) +
    labs(x = "Easting", y = "Northing",
         color = expression("E("~psi~")"))

