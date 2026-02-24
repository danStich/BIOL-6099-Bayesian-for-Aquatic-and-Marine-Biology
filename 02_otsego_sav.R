# Library loads ----
library(tidyverse)
library(rstanarm)
library(rstan)
library(sf)
library(RColorBrewer)

# Data ----
# . Read in the data ----
# Sample data
plants <- read.csv("data/plant_abundances.csv")

# . Data manipulation ----
# Stack data into long form
sav_data <- pivot_longer(data = plants,
                          cols = c(8:ncol(plants)),
                          names_to = "species",
                          values_to = "abundance")

# . Presence/absence ----
# Create a new column for presence/absence 
# from the "abundance" column initializing it as all 1s
sav_data$detected <- 1

# Now, mark detected zero wherever abundance was Z for zero
sav_data$detected[sav_data$abundance == "Z"] <- 0


# . Survey variable formatting ----
sav_data$survey <- factor(sav_data$Survey,
                          labels = c("June", "July", "August"))


# . Data subsetting ----
sav_data <- sav_data %>%
  filter(species == "Vallisneria_americana")


# . Standardizing continuous covariates ----
# Convergence is good, but the coefficient estimates are huge
# if we don't standardize covariates because the scale of 
# Latitude and Longitude is far away from zero. In practice, 
# it is best to "scale and center", or "standardize" 
# covariates like this in Bayesian, and in general really...
sav_data$s_latitude <- scale(sav_data$Latitude)
sav_data$s_longitude <- scale(sav_data$Longitude)

# Model fitting with rstanarm ----
# . Model fit ----
# Then we can re-run it just like this:
sav_mod <- stan_glm(detected ~ survey + s_latitude * s_longitude, 
                    data = sav_data,
                    family = binomial,
                    chains = 3,
                    iter = 2000,
                    warmup = 1000,
                    cores = 3
)

# You can get the model code like this
# rstan::get_stanmodel(sav_mod$stanfit)


# You can look at the summary like this
summary(sav_mod)

# . Predictions ----
# The major convenience of rstanarm is that it pipes directly into our
# usual linear modeling workflow in R. Make predictions
sav_mod_preds <- predict(sav_mod, type = "response")

# Then smoosh it together with the raw data
sav_preds <- data.frame(sav_data, 
                        fit = sav_mod_preds)

# . Plotting results ----
# .. A color palette for sweet, sweet predictions ----
myPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")),
                              alpha = 0.5)
sc <- scale_colour_gradientn(colours = myPalette(10), 
                             limits=c(0, 1),
                             oob = scales::squish)

# .. A fittingly awesome graph ----
ggplot(sav_preds, aes(x = Longitude, y = Latitude, color = fit)) +
  geom_point(size = 3) +
  coord_sf() +
  sc +
  facet_wrap(~survey)



# Model fitting with rstan ----
# . Covariates ----
X <- model.matrix(lm(detected ~ survey + s_latitude * s_longitude,
                     data = sav_data))


# . Packaging data for rstan ----
data_list <- list(
  N = nrow(sav_data),   # Number of rows in data
  X = X,                # Model matrix with covariates
  K = ncol(X),          # Number of model parameters (betas)
  y = sav_data$detected # observed response
)

# . Model fit ----
fit <- stan(
  file = "models/otsego_sav_logistic.stan",
  data = data_list,
  iter = 2000,
  chains = 3
)

# . Model summary ----
# This is different from rstanarm
print(fit)

# . Predictions ----
posts <- extract(fit)

beta_fit <- apply(posts$beta, 2, mean)
beta_lwr <- apply(posts$beta, 2, quantile, 0.025)
beta_upr <- apply(posts$beta, 2, quantile, 0.025)

sav_preds_stan <- sav_data
sav_preds_stan$fit <- invlogit(X %*% beta_fit)
sav_preds_stan$lwr <- invlogit(X %*% beta_lwr)
sav_preds_stan$upr <- invlogit(X %*% beta_upr)

# .. A color palette for sweet, sweet predictions ----
myPalette <- colorRampPalette(rev(brewer.pal(10, "Spectral")),
                              alpha = 0.5)
sc <- scale_colour_gradientn(colours = myPalette(10), 
                             limits=c(0, 1),
                             oob = scales::squish)

# .. A fittingly awesome graph ----
ggplot(sav_preds_stan, aes(x = Longitude, y = Latitude, color = fit)) +
  geom_point(size = 3) +
  coord_sf() +
  sc +
  facet_wrap(~survey)





