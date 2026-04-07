# Package load ----
library(tidyverse)
library(reshape)
library(R2jags)


# Data ----
whs <- read.csv("data/whs.csv")

# Get sites with more than two samples
counts <- whs %>% 
  group_by(water) %>% 
  summarize(n_sample = n()) %>% 
  filter(n_sample >= 3)

# Retain only sites with n >= 3 in sample data
whs_occ <- whs %>% 
  filter(water %in% counts$water)

# Create a wide-form capture history
whs_occ_mat <- reshape::cast(data = whs_occ,
                        formula = water ~ year,
                        value = "Detected",
                        drop = FALSE,
                        fill = NA, 
                        add.missing = TRUE,
                        fun.aggregate = max)

# Extract only the capture history columns
# bc first column has waterbody name
whs_ch <- whs_occ_mat[, 2:ncol(whs_occ_mat)]


# Data formatting ----
# Capture matrix
y <- as.matrix(whs_ch)

# Initial values for occupancy state
z_init <- apply(y, 1, max, na.rm = TRUE)
z_init[z_init == 0] <- 1

# Package the data for jags
jags_data <- list(
  y = y,
  n_site = nrow(y),
  n_year = ncol(y)
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
                 model.file = "models/04_single_species_occ_jags.jags",
                 n.chains = 3,
                 n.iter = 1000,
                 n.burnin = 500,
                 n.thin = 10)

print(jags_fit, digits = 3)


# For Class ----
# Color ramp
myPalette <- colorRampPalette(rev(brewer.pal(5, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(10), limits=c(10, 25))

