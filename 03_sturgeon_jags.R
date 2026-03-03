# Libraries ----
library(tidyverse)
library(lubridate)
library(reshape)
library(R2jags)


# Data ----
# . Data read ----
# Read data file to reproduce the analysis
juv_haverstraw <- read.csv("data/juv_haverstraw.csv")

juv_haverstraw <- juv_haverstraw %>%
  filter(!is.na(BATCH))

# . Covariates for model ----
juv_haverstraw$year <- as.character(juv_haverstraw$YEAR) 
juv_haverstraw$river_mile <- as.character(juv_haverstraw$RM)

X <- model.matrix(
  lm(NUMBER ~ year + river_mile, data = juv_haverstraw)
)

# Analysis ----
# . Jags data ----
jags_data <- list(
  y = juv_haverstraw$NUMBER,
  nobs = nrow(juv_haverstraw),
  X = X,
  K = ncol(X)
)


# . Parameters monitored----
params <- c("lambda", "beta")

# . Initial values ----
inits = function(){
  list(
    beta = rnorm(ncol(X), 0, 1)
  )
}  

# . Compile model ----
fit <- jags(jags_data, inits = inits,
            parameters.to.save = params,
            model.file = "models/juv_sturgeon_poisson",
            n.chains = 3, n.iter = 250,
            n.burnin = 100, n.thin = 10)

# Results ----
# . Print summary ----
print(fit, digits=3)

# . Extract posteriors ----
posts <- fit$BUGSoutput$sims.list

# Have a quick look
glimpse(posts)


# . Abundance estimates ----
preds <- data.frame(
  fit = apply(posts$lambda, 2, mean),
  lwr = apply(posts$lambda, 2, quantile, 0.025),
  upr = apply(posts$lambda, 2, quantile, 0.975)
)

# . Put predictions together with raw data ----
juv_preds <- data.frame(juv_haverstraw, preds)


# . Plot it out ----
ggplot(juv_preds, aes(x = YEAR, y = fit, color = river_mile,
                      fill = river_mile)) +
  geom_line() +
  geom_ribbon(aes(xmin = YEAR, 
                  ymin = lwr, 
                  ymax = upr, 
                  color = NULL), 
              alpha = 0.5) +
  ylab("Fish per net") +
  stat_summary(
    aes(x = YEAR, y = NUMBER),
    geom = "point",
    fun = "mean",
    col = "black",
    size = 3,
    inherit.aes = FALSE
  )
