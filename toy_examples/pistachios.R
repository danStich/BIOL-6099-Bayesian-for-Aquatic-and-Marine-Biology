# Library load
library(tidyverse)
library(rstanarm)

# Stan options
options(mc.cores = parallel::detectCores())

# Data read
# replicate: handful (sample)
# shelled: number of nuts successfully opened
# total: number of nuts total
pistachios <- read.csv("data/pistachios.csv")

# Data manipulation: calculate number of failures per handful
pistachios$failed <- pistachios$total - pistachios$shelled

# Run the model using successes and failures as binomial response.
# Here, we are just estimating the "mean" probability of shelling
# a given nut...should have kept track of who was opening them...
pistachio_model <- stan_glm(formula = cbind(shelled, failed) ~ 1,
         family = binomial(link = logit),
         data = pistachios
         )

# Have a look at these wild results
summary(pistachio_model, digits = 3)

# Get posterior probability that a nut can be easily shelled
# Postachio...haha
postachio_on_logit <- data.frame(pistachio_model)
names(postachio_on_logit) <- "p_success"

# Invert the logit
postachio <- postachio_on_logit %>%
  mutate(p_success = invlogit(p_success))

# Make a histogram
hist(postachio$p_success)

# Calculate mean and credible interval
mean(postachio$p_success)
quantile(postachio$p_success, c(0.025, 0.975))

# In the average handful, about 15% of the pistachios (11-21%)
# cannot be easily shelled. But, it looks like we'd want to collect
# more data to zero this in a bit
ggplot(postachio, aes(p_success)) + 
  geom_density() +
  geom_dotplot(mapping = aes(x = proportion),
               data = pistachios
               ) +
  xlab("Probability that a pistachio can be shelled") +
  ylab("Frequency") +
  scale_y_continuous(limits = c(0, 25), expand = c(0.01, 0)) +
  ggtitle("Posterior distribution after n = 15 samples") +
  theme_bw() +
  theme(
    axis.title.x = element_text(vjust = -1),
    axis.title.y = element_text(vjust = 3)
  )
