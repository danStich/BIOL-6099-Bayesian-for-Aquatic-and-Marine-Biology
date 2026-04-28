library(nimble)

dyn_occ <- nimbleCode({
  
  #### ---------------------------
  #### Priors
  #### ---------------------------
  
  # Linear predictors of colonization and extinction
  for (i in 1:n_site) {
    for (t in 1:n_year) {
      
      # Colonization
      logit(gamma[i, t]) <- beta0_gamma +
        inprod(X_smooth[i, ], beta_gamma[1:K_smooth]) +
        inprod(Xy_smooth[t, ], alpha_gamma[1:Ky_smooth])
      
      # Extinction
      logit(eps[i, t]) <- beta0_eps +
        inprod(X_smooth[i, ], beta_eps[1:K_smooth]) +
        inprod(Xy_smooth[t, ], alpha_eps[1:Ky_smooth])
    }
  }
  
  # Priors on spatial parameters (non-centered!)
  for (k in 1:K_smooth) {
    
    # Occupancy
    beta_psi[k] <- sigma_psi * z_psi[k]
    z_psi[k] ~ dnorm(0, 1)
    
    # Colonization
    beta_gamma[k] <- sigma_gamma * z_gamma[k]
    z_gamma[k] ~ dnorm(0, 1)
    
    # Extinction
    beta_eps[k] <- sigma_eps * z_eps[k]
    z_eps[k] ~ dnorm(0, 1)
  }
  
  # Priors for intercepts
  beta0_psi ~ dnorm(0, 1)
  beta0_gamma ~ dnorm(0, 1)
  beta0_eps ~ dnorm(0, 1)
  
  # Priors for variance components
  sigma_psi ~ dunif(0, 5)
  sigma_gamma ~ dunif(0, 5)
  sigma_eps ~ dunif(0, 5)
  
  # Temporal smooth terms
  for (t in 1:Ky_smooth) {
    alpha_gamma[t] ~ dnorm(0, 1)
    alpha_eps[t] ~ dnorm(0, 1)
  }
  
  # Detection
  logit(p) <- p_lp
  p_lp ~ dnorm(0, 1)
  
  
  #### ---------------------------
  #### State process (marginalized)
  #### ---------------------------
  
  for (i in 1:n_site) {
    
    # Occupancy in first year
    logit(psi[i, 1]) <- beta0_psi +
      inprod(X_smooth[i, ], beta_psi[1:K_smooth])
    
    # Dynamic occupancy recursion
    for (t in 2:n_year) {
      
      psi[i, t] <- psi[i, t - 1] * (1 - eps[i, t]) +
        (1 - psi[i, t - 1]) * gamma[i, t]
    }
  }
  
  
  #### ---------------------------
  #### Observation model
  #### ---------------------------
  
  for (i in 1:n_site) {
    for (t in 1:n_year) {
      y[i, t] ~ dbin(p * psi[i, t], K[i, t])
    }
  }
  
})
