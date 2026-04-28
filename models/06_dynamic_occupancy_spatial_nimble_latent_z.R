library(nimble)

dyn_occ_latent <- nimbleCode({
  
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
  
  # Spatial smooth coefficients (non-centered)
  for (k in 1:K_smooth) {
    
    beta_psi[k]   <- sigma_psi   * z_psi[k]
    beta_gamma[k] <- sigma_gamma * z_gamma[k]
    beta_eps[k]   <- sigma_eps   * z_eps[k]
    
    z_psi[k]   ~ dnorm(0, 1)
    z_gamma[k] ~ dnorm(0, 1)
    z_eps[k]   ~ dnorm(0, 1)
  }
  
  # Intercepts
  beta0_psi   ~ dnorm(0, 1)
  beta0_gamma ~ dnorm(0, 1)
  beta0_eps   ~ dnorm(0, 1)
  
  # Variance components
  sigma_psi   ~ dunif(0, 5)
  sigma_gamma ~ dunif(0, 5)
  sigma_eps   ~ dunif(0, 5)
  
  # Temporal smooths
  for (t in 1:Ky_smooth) {
    alpha_gamma[t] ~ dnorm(0, 1)
    alpha_eps[t]   ~ dnorm(0, 1)
  }
  
  # Detection
  logit(p) <- p_lp
  p_lp ~ dnorm(0, 1)
  
  
  #### ---------------------------
  #### State process (latent z)
  #### ---------------------------
  
  for (i in 1:n_site) {
    
    # Initial occupancy
    logit(psi[i, 1]) <- beta0_psi +
      inprod(X_smooth[i, ], beta_psi[1:K_smooth])
    z[i, 1] ~ dbern(psi[i, 1])
    
    # Dynamic occupancy
    for (t in 2:n_year) {
      
      muZ[i, t] <- z[i, t-1] * (1 - eps[i, t]) +
        (1 - z[i, t-1]) * gamma[i, t]
      
      z[i, t] ~ dbern(muZ[i, t])
    }
  }
  
  
  #### ---------------------------
  #### Observation process
  #### ---------------------------
  
  for (i in 1:n_site) {
    for (t in 1:n_year) {
      y[i, t] ~ dbin(z[i, t] * p, K[i, t])
    }
  }

  
})