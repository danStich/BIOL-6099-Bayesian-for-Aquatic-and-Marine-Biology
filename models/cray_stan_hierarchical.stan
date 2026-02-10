// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
  int<lower=0> pop[N];
  int<lower=0> n_pop;
  real p_a;
  real p_b;
}

// The parameters accepted by the model.
parameters {
  vector[n_pop] alpha;
  real alpha_hyper;
  vector[n_pop] beta;
  real beta_hyper;
  real<lower=0> sigma;
}

// The model to be estimated. We model the response
// 'y' to be normally distributed with mean 'y_hat'
// and standard deviation 'sigma'.
model {
  vector[N] y_hat;
  y_hat = alpha[pop] + beta[pop] .* x;
  y ~ normal(y_hat, sigma);
  alpha ~ normal(alpha_hyper, 1);
  alpha_hyper ~ normal(p_a, 1);
  beta ~ normal(beta_hyper, 1);
  beta_hyper ~ normal(p_b, 1);
}

