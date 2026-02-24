data {
  int<lower=0> N; // Number of observations
  int<lower=0, upper=1> y[N]; // Binary outcome variable
  int<lower=0> K; // Number of columns in design matrix
  matrix[N, K] X; // Predictor variables
}
parameters {
  vector[K] beta; // Coefficients
}
model {
  y ~ bernoulli_logit(X * beta); // Logistic regression likelihood
  beta ~ normal(0, 1); // Prior for coefficients
}
