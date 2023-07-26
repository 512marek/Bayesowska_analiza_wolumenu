data {
  int<lower=0> N;
  vector[N] Q_cuc;
  vector[N] P_cuc;
  vector[N] P_tom;
  vector[N] P_oni;
}

parameters {
  real beta_P_cuc;
  real beta_P_tom;
  real beta_P_oni;
  real intercept;
  real<lower=0> sigma;
}


model {
  beta_P_cuc ~ logistic(-0.57, 0.15);
  beta_P_tom ~ normal(-0.04, 0.36);
  beta_P_oni ~ normal(-0.11, 0.36);
  intercept ~ uniform(-5, 25);
  sigma ~ inv_gamma(3, 1);
  Q_cuc ~ normal(intercept + beta_P_cuc * P_cuc + beta_P_tom * P_tom + beta_P_oni * P_oni, sigma);
}

