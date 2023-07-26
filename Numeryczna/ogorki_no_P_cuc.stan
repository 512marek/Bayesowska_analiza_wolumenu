data {
  int<lower=0> N;
  vector[N] Q_cuc;
  vector[N] P_tom;
  vector[N] P_oni;
}

parameters {
  real beta_P_tom;
  real beta_P_oni;
  real intercept;
  real<lower=0> sigma;
}


model {
  beta_P_tom ~ normal(-0.04, 0.36);
  beta_P_oni ~ normal(-0.11, 0.36);
  intercept ~ uniform(-25, 25);
  sigma ~ inv_gamma(3, 1);
  Q_cuc ~ normal(intercept + beta_P_tom * P_tom + beta_P_oni * P_oni, sigma);
}

