rm(list = ls())
setwd("C:/Users/Marek/Documents/R_Scripts/Ekonometria_bayesowska/Praca_domowa2")
data <- read.csv("C:/Users/Marek/Documents/R_Scripts/Ekonometria_bayesowska/Praca_domowa2/ogorki.csv")

data <- data[, c("Q_cuc", "P_cuc","P_tom", "P_oni")]
data_log <- log(data)
ols <- lm(Q_cuc ~. , data = data_log)
summary(ols)

Q_cuc <- data_log$Q_cuc
P_cuc <- data_log$P_cuc
P_tom <- data_log$P_tom
P_oni <- data_log$P_oni
N <- nrow(data_log)

library(rstan)
library(bridgesampling)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

n_iter <- 1000
n_chains <- 4

ogorki_fit_full <- stan(file = "ogorki_full.stan",
                   data = c("N", "Q_cuc", "P_cuc", "P_tom", "P_oni"),
                   iter = n_iter,
                   chains = n_chains)
print(ogorki_fit_full)
stan_dens(ogorki_fit_full)

full_sampled <- bridge_sampler(ogorki_fit_full, silent = TRUE)

ogorki_fit_no_P_cuc <- stan(file = "ogorki_no_P_cuc.stan",
                   data = c("N", "Q_cuc", "P_tom", "P_oni"),
                   iter = n_iter,
                   chains = n_chains)

print(ogorki_fit_no_P_cuc)
stan_dens(ogorki_fit_no_P_cuc)

no_P_cuc_sampled <- bridge_sampler(ogorki_fit_no_P_cuc, silent = TRUE)

ogorki_fit_no_P_tom <- stan(file = "ogorki_no_P_tom.stan",
                   data = c("N", "Q_cuc", "P_cuc", "P_oni"),
                   iter = n_iter,
                   chains = n_chains)

print(ogorki_fit_no_P_tom)
stan_dens(ogorki_fit_no_P_tom)

no_P_tom_sampled <- bridge_sampler(ogorki_fit_no_P_tom, silent = TRUE)

ogorki_fit_no_P_oni <- stan(file = "ogorki_no_P_oni.stan",
                   data = c("N", "Q_cuc", "P_cuc", "P_tom"),
                   iter = n_iter,
                   chains = n_chains)

print(ogorki_fit_no_P_oni)
stan_dens(ogorki_fit_no_P_oni)

no_P_oni_sampled <- bridge_sampler(ogorki_fit_no_P_oni, silent = TRUE)



bf(full_sampled, no_P_cuc_sampled)
bf(full_sampled, no_P_tom_sampled)
bf(full_sampled, no_P_oni_sampled)


N_new <- 1

P_cuc_fct <- 3.321432
P_tom_fct <- 3.990734
P_oni_fct <- 2.921797 #obs nr 20

ogorki_fit_fct <- stan(file = "ogorki_full_prediction.stan",
                        data = c("N", "Q_cuc", "P_cuc", "P_tom", "P_oni",
                                "P_cuc_fct", "P_tom_fct", "P_oni_fct"),
                        iter = n_iter,
                        chains = n_chains)

print(ogorki_fit_fct)
stan_dens(ogorki_fit_fct, "Q_cuc_new")

