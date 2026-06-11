library(tidyverse)
library(parallel)

rm(list = ls())

beta = c(-3, 0.5, -0.3)
H = length(beta)
Q = 15
nu0_true = 1000
a_phi_true = 4
b_phi_true = 2
mseed = 434
b0 = c(-2, rep(0, H - 1))
B0 = diag(1, H)
r0 = 10
a_phi = 2
b_phi = 2
nSIM = 500

n_chains = 4
n_parallel = 5   # 5 x 4 = 20 cores

set.seed(mseed)
n = sample(15:40, Q, replace = TRUE)
phi_true = rgamma(Q, a_phi_true, b_phi_true)

#--------# worker function #----------#
run_one = function(r, r0_true) {
  library(rstan)
  source("functs/helper.R")
  rstan_options(auto_write = TRUE)
  options(mc.cores = n_chains)
  
  data_sim = synt_dataH(n, Q, beta, phi_true, r0_true, nu0_true, mseed + r)
  nu0 = 1.5 * data_sim$k
  stan_data = list(Q = Q, H = H, k = data_sim$k, n_q = data_sim$n, m = data_sim$m,
                    X = data_sim$X, b0 = b0, B0 = B0, a_phi = a_phi, b_phi = b_phi,
                    r0 = r0, nu0 = nu0)
  
  fit_sim = stan(file = "functs/mosaic_NegBin.stan", data = stan_data, seed = mseed + r,
                  chains = n_chains, iter = 2000, warmup = 1000)

  post = rstan::extract(fit_sim)
  beta_lo = apply(post$beta, 2, quantile, probs = 0.025)
  beta_hi = apply(post$beta, 2, quantile, probs = 0.975)
  phi_lo = apply(post$phi,  2, quantile, probs = 0.025)
  phi_hi = apply(post$phi,  2, quantile, probs = 0.975)
  list(cover_beta = as.integer((beta >= beta_lo) & (beta <= beta_hi)),
       cover_phi  = as.integer((phi_true >= phi_lo) & (phi_true <= phi_hi)))
}

#------# generic runner #--------#
run_scenario = function(r0_true) {
  cl = makeCluster(n_parallel, type = "PSOCK")
  clusterExport(cl, varlist = c("n", "Q", "H", "beta", "phi_true",
                                "r0_true", "nu0_true", "b0", "B0",
                                "r0", "a_phi", "b_phi", "mseed", "n_chains"),
                envir = environment())
  res_list = parLapply(cl, 1:nSIM, run_one, r0_true = r0_true)
  stopCluster(cl)
  
  cover_beta = do.call(rbind, lapply(res_list, `[[`, "cover_beta"))
  cover_phi = do.call(rbind, lapply(res_list, `[[`, "cover_phi"))
  list(cover_beta = cover_beta, cover_phi = cover_phi)
}

# run both scenarios
results_5000 = run_scenario(r0_true = 5000)
results_3000 = run_scenario(r0_true = 3000)


