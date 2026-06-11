library(tidyverse)
library(rstan)
library(Rcpp)
# remotes::install_github("LGhilotti/ProductFormFA")
# options(rgl.useNULL = TRUE)
library(ProductFormFA)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)
rm(list = ls())

source("functs/helper.R")
source("functs/rarefaction.R")

# Nsim = 25
Nsim = 3
beta = c(-3, 0.5, -0.3)
H = length(beta)
# Q_list = c(5, 8, 10, 15, 20, 25, 30)  
Q_list = c(2, 3)  
a_phi_true = 4
b_phi_true = 2
r0_true = 1000
nu0_true = 5000
b0 = c(-2, rep(0, H - 1))
B0 = diag(1, H)
r0 = 10

#-------------# model MOSAIC helper #---------------#

run_mosaic = function(sim, Q, base_seed = 123) {
  message("Running MOSAIC simulation ", sim, " / ", Nsim, "  [Q = ", Q, "]")
  sim_seed = base_seed + sim
  set.seed(sim_seed)
  n = sample(15:40, Q, replace = TRUE)
  
  data_sim = synt_dataH(n = n, Q = Q, beta = beta, a_phi = a_phi_true, 
                        b_phi = b_phi_true, r0 = r0_true, nu0 = nu0_true, 
                        mseed = sim_seed)
  tt_split = split_train_test(data = data_sim, train_prop = 0.8, mseed = sim_seed,
                              test_n = 10)
  accum_groups = accum_by_group(Y = data_sim$Y, n = data_sim$n, lab_group = 1:Q,
                                method = "random", permutations = 500, 
                                mseed = sim_seed)
  accum_test <- accum_groups %>%
    group_by(group) %>%
    filter(x > tt_split$n_train[group]) %>%
    ungroup()
  
  nu0 = 1.5*tt_split$k_train
  stan_data = list(Q = Q, H = H, k = tt_split$k_train, n_q = tt_split$n_train,
                    m = tt_split$m_train, X = data_sim$X, b0 = b0, B0 = B0,
                    a_phi = 2, b_phi = 2, r0 = r0, nu0 = nu0)
  # fit_sim = stan(file = "functs/MOSAIC_NegBin.stan", data = stan_data, seed = sim_seed,
  #                 chains = 4, iter = 2000, warmup = 1000)
  fit_sim = stan(file = "functs/MOSAIC_NegBin.stan", data = stan_data, seed = sim_seed,
                chains = 4, iter = 200, warmup = 100)
  post_param = rstan::extract(fit_sim, pars = c("alpha_out", "theta_out"))
  
  pred_df = do.call(rbind, lapply(1:Q, function(q) {
    g1 = rar_pred_Kgroup(alphaq = post_param$alpha_out[, q], 
                         thetaq = post_param$theta_out[, q],
                              r0 = r0, nu0 = nu0, nq = tt_split$n_train[q],
                              kq = tt_split$kq_train[q],
                              sq = tt_split$n_test[q], seed = sim_seed)
    true_q = accum_test[accum_test$group == q, ]
    data.frame(sim = sim, Q = Q, q = q, s = seq_len(tt_split$n_test[q]),
               mean = apply(g1, 1, mean), true = true_q$n_feat) %>%
      mutate(MSE = (mean - true)^2 / true)
  }))
  
  pred_df %>%
    group_by(sim, Q, s) %>%
    summarise(MOSAIC = mean(MSE, na.rm = TRUE), .groups = "drop")
}

#--------# model BB (ProductFormFA) helper #--------#
run_BB = function(sim, Q, base_seed = 123) {
  message("Running BB simulation ", sim, " / ", Nsim, "  [Q = ", Q, "]")
  sim_seed = base_seed + sim
  set.seed(sim_seed)
  n = sample(15:40, Q, replace = TRUE)
  
  data_sim = synt_dataH(n = n, Q = Q, beta = beta, a_phi = a_phi_true,
                        b_phi = b_phi_true, r0 = r0_true, nu0 = nu0_true,
                        mseed = sim_seed)
  tt_split = split_train_test(data = data_sim, train_prop = 0.8, mseed = sim_seed,
                              test_n = 10)
  accum_groups = accum_by_group(Y = data_sim$Y, n = data_sim$n, lab_group = 1:Q,
                                method = "random", permutations = 500,
                                mseed = sim_seed)
  accum_test = accum_groups %>%
    group_by(group) %>%
    filter(x > tt_split$n_train[group]) %>%
    ungroup()
  
  nu0 = 1.5*tt_split$k_train
  eb_init = list(alpha = -1, s = 1)
  eb_known = list(mu0 = nu0, var_fct = 10)
  eb_params_obj = eb_params(model = "NegBinBB", init = eb_init, known = eb_known)

  eb_df = do.call(rbind, lapply(1:Q, function(q) {
    eb_fit = GibbsFA_eb(feature_matrix = tt_split$Y_train[[q]],
                         model = "NegBinBB", type = "EFPF",
                         eb_params = eb_params_obj)
    true = accum_test %>% filter(group == q) %>% arrange(x) %>% pull(n_feat)
    mean_pred = unname(unlist(extrapolation(object = eb_fit,
                                             M = tt_split$n_test[q])$mu0_post))
    mean_pred = mean_pred + tt_split$kq_train[q]
    data.frame(sim = sim, Q = Q, q = q,
               s = seq_len(tt_split$n_test[q]),
               mean = mean_pred, true = true) %>%
      mutate(MSE = (mean - true)^2 / true)
  }))
  
  eb_df %>%
    group_by(sim, Q, s) %>%
    summarise(BB = mean(MSE, na.rm = TRUE), .groups = "drop")
}

#-----------------# run simulations #------------------#

res_MOSAIC = bind_rows(lapply(Q_list, function(Q) {
  bind_rows(lapply(1:Nsim, function(sim) run_mosaic(sim, Q))) %>%
    group_by(Q, s) %>%
    summarise(MOSAIC = mean(MOSAIC, na.rm = TRUE), .groups = "drop")
}))

res_BB = bind_rows(lapply(Q_list, function(Q) {
  bind_rows(lapply(1:Nsim, function(sim) run_BB(sim, Q))) %>%
    group_by(Q, s) %>%
    summarise(BB = mean(BB, na.rm = TRUE), .groups = "drop")
}))


# save(res_pBBI, file = "resDCC/simSpecific_QmulPBBI.RData")
# save(res_BB,   file = "resDCC/simSpecific_Qmul_BB.RData")

