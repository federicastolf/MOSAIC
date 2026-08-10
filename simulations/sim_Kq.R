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

#--------------------# fixed settings (shared) #--------------------#
Nsim = 25
beta = c(-3, 0.5, -0.3)
H = length(beta)
Q = 15
a_phi_true = 4
b_phi_true = 2
r0_true = 1000
nu0_true = 5000
b0 = c(-2, rep(0, H - 1))
B0 = diag(1, H)
r0 = 10

#-------# data-generation helper #-----------#
make_data = function(n, Q, beta, a_phi_true, b_phi_true, r0_true, nu0_true, 
                      sim_seed, misspecified) {
  if (misspecified) {
    synt_dataH(n = n, Q = Q, beta = beta, a_phi = a_phi_true, b_phi = b_phi_true, 
               r0 = r0_true, nu0 = nu0_true, mseed = sim_seed, tau = 0.1)
  } else {
    synt_dataH(n = n, Q = Q, beta = beta, a_phi = a_phi_true, b_phi = b_phi_true,
               r0 = r0_true, nu0 = nu0_true, mseed = sim_seed)
  }
}

#-----------# MOSAIC model helper #-------------#

run_mosaic = function(sim, misspecified = FALSE, base_seed = 123) {
  scenario_tag = if (misspecified) "B (misp)" else "A"
  message("Running MOSAIC simulation ", sim, " / ", Nsim, " 
          [Scenario ", scenario_tag, "]")
  
  sim_seed = base_seed + sim
  set.seed(sim_seed)
  n = sample(15:40, Q, replace = TRUE)
  data_sim = make_data(n, Q, beta, a_phi_true, b_phi_true, r0_true, nu0_true, 
                       sim_seed, misspecified)
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
  stan_data = list(Q = Q, H = H, k = tt_split$k_train, n_q = tt_split$n_train,
                    m = tt_split$m_train, X = data_sim$X, b0 = b0, B0 = B0,
                    a_phi = 2, b_phi = 2, r0 = r0, nu0 = nu0)
  fit_sim = stan(file = "functs/mosaic_NegBin.stan", data = stan_data, seed = sim_seed,
                    chains = 4, iter = 2000, warmup = 1000)
  post_param = rstan::extract(fit_sim, pars = c("alpha_out", "theta_out"))
  
  # p_0(alpha, theta, n) = prod_q p_{0,q}(alpha_q, theta_q, n_q)
  p0_all = Reduce(`*`, lapply(1:Q, function(q)
    d_n(n_q = tt_split$n_train[q], a_q = -post_param$alpha_out[, q],
        b_q = post_param$alpha_out[, q] + post_param$theta_out[, q])))
  
  
  pred_df = do.call(rbind, lapply(1:Q, function(q) {
    g1 = rar_pred_Kgroup(alphaq = post_param$alpha_out[, q], 
                         thetaq = post_param$theta_out[, q], r0 = r0, nu0 = nu0, 
                         nq = tt_split$n_train[q], kq = tt_split$kq_train[q], 
                         sq = tt_split$n_test[q], seed = sim_seed,
                         p0_all = p0_all, k = tt_split$k_train)
    true_q = accum_test[accum_test$group == q, ]
    data.frame(sim = sim, q = q, s = seq_len(tt_split$n_test[q]),
               mean = apply(g1, 1, mean), true = true_q$n_feat) %>%
      mutate(MSE = (mean - true)^2 / true)
  }))
  
  pred_df %>%
    group_by(sim, s) %>%
    summarise(MOSAIC = mean(MSE, na.rm = TRUE), .groups = "drop")
}

#---------# model BB (ProductFormFA)  helper #---------#

run_BB = function(sim, misspecified = FALSE, base_seed = 123) {
  scenario_tag <- if (misspecified) "B (misp)" else "A"
  message("Running BB simulation ", sim, " / ", Nsim, "  [Scenario ", scenario_tag, "]")
  
  sim_seed = base_seed + sim
  set.seed(sim_seed)
  n = sample(15:40, Q, replace = TRUE)
  
  data_sim = make_data(n, Q, beta, a_phi_true, b_phi_true, r0_true, nu0_true, 
                       sim_seed, misspecified)
  tt_split = split_train_test(data = data_sim, train_prop = 0.8, mseed = sim_seed,
                              test_n = 10)
  accum_groups = accum_by_group(Y = data_sim$Y, n = data_sim$n, lab_group = 1:Q, 
                                method = "random", permutations = 500, mseed = sim_seed)
  
  accum_test <- accum_groups %>%
    group_by(group) %>%
    filter(x > tt_split$n_train[group]) %>%
    ungroup()
  
  nu0 = 1.5 * tt_split$k_train
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
    data.frame(sim = sim, q = q, s = seq_len(tt_split$n_test[q]),
               mean = mean_pred, true = true) %>%
      mutate(MSE = (mean - true)^2 / true)
  }))
  
  eb_df %>%
    group_by(sim, s) %>%
    summarise(BB = mean(MSE, na.rm = TRUE), .groups = "drop")
}

#---------------------------# run simulations #--------------------------------#
# scenario A
res_MOSAIC = bind_rows(lapply(1:Nsim, run_mosaic, misspecified = FALSE))
res_BB = bind_rows(lapply(1:Nsim, run_BB, misspecified = FALSE))

# scenario B
res_MOSAIC_misp = bind_rows(lapply(1:Nsim, run_mosaic, misspecified = TRUE))
res_BB_misp = bind_rows(lapply(1:Nsim, run_BB, misspecified = TRUE))

save(res_MOSAIC, res_MOSAIC_misp, res_BB, res_BB_misp)

#-----------------# combined plot (Scenarios A & B) #--------------------------#

res_all = res_MOSAIC %>%
  left_join(res_BB, by = c("sim", "s")) %>%
  mutate(scenario = "A")
res_all_misp = res_MOSAIC_misp %>%
  left_join(res_BB_misp, by = c("sim", "s")) %>%
  mutate(scenario = "B")

res_plot = bind_rows(res_all, res_all_misp)

MSEplot = res_plot %>%
  pivot_longer(cols = c(MOSAIC, BB), names_to = "model", values_to = "mse") %>%
  ggplot(aes(x = factor(s), y = mse, fill = model)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  # ylim(c(0.03,0.5)) +
  facet_wrap(~ scenario) +
  scale_fill_manual(values = c("MOSAIC" = "firebrick2", "BB" = "steelblue")) +
  labs(x = "s", y = "MSE", fill = NULL) +
  theme_light()
MSEplot

# ggsave(plot = MSEplot, filename = "figures/MSE_sim.png", width = 8, height = 3.5)
