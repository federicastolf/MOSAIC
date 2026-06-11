library(tidyverse)
library(rstan)
library(Rcpp)
library(latex2exp)
library(gridExtra)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

rm(list=ls())

source("functs/helper.R")
source("functs/rarefaction.R")

#---# generate data #---#

beta = c(-3, 0.5, -0.3)
H = length(beta)
Q = 3
a_phi_true = 4
b_phi_true = 2
mseed = 134
n = c(100,80,80)
r0_true = 1000
nu0_true = 3000
N = 3000

data_sim = synt_dataH(n = n, Q = Q, beta = beta, a_phi = a_phi_true,
                      b_phi = b_phi_true, mseed = mseed, N = N)
tt_split = split_train_test(data_sim, train_prop = 0.6, mseed = mseed) 

#---# fit model #---#

b0 = c(-2,rep(0,H-1))
B0 = diag(1, H)
r0 = 10 
nu0 = 1.5*tt_split$k_train 
a_phi = 2
b_phi = 2

stan_data = list(Q = Q, H = H, k = tt_split$k_train, n_q = tt_split$n_train, 
                 m = tt_split$m_train, X = data_sim$X, b0 = b0, B0 = B0, 
                 a_phi = a_phi, b_phi = b_phi, r0 = r0, nu0 = nu0)
fit_sim = stan(file = "functs/mosaic_NegBin.stan", data = stan_data, seed = mseed,
                 chains = 4, iter = 2000, warmup = 1000)

post_param = rstan::extract(fit_sim, pars = c("alpha_out", "theta_out"))


#--------# group shared accumulation curve simplified for group 2-3 #----------#

# empirical
p_23 = 2
q_23 = 3

emp23_train = simplified_accum(n_pq = tt_split$n_train[c(p_23,q_23)], 
                               Yp = tt_split$Y_train[[p_23]], 
                               Yq = tt_split$Y_train[[q_23]], nperm = 300)
emp23_test = simplified_accum_anchor(ntest_pq = tt_split$n_test[c(p_23,q_23)], 
                                     Yp_test = tt_split$Y_test[[p_23]], 
                                     Yq_test = tt_split$Y_test[[q_23]], 
                                     Yp_train = tt_split$Y_train[[p_23]], 
                                     Yq_train = tt_split$Y_train[[q_23]], nperm = 300)
emp23_data = data.frame(s=c(1:data_sim$n[p_23]), y = c(emp23_train, emp23_test))

# predicted
c_train = shared_species_pairs(tt_split$Y_train)
c_train
pred23 = rar_pred_shared(k = tt_split$k_train, k_p = tt_split$kq_train[p_23], 
                         k_q = tt_split$kq_train[q_23], c_pq = c_train["2_3"], 
                         r0 = r0, nu0 = nu0, alpha = post_param$alpha_out, 
                         theta = post_param$theta_out, n = tt_split$n_train, 
                         s_pq = tt_split$n_test[c(p_23,q_23)], p = p_23, q = q_23, 
                         full = FALSE)

sum_pred23 = data.frame(
  s = c(tt_split$n_train[p_23]:((tt_split$n_train[p_23] + tt_split$n_test[p_23]))), 
  mean = c(tail(emp23_train,n=1), rowMeans(pred23)), 
  lo = c(tail(emp23_train,n=1), apply(pred23, 1, quantile, probs = 0.025)), 
  hi = c(tail(emp23_train,n=1), apply(pred23, 1, quantile, probs = 0.975))
)


# plot
pcurve_23 = ggplot() +
  geom_point(data = emp23_data, aes(x = s, y = y), size = 0.3, color = "grey20") +
  geom_ribbon(data = sum_pred23, aes(x = s, ymin = lo, ymax = hi), alpha = 0.4, 
              fill = "steelblue") +
  geom_line(data = sum_pred23, aes(x = s, y = mean), linewidth = 1, 
            color = "steelblue") +
  geom_vline(xintercept = tt_split$n_train[p_23], linetype = 2, color = "gray50") +
  labs(x = "Sampling effort", y = TeX("$C_{n_2,n_3}$")) +
  facet_wrap(~ "2-3") +                        # variabile fittizia con label "2-3"
  theme_minimal() +
  theme(strip.background = element_rect(fill = "grey70", color = NA),
        strip.text = element_text(size = 18), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18), axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12))

pcurve_23

# ggsave(plot = pcurve_23, filename = "figures/sim_simplified_main.png", 
#        width = 9, height = 5.5)


#---------# group shared accumulation surface - test prediction #--------------#

pairs_list = list(c(1,2), c(1,3), c(2,3))
pair_names = c("1_2", "1_3", "2_3")

surface_results = mapply(function(pq, name) {
  compute_surface_pair(
    p = pq[1], q = pq[2],
    tt_split = tt_split,
    c_train = c_train,
    post_param = post_param,
    r0 = r0, nu0 = nu0,
    pair_key = name,
    nperm = 300
  )
}, pairs_list, pair_names, SIMPLIFY = FALSE)
names(surface_results) = pair_names

# unpack results
dfemp_12_test = surface_results[["1_2"]]$dfemp
dfemp_13_test = surface_results[["1_3"]]$dfemp
dfemp_23_test = surface_results[["2_3"]]$dfemp
dfpred_12 = surface_results[["1_2"]]$dfpred
dfpred_13 = surface_results[["1_3"]]$dfpred
dfpred_23 = surface_results[["2_3"]]$dfpred


#---# plot posterior mean vs empirical #---#

shared_limits <- range(c(dfpred_13$mean, dfemp_13_test$shared, dfpred_12$mean, 
                         dfemp_12_test$shared, dfpred_23$mean, dfemp_23_test$shared))

df_surface_main <- bind_rows(
  dfemp_12_test %>% mutate(value = shared, pair = "1-2", type = "Empirical"),
  dfemp_13_test %>% mutate(value = shared, pair = "1-3", type = "Empirical"),
  dfemp_23_test %>% mutate(value = shared, pair = "2-3", type = "Empirical"),
  dfpred_12     %>% mutate(value = mean,   pair = "1-2", type = "Predicted"),
  dfpred_13     %>% mutate(value = mean,   pair = "1-3", type = "Predicted"),
  dfpred_23     %>% mutate(value = mean,   pair = "2-3", type = "Predicted")
)

df_surface_main$type <- factor(df_surface_main$type, 
                               levels = c("Empirical", "Predicted"))
df_surface_main$pair <- factor(df_surface_main$pair, levels = c("1-2", "1-3", "2-3"))
 
# plot_share_main = ggplot(df_surface_main, aes(x = s, y = t, fill = value)) +
#   geom_raster(interpolate = FALSE) +
#   scale_fill_gradientn(colors = c("white", "lightblue", "steelblue3", "darkblue",
#                                   "black"), limits = shared_limits) +
#   facet_grid(type ~ pair) +                  
#   labs(fill = TeX("$C_{s_p,s_q}^{(n_p,n_q)}$"), x = TeX("$s_p$"), y = TeX("$s_q$")) +
#   theme_minimal() +
#   theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
#         axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
#         legend.title = element_text(size = 13), legend.text = element_text(size = 10),
#         plot.title = element_text(hjust = 0.5, size = 13),
#         strip.background.x = element_rect(fill = "grey70", color = NA),  
#         strip.background.y = element_blank(), strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12))
# plot_share_main

plot_share_main_o2 = ggplot(df_surface_main, aes(x = s, y = t, z = value)) +
  geom_contour_filled(bins = 15) +
  scale_fill_manual(values = colorRampPalette(c("white", "lightblue", "steelblue3", 
                                                "darkblue","black"))(15)) +
  facet_grid(type ~ pair) +
  labs(fill = TeX("$C_{s_p,s_q}^{(n_p,n_q)}$"), x = TeX("$s_p$"), y = TeX("$s_q$")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.x  = element_text(size = 10), axis.text.y  = element_text(size = 10),
        legend.title = element_text(size = 13), legend.text  = element_text(size = 10),
        strip.background.x = element_rect(fill = "grey70", color = NA),
        strip.background.y = element_blank(),
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), panel.grid = element_blank())
plot_share_main_o2

# ggsave(plot = plot_share_main_o2, filename = "figures/sim_share_main_o2.png",
#        width = 10, height = 5)


#---# plot CI model #---#

df_surface_supp <- bind_rows(
  dfpred_12 %>% mutate(value = lower, pair = "1-2", type = "Lower"),
  dfpred_13 %>% mutate(value = lower, pair = "1-3", type = "Lower"),
  dfpred_23 %>% mutate(value = lower, pair = "2-3", type = "Lower"),
  dfpred_12 %>% mutate(value = mean,  pair = "1-2", type = "Mean"),
  dfpred_13 %>% mutate(value = mean,  pair = "1-3", type = "Mean"),
  dfpred_23 %>% mutate(value = mean,  pair = "2-3", type = "Mean"),
  dfpred_12 %>% mutate(value = upper, pair = "1-2", type = "Upper"),
  dfpred_13 %>% mutate(value = upper, pair = "1-3", type = "Upper"),
  dfpred_23 %>% mutate(value = upper, pair = "2-3", type = "Upper")
)

df_surface_supp$type = factor(df_surface_supp$type, levels = c("Lower", "Mean",
                                                               "Upper"))
df_surface_supp$pair = factor(df_surface_supp$pair, levels = c("1-2", "1-3", "2-3"))

shared_limits_supp = range(c(dfpred_12$lower, dfpred_12$upper,dfpred_13$lower, 
                             dfpred_13$upper, dfpred_23$lower, dfpred_23$upper))

# supp_share = ggplot(df_surface_supp, aes(x = s, y = t, fill = value)) +
#   geom_raster(interpolate = FALSE) +
#   scale_fill_gradientn(colors = c("white", "lightblue", "steelblue3", "darkblue", 
#                                   "black"), limits = shared_limits_supp) +
#   facet_grid(type ~ pair) +
#   labs(fill = TeX("$C_{s_p,s_q}^{(n_p,n_q)}$"), x = TeX("$s_p$"), y = TeX("$s_q$")) +
#   theme_minimal() +
#   theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
#         axis.text.x = element_text(size = 10), axis.text.y  = element_text(size = 10),
#         legend.title = element_text(size = 13), legend.text  = element_text(size = 10),
#         plot.title = element_text(hjust = 0.5, size = 13),
#         strip.background.x = element_rect(fill = "grey70", color = NA),
#         strip.background.y = element_blank(),          
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_blank())       
# supp_share

supp_share = ggplot(df_surface_supp, aes(x = s, y = t, z = value)) +
  geom_contour_filled(bins = 15) +
  scale_fill_manual(values = colorRampPalette(c("white", "lightblue", "steelblue3", 
                                                "darkblue","black"))(15)) +
  facet_grid(type ~ pair) +
  labs(fill = TeX("$C_{s_p,s_q}^{(n_p,n_q)}$"), x = TeX("$s_p$"), y = TeX("$s_q$")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.x  = element_text(size = 10), axis.text.y  = element_text(size = 10),
        legend.title = element_text(size = 13), legend.text  = element_text(size = 10),
        strip.background.x = element_rect(fill = "grey70", color = NA),
        strip.background.y = element_blank(),
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), panel.grid = element_blank())
supp_share

# ggsave(plot = supp_share, filename = "figures/sim_share_suppUQ.png",
#        width = 10, height = 7)

