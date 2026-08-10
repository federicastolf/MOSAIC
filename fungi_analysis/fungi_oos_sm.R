library(tidyverse)
library(rstan)
library(parallel)
library(latex2exp)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

rm(list=ls())

load("data/fungiGSSP.Rdata")
source("functs/helper.R")
source("functs/rarefaction.R")

# ----------------------------------# fit model #------------------------------#

Q = nrow(Xgroup)
Xcov = Xgroup
Xcov$temp_mean2 = Xcov$temp_mean^2
Xcov = Xcov %>% dplyr::select(-c(cz, site, n, lon))
Xcov = scale(Xcov)
Xcov = cbind(rep(1, nrow(Xcov)), Xcov)
head(Xcov)
H = ncol(Xcov)
  
# split train-test
data_all = list("Y" = Y, "n" = Xgroup$n)
mseed = 434
tt_split = split_train_test(data_all, train_prop = 0.7, mseed = mseed) 

# hyperparameters
b0 = rep(0, H)
B0 = diag(1, H)
r0 = 300 
nu0 = 1.3*tt_split$k_train
a_phi = 2
b_phi = 0.5 

# fit stan model 
stan_data = list(Q = Q, H = H, k = tt_split$k_train, n_q = tt_split$n_train,
                 m = tt_split$m_train, X = Xcov, b0 = b0, B0 = B0, a_phi = a_phi, 
                 b_phi = b_phi, r0 = r0, nu0 = nu0)
fit_fungi_train = stan(file = "functs/mosaic_NegBin.stan", data = stan_data, seed = 434,
                 chains = 4, iter = 2000, warmup = 1000)

#load("/Users/stolffederica/Library/CloudStorage/Dropbox/ActiveResearch/Pex_IBP/code/code_pBB/fungi_analysis/resDCC/stan_fungi300_train2.Rdata")

post_param = rstan::extract(fit_fungi_train, pars = c("alpha_out", "theta_out"))
c_train = shared_species_pairs(tt_split$Y_train)

#--------------------------# prediction on shared species #--------------------#

pairs_list = list(c(11,21), c(4,10), c(28,30))
pair_names2 = c("HAI-SAN", "BAV-FTS", "VII-YAK")
pair_names = c("11_21", "4_10", "28_30") # put before lower number between the two

surface_results = mapply(function(pq, name) {
  compute_surface_pair(p = pq[1], q = pq[2], tt_split = tt_split, c_train = c_train,
                       post_param = post_param, r0 = r0, nu0 = nu0, pair_key = name,
                       nperm = 100)
}, pairs_list, pair_names, SIMPLIFY = FALSE)
names(surface_results) = pair_names2

# unpack results
dfemp_12_test = surface_results[["HAI-SAN"]]$dfemp
dfemp_13_test = surface_results[["BAV-FTS"]]$dfemp
dfemp_23_test = surface_results[["VII-YAK"]]$dfemp
dfpred_12 = surface_results[["HAI-SAN"]]$dfpred
dfpred_13 = surface_results[["BAV-FTS"]]$dfpred
dfpred_23 = surface_results[["VII-YAK"]]$dfpred

# plot posterior mean vs empirical 

# shared_limits <- range(c(dfpred_13$mean, dfemp_13_test$shared, dfpred_12$mean, 
#                          dfemp_12_test$shared, dfpred_23$mean, dfemp_23_test$shared))

df_surface_main = bind_rows(
  dfemp_12_test %>% mutate(value = shared, pair = "HAI-SAN", type = "Empirical"),
  dfemp_13_test %>% mutate(value = shared, pair = "BAV-FTS", type = "Empirical"),
  dfemp_23_test %>% mutate(value = shared, pair = "VII-YAK", type = "Empirical"),
  dfpred_12     %>% mutate(value = mean,   pair = "HAI-SAN", type = "Predicted"),
  dfpred_13     %>% mutate(value = mean,   pair = "BAV-FTS", type = "Predicted"),
  dfpred_23     %>% mutate(value = mean,   pair = "VII-YAK", type = "Predicted")
)

df_surface_main$type = factor(df_surface_main$type, 
                              levels = c("Empirical", "Predicted"))
df_surface_main$pair = factor(df_surface_main$pair, 
                               levels = c("HAI-SAN", "BAV-FTS", "VII-YAK"))

breaks_custom = c(
    seq(306, 390,  length.out = 8),   # first range 
    seq(128, 145,  length.out = 5),   # third range
    seq(630, 722, length.out = 16)   # second range
  )
breaks_custom = sort(unique(breaks_custom))

plot_fungi_oos = ggplot(df_surface_main, aes(x = s, y = t, z = value)) +
  geom_contour_filled(breaks = breaks_custom) +
  scale_fill_manual(
    values = colorRampPalette(c("lightblue1", "steelblue2","royalblue3","darkblue", 
                                "black"))(length(breaks_custom) - 1)) +
  facet_grid(type ~ pair) +
  labs(fill = TeX("$C_{s_q,s_r}^{(n_q,n_r)}$"), x = TeX("$s_q$"), y = TeX("$s_r$")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.x  = element_text(size = 10), axis.text.y  = element_text(size = 10),
        legend.title = element_text(size = 13), legend.text  = element_text(size = 10),
        strip.background.x = element_rect(fill = "grey70", color = NA),
        strip.background.y = element_blank(),
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12), panel.grid = element_blank())
plot_fungi_oos 

# ggsave(plot = plot_fungi_oos, filename = "fungi_analysis/figures/plot_fungi_oos.png",
#        width = 12, height = 5.15)

