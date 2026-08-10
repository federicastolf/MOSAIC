
library(tidyverse)
library(rstan)
library(Rcpp)
library(parallel)
library(pheatmap)
library(gridExtra)
library(latex2exp)
library(ggpubr)   

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

rm(list = ls())

source("functs/gibbs_mosaic.R")
source("functs/helper.R")
source("functs/rarefaction.R")


################################################################################
#----------------------------#  data generation #------------------------------#

Q = 4                       # number of sites
N = 3000                    # total number of species 
beta_true = c(-3, 0.5, -0.3)  # regression coefficients (intercept + 2 covariates)
H = length(beta_true)
a_phi_true = 4              # gamma shape for phi_q
b_phi_true = 2              # gamma rate  for phi_q
n = c(50, 40, 35, 30)       # number of samples at each site
mseed = 134

data_sim = synt_dataH(n = n, Q = Q, beta = beta_true, a_phi = a_phi_true,
                      b_phi = b_phi_true, mseed = mseed, N = N)

Y = data_sim$Y[, colSums(data_sim$Y) > 0, drop = FALSE]
k = ncol(Y)                 # number of observed species 
m = data_sim$m              # k x Q matrix of species frequencies
site_names = paste0("S", 1:Q)

Xcov = data_sim$X # design matrix
colnames(Xcov) = c("Intercept", "x1", "x2")


################################################################################
#----------------------------#  fit the model #--------------------------------#

b0 = c(-2, rep(0, H - 1))
B0 = diag(1, H)
r0 = 300
nu0 = 1.3*k
a_phi = 2
b_phi = 2

# stan fit
stan_data = list(Q = Q, H = H, k = k, n_q = n, m = m, X = Xcov, b0 = b0, B0 = B0,
                 a_phi = a_phi, b_phi = b_phi, r0 = r0, nu0 = nu0)
fit_sim = stan(file = "functs/mosaic_NegBin.stan", data = stan_data, seed = mseed,
               chains = 4, iter = 2000, warmup = 1000)

# posterior of the regression coefficients 
zeta_post = rstan::extract(fit_sim, pars = c("beta"))[[1]]
colnames(zeta_post) = colnames(Xcov)
summ = data.frame(mean = colMeans(zeta_post),
                  lo   = apply(zeta_post, 2, quantile, probs = 0.025),
                  hi   = apply(zeta_post, 2, quantile, probs = 0.975))
summ   

#----# posterior distribution of N #----#

post_sim = rstan::extract(fit_sim, pars = c("alpha_out", "theta_out"))
Npost = genNNegBinom(post_sim, n, r0, nu0, k, mseed)
print(quantile(Npost, probs = c(0.025, 0.5, 0.975)))


#---# Data-augmentation steps  #------#
alpha_post = post_sim$alpha_out
theta_post = post_sim$theta_out
burnin = 1000
thin = 2
sigma6 = rep(0.6, Q)
occ6 = gibbs_mosaic_parallel(Y = Y, alpha_post = alpha_post,
                                 theta_post = theta_post, Npost = Npost, n = n,
                                 burnin = burnin, thin = thin, sigma = sigma6,
                                 mseed = mseed, n_cores = 4)


################################################################################
#---------------------------# compute beta-diversity #-------------------------#

betaDiv_post6 = BetaDiversity_Mat(occ6$eta)

# posterior median of the beta-similarity, arranged as a Q x Q matrix
beta_sim_matrix = function(betaDiv_post, labels = site_names) {
  betaSIM_median = apply(1 - betaDiv_post, 1, median)
  mat = diag(1, length(labels))
  mat[lower.tri(mat)] = betaSIM_median
  mat[upper.tri(mat)] = t(mat)[upper.tri(mat)]
  colnames(mat) = rownames(mat) = labels
  mat
}

matbSIM_median = beta_sim_matrix(betaDiv_post6)
hc = hclust(as.dist(1 - matbSIM_median))
ord_names = rownames(matbSIM_median)[hc$order]
C_ordered_plot = matbSIM_median[ord_names, ord_names]
diag(C_ordered_plot) = NA
brks_cmp = seq(0, 0.3, length.out = 101)
BD_plotmean = pheatmap(C_ordered_plot, treeheight_row = 0, treeheight_col = 0,
                       color = colorRampPalette(c("white", "steelblue1", "dodgerblue4"))(100),
                       cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10,
                       show_rownames = TRUE, labels_col = ord_names,
                       main = expression("Proposed " * beta * " similarity"),
                       breaks = brks_cmp, na_col = "blue4",
                       annotation_names_row = FALSE, annotation_legend = FALSE,
                       legend = FALSE)

################################################################################
#----------------------# rarefaction surfaces and curve #----------------------#

Ylist = convertY_tolist(Y, n)
Cpq = shared_species_pairs(Ylist, labels = site_names)   # observed shared species
Kq = compute_Kq(Y_list = Ylist, labels = site_names)     # observed richness
Cpq
Kq

#-------------------# out-of-sample shared-species surfaces #------------------#

s_extra = 20    # additional samples per site (the surface is s_extra x s_extra)
pairs_oos = list(list(p = 1, q = 4, name = "S1-S4"),
                 list(p = 2, q = 3, name = "S2-S3"),
                 list(p = 3, q = 4, name = "S3-S4"))

pred_surfaces = lapply(pairs_oos, function(pr) {
  pair_key = paste(site_names[pr$p], site_names[pr$q], sep = "_")
  arr = rar_pred_shared(k = k, k_p = Kq[pr$p], k_q = Kq[pr$q],
                        c_pq = Cpq[pair_key], r0 = r0, nu0 = nu0,
                        alpha = post_sim$alpha_out, theta = post_sim$theta_out,
                        n = n, s_pq = c(s_extra, s_extra),
                        p = pr$p, q = pr$q, full = TRUE)
  df = expand.grid(s = 1:s_extra, t = 1:s_extra)
  df$pair  = pr$name
  df$observed = as.numeric(Cpq[pair_key])   
  df$mean  = c(apply(arr, c(1, 2), mean))
  df$lower = c(apply(arr, c(1, 2), quantile, probs = 0.025))
  df$upper = c(apply(arr, c(1, 2), quantile, probs = 0.975))
  df
})

df_oos = do.call(rbind, pred_surfaces)
df_oos$pair = factor(df_oos$pair, levels = sapply(pairs_oos, `[[`, "name"))

# posterior mean surface and ci
df_plot = bind_rows(
  df_oos %>% mutate(value = lower, type = "Lower"),
  df_oos %>% mutate(value = mean,  type = "Mean"),
  df_oos %>% mutate(value = upper, type = "Upper")
)
df_plot$type = factor(df_plot$type, levels = c("Lower", "Mean", "Upper"))

breaks_custom = c(
  seq(78,  107, length.out = 10),   # S1-S4 and S3-S4
  seq(330, 436, length.out = 15)    # S2-S3
)
breaks_custom = sort(unique(breaks_custom))

plot_oos_surface = ggplot(df_plot, aes(x = s, y = t, z = value)) +
  geom_contour_filled(breaks = breaks_custom) +
  scale_fill_manual(
    values = colorRampPalette(c("lightblue1", "steelblue2", "royalblue3",
                                "darkblue", "black"))(length(breaks_custom) - 1)) +
  facet_grid(type ~ pair) +
  labs(fill = TeX("$C_{s_q,s_r}^{(n_q,n_r)}$"),
       x = TeX(" $s_q$"), y = TeX("$s_r$")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13),
        axis.text.x  = element_text(size = 10), axis.text.y  = element_text(size = 10),
        legend.title = element_text(size = 13), legend.text  = element_text(size = 10),
        strip.background.x = element_rect(fill = "grey70", color = NA),
        strip.background.y = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12), panel.grid = element_blank())
plot_oos_surface

#---------#  shared-species prediction curves with one-side test #------------#

pi_post = Map(function(g_list, e_list) {
  Map(function(g_vec, e_vec) g_vec * e_vec, g_list, e_list)
}, occ6$gamma, occ6$eta)

# helper: for a pair (q, p) with n_q <= n_p, the curve is
#   - empirical up to n_q      (both sites observed)
#   - "semi-test" up to n_p    (site q augmented, site p observed)
#   - fully predicted afterwards
build_shared_data = function(q, p, pair_name, s_aug, B) {
  n_obs_q = n[q] ; n_obs_p = n[p]
  u = n_obs_p - n_obs_q
  stopifnot(u >= 0)          # site q must be the one with the smaller effort
  s_aug_q = u + s_aug
  s_ext = s_aug
  stopifnot(B*s_aug_q <= length(pi_post[[q]]))
  stopifnot(B*s_ext <= length(pi_post[[p]]))
  
  all_idx_q = 1:(B * s_aug_q)
  all_idx_p = 1:(B * s_ext)
  max_lenq = max(sapply(pi_post[[q]][all_idx_q], length))
  max_lenp = max(sapply(pi_post[[p]][all_idx_p], length))
  max_len = max(max_lenq, max_lenp, ncol(Ylist[[q]]), ncol(Ylist[[p]]))
  Yq_base = cbind(Ylist[[q]], matrix(0, nrow(Ylist[[q]]), max_len - ncol(Ylist[[q]])))
  Yp_base = cbind(Ylist[[p]], matrix(0, nrow(Ylist[[p]]), max_len - ncol(Ylist[[p]])))
  
  curves = lapply(1:B, function(b) {
    idx_q = ((b - 1)*s_aug_q + 1):(b*s_aug_q)
    Yaugq = matrix(0, nrow = s_aug_q, ncol = max_len)
    for (i in seq_along(idx_q)) {
      pii = pi_post[[q]][[idx_q[i]]]
      Yaugq[i, 1:length(pii)] = rbinom(length(pii), 1, pii)
    }
    Yaugq_all = rbind(Yq_base, Yaugq)
    idx_p = ((b - 1)*s_ext + 1):(b*s_ext)
    Yaugp = matrix(0, nrow = s_ext, ncol = max_len)
    for (i in seq_along(idx_p)) {
      piip = pi_post[[p]][[idx_p[i]]]
      Yaugp[i, 1:length(piip)] = rbinom(length(piip), 1, piip)
    }
    Yaugp_all = rbind(Yp_base, Yaugp)
    simplified_accum(n_pq = c(nrow(Yaugp_all), nrow(Yaugq_all)),
                     Yp = Yaugp_all, Yq = Yaugq_all, nperm = 50)
  })
  
  mat_curves = do.call(cbind, curves)
  mean_curve = rowMeans(mat_curves)
  lo_curve = apply(mat_curves, 1, quantile, 0.025)
  hi_curve = apply(mat_curves, 1, quantile, 0.975)
  n_total = nrow(mat_curves)
  
  make_df = function(idx, seg)
    data.frame(pair = pair_name, x = idx, mean = mean_curve[idx], lo = lo_curve[idx],
               hi = hi_curve[idx], segment = seg)
  
  list(emp = make_df(1:n_obs_q, "emp"), semi = make_df(n_obs_q:n_obs_p, "semi"),
       pred = make_df(n_obs_p:n_total, "pred"),
       vlines = data.frame(pair = pair_name, xint = c(n_obs_q, n_obs_p)))
}

# each pair is given as (q, p) with n_q <= n_p
pairs_list = list(list(q = 4, p = 1, name = "S4-S1"),
                  list(q = 3, p = 2, name = "S3-S2"),
                  list(q = 4, p = 3, name = "S4-S3"))
s_aug = 20
B = 20

results = lapply(pairs_list, function(pr) {
  cat("Processing pair:", pr$name, "\n")
  build_shared_data(pr$q, pr$p, pr$name, s_aug = s_aug, B = B)
})

all_emp = do.call(rbind, lapply(results, `[[`, "emp"))
all_semi = do.call(rbind, lapply(results, `[[`, "semi"))
all_pred = do.call(rbind, lapply(results, `[[`, "pred"))
all_vline = do.call(rbind, lapply(results, `[[`, "vlines"))

pair_levels = sapply(pairs_list, `[[`, "name")
all_emp$pair = factor(all_emp$pair, pair_levels)
all_semi$pair = factor(all_semi$pair, pair_levels)
all_pred$pair = factor(all_pred$pair, pair_levels)
all_vline$pair = factor(all_vline$pair, pair_levels)

pall_semi = ggplot() +
  geom_vline(data = all_vline, aes(xintercept = xint), linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_point(data = all_emp, aes(x, mean), size = 0.4, color = "black") +
  geom_ribbon(data = all_semi, aes(x, ymin = lo, ymax = hi), fill = "orange",
              alpha = 0.4) +
  geom_line(data = all_semi, aes(x, mean), linewidth = 0.7, color = "orange") +
  geom_ribbon(data = all_pred, aes(x, ymin = lo, ymax = hi), fill = "steelblue",
              alpha = 0.4) +
  geom_line(data = all_pred, aes(x, mean), linewidth = 0.7, color = "steelblue") +
  facet_wrap(~ pair, nrow = 1, scales = "free_x") +
  labs(x = "Sampling effort", y = TeX("$C_{n_q, n_r}$")) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey70", color = NA),
    strip.text = element_text(size = 19), axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 17), axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )
pall_semi

#----------# site-specific rarefaction curves for different sigma #---------#

sites_info = list(list(site = "S1", q = 1), list(site = "S2", q = 2),
                  list(site = "S3", q = 3))

sigma_vals = c(0.1, 0.5, 0.9)
sigma_colors = setNames(
  colorRampPalette(c("steelblue3", "darkgoldenrod", "firebrick"))(length(sigma_vals)),
  as.character(sigma_vals)
)

# helper:
make_site_plot = function(site_name, q, post_param, n, sigma_vals,
                          sigma_colors, r0, nu0, mseed, ylim = NULL,
                          show_legend = FALSE, show_y_axis = TRUE) {
  B = nrow(post_param$alpha_out)
  pred_list = lapply(sigma_vals, function(sig) {
    pred = rar_pred_Kgroup(alphaq = post_param$alpha_out[, q],
                           thetaq = post_param$theta_out[, q], r0 = r0, nu0 = nu0,
                           nq = 0, kq = 0, sq = n[q], seed = mseed,
                           sigmaq = sig, p0_all = rep(1, B), k = 0)
    df = data.frame(s = 1:nrow(pred), mean = rowMeans(pred),
                    lo = apply(pred, 1, quantile, probs = 0.025),
                    hi = apply(pred, 1, quantile, probs = 0.975),
                    sigma = as.factor(sig))
    rbind(data.frame(s = 0, mean = 0, lo = 0, hi = 0, sigma = as.factor(sig)), df)
  })
  df_pred = do.call(rbind, pred_list)
  df_pred$facet = site_name   # dummy variable for the facet strip
  
  p = ggplot() +
    geom_ribbon(data = df_pred, aes(x = s, ymin = lo, ymax = hi, fill = sigma),
                alpha = 0.2) +
    geom_line(data = df_pred, aes(x = s, y = mean, color = sigma), linewidth = 0.9) +
    scale_color_manual(values = sigma_colors, name = expression(sigma)) +
    scale_fill_manual(values  = sigma_colors, name = expression(sigma)) +
    coord_cartesian(ylim = ylim) +
    labs(x = "Sampling effort", y = TeX("$L_{q,n_q}$")) +
    facet_wrap(~ facet) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 10),
      axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 12), legend.text = element_text(size = 10),
      strip.background = element_rect(fill = "grey70", color = NA),
      strip.text = element_text(size = 12)
    )
  if (!show_legend) p = p + theme(legend.position = "none")
  if (!show_y_axis)
    p = p + theme(axis.title.y = element_blank(),
                  axis.text.y  = element_blank(),
                  axis.ticks.y = element_blank())
  p
}

B = nrow(alpha_post)
global_ylim = local({
  y_pred_hi = unlist(lapply(sites_info, function(s) {
    unlist(lapply(sigma_vals, function(sig) {
      pred = rar_pred_Kgroup(
        alphaq = post_sim$alpha_out[, s$q], thetaq = post_sim$theta_out[, s$q],
        r0 = r0, nu0 = nu0, nq = 0, kq = 0, sq = n[s$q], seed = mseed,
        sigmaq = sig, p0_all = rep(1, B), k = 0
      )
      apply(pred, 1, quantile, probs = 0.975)
    }))
  }))
  c(0, max(y_pred_hi, na.rm = TRUE) * 1.05)
})

plot_list = lapply(seq_along(sites_info), function(i) {
  s = sites_info[[i]]
  make_site_plot(site_name = s$site, q = s$q,
                 post_param = post_sim, n = n, sigma_vals = sigma_vals,
                 sigma_colors = sigma_colors, r0 = r0, nu0 = nu0, mseed = mseed,
                 ylim = global_ylim,
                 show_legend = (i == length(sites_info)),
                 show_y_axis = (i == 1))
})

shared_legend = get_legend(plot_list[[length(plot_list)]])
plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] +
  theme(legend.position = "none")

psigma_group = grid.arrange(
  arrangeGrob(grobs = lapply(plot_list, ggplotGrob), nrow = 1),
  shared_legend,
  ncol = 2,
  widths = c(0.9, 0.1)
)



