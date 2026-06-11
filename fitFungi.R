library(tidyverse)
library(rstan)
library(parallel)
library(pheatmap)
library(gridExtra)
library(latex2exp)
library(ggpubr)  # for get_legend

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

rm(list=ls())

load("data/fungiGSSP.Rdata")
source("functs/gibbs_mosaic.R")
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
n = Xgroup$n
k = ncol(Y)

m = matrix(0, nrow = ncol(Y), ncol = Q) 
ends_n = cumsum(n)
starts_n = c(1, head(ends_n, -1) + 1)
for (q in 1:Q) {
  Yq = Y[starts_n[q]:ends_n[q], ]
  m[, q] = colSums(Yq)
}

# hyperparameters
b0 = rep(0, H)
B0 = diag(1, H)
r0 = 300 
nu0 = 1.3*k 
a_phi = 2
b_phi = 2

# fit stan model 
stan_data = list(Q = Q, H = H, k = k, n_q = n, m = m, X = Xcov, b0 = b0, B0 = B0, 
                 a_phi = a_phi, b_phi = b_phi, r0 = r0, nu0 = nu0)
fit_fungi = stan(file = "functs/mosaic_NegBin.stan", data = stan_data, seed = 434,
                 chains = 4, iter = 2000, warmup = 1000)

# load("/Users/stolffederica/Library/CloudStorage/Dropbox/ActiveResearch/Pex_IBP/code/code_pBB/fungi_analysis/resDCC/stan_fungi4_gssp300.Rdata")

#----# posterior of zeta #----#
zeta_post = rstan::extract(fit_fungi, pars = c("beta"))[[1]]
colnames(zeta_post) = colnames(Xcov)
colnames(zeta_post)[1] = "Intercept"
summ = data.frame(zeta = colnames(zeta_post), mean = colMeans(zeta_post),
                   lo   = apply(zeta_post, 2, quantile, probs = 0.025),
                   hi   = apply(zeta_post, 2, quantile, probs = 0.975))
summ$zeta = factor(summ$zeta)
summ

#------# generate posterior N #---------#
post_fungi = rstan::extract(fit_fungi, pars = c("alpha_out", "theta_out"))
Npost = genNNegBinom(post_fungi, n, r0, nu0, k, 434)
mean(Npost)
quantile(Npost, probs = c(0.025, 0.975))
alpha_post = post_fungi$alpha_out
theta_post = post_fungi$theta_out


################################################################################
#--# fit DA gibbs and compute beta-diversity for different values of sigma #--#

burnin = 1000
thin = 2
mseed = 343

# sigma = 0.3
sigma3 = rep(0.3, Q)
occ_fungi3 = gibbs_mosaic_parallel(Y = Y, alpha_post = alpha_post, 
                                 theta_post = theta_post, Npost = Npost, n = n,
                                 burnin = burnin, thin = thin, sigma = sigma3,
                                 mseed = mseed, n_cores = 34)
betaDiv_post3 = BetaDiversity_Mat(occ_fungi3$eta)
rm(occ_fungi3)

# sigma = 0.6
sigma6 = rep(0.6, Q)
occ_fungi6 = gibbs_mosaic_parallel(Y = Y, alpha_post = alpha_post, 
                                 theta_post = theta_post, Npost = Npost, n = n, 
                                 burnin = burnin, thin = thin, sigma = sigma6,
                                 mseed = mseed, n_cores = 34)
betaDiv_post6 = BetaDiversity_Mat(occ_fungi6$eta)
rm(occ_fungi6)

# sigma = 0.9
sigma9 = rep(0.9, Q)
occ_fungi9 = gibbs_mosaic_parallel(Y = Y, alpha_post = alpha_post, 
                                 theta_post = theta_post, Npost = Npost, n = n, 
                                 burnin = burnin, thin = thin, sigma = sigma9, 
                                 mseed = mseed, n_cores = 34)
betaDiv_post9 = BetaDiversity_Mat(occ_fungi9$eta)
rm(occ_fungi9)

load("/Users/stolffederica/Library/CloudStorage/Dropbox/ActiveResearch/Pex_IBP/code/code_pBB/fungi_analysis/resDCC/betaDiv300_gssp.Rdata")

# helper plot comparison between different values of sigma
make_beta_heatmap = function(betaDiv_post, title, ord_names = NULL) {
  betaSIM_median = apply(1 - betaDiv_post, 1, median)
  matbSIM_median = diag(1, Q)
  matbSIM_median[lower.tri(matbSIM_median)] = betaSIM_median
  matbSIM_median[upper.tri(matbSIM_median)] = t(matbSIM_median)[upper.tri(matbSIM_median)]
  colnames(matbSIM_median) = rownames(matbSIM_median) = names(table(Xgroup$site))
  
  # order by cz clustering
  X_CZ = Xgroup %>%
    dplyr::select(site, cz) %>%
    distinct(site, cz)
  X_CZ = X_CZ[match(rownames(matbSIM_median), X_CZ$site), ]
  cz_levels = unique(X_CZ$cz)
  
  if (is.null(ord_names)) {
    ord_names = unlist(lapply(cz_levels, function(z) {
      sites_z = X_CZ$site[X_CZ$cz == z]
      submat = matbSIM_median[sites_z, sites_z, drop = FALSE]
      if (length(sites_z) > 1) {
        hc = hclust(as.dist(1 - submat))
        sites_z[hc$order]
      } else {
        sites_z
      }
    }))
  }
  C_ordered = matbSIM_median[ord_names, ord_names]
  annotation_row = data.frame(CZ = X_CZ$cz[match(ord_names, X_CZ$site)])
  rownames(annotation_row) = ord_names
  lb_site = ord_names
  C_ordered_plot = C_ordered
  diag(C_ordered_plot) = NA
  
  p = pheatmap(C_ordered_plot, treeheight_row = 0, treeheight_col = 0,
               color = colorRampPalette(c("white", "steelblue1", "dodgerblue4"))(100),
               cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 14, 
               show_rownames = FALSE, show_colnames = FALSE, labels_col = lb_site,
               main = title, breaks = seq(0.003, 0.9, length.out = 101), 
               na_col = "blue4", annotation_names_row = FALSE, 
               annotation_legend = FALSE, legend = FALSE)
  return(list(plot = p, ord_names = ord_names))
}

p9 = make_beta_heatmap(betaDiv_post9, expression(sigma == 0.9))
ref_ord = p9$ord_names
p3 = make_beta_heatmap(betaDiv_post3, expression(sigma == 0.3), ord_names = ref_ord)
p6 = make_beta_heatmap(betaDiv_post6, expression(sigma == 0.6), ord_names = ref_ord)

comp_beta = grid.arrange(p3$plot$gtable, p6$plot$gtable, p9$plot$gtable, nrow = 1)

# ggsave(plot = comp_beta, filename = "fungi_analysis/figures/betadiv_sigma.png",
#        width = 17, height = 5.5)


#-------# comparison with Jaccard #----------#

betaSIM_median = apply(1 - betaDiv_post6, 1, median)
matbSIM_median = diag(1, Q)
matbSIM_median[lower.tri(matbSIM_median)] = betaSIM_median
matbSIM_median[upper.tri(matbSIM_median)] = t(matbSIM_median)[upper.tri(matbSIM_median)]
colnames(matbSIM_median) = rownames(matbSIM_median) = names(table(Xgroup$site))

ord_names = ref_ord
X_CZ = Xgroup %>%
  dplyr::select(site, cz) %>%
  distinct(site, cz)
annotation_row = data.frame(CZ = X_CZ$cz[match(ord_names, X_CZ$site)])
rownames(annotation_row) = ord_names
lb_site = ord_names

C_ordered_plot = matbSIM_median[ord_names, ord_names]
diag(C_ordered_plot) = NA

BD_plotmean = pheatmap(C_ordered_plot, treeheight_row = 0, treeheight_col = 0,
                       color = colorRampPalette(c("white", "steelblue1", "dodgerblue4"))(100),
                       cluster_rows = FALSE, cluster_cols = FALSE,  fontsize = 10,
                       show_rownames = TRUE, labels_col = lb_site, 
                       main = expression("Proposed " * beta * " similarity"),
                       breaks = seq(0.003, 0.9, length.out = 101),
                       na_col = "blue4", annotation_names_row = FALSE, 
                       annotation_row = annotation_row, annotation_legend = FALSE, 
                       legend = FALSE)

# Jaccard
lab_group = Xgroup$site
Jac = 1 - Betadiss_data(Y, Xgroup$n, "Jaccard")
Jac_mat = diag(1, length(lab_group))
Jac_mat[lower.tri(Jac_mat)] = Jac
Jac_mat[upper.tri(Jac_mat)] = t(Jac_mat)[upper.tri(Jac_mat)]
colnames(Jac_mat) = rownames(Jac_mat) <- lab_group

Jac_mat1 = Jac_mat[ord_names, ord_names]
diag(Jac_mat1) = NA

pJac = pheatmap(Jac_mat1, treeheight_row = 0, treeheight_col = 0,
                color = colorRampPalette(c("white", "steelblue1", "dodgerblue4"))(100),
                cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10, 
                legend = FALSE, show_rownames = TRUE, labels_col = lb_site,
                breaks = seq(0.003, 0.9, length.out = 101), na_col = "blue4", 
                main = expression("Jaccard similarity"), 
                annotation_names_row = FALSE, annotation_row = annotation_row,
                annotation_legend = FALSE)

comp_betadiv = grid.arrange(BD_plotmean$gtable, pJac$gtable, nrow = 1)

# ggsave(plot = comp_betadiv, filename = "fungi_analysis/figures/comp_betadiv.png",
#        width = 13, height = 6.5)


#----# beta-diversity density plot #----#

site_names = Xgroup$site  
pair_index = data.frame(row = combn(1:Q, 2)[1,], col = combn(1:Q, 2)[2,],
                        site1 = combn(site_names, 2)[1,],
                        site2 = combn(site_names, 2)[2,])

selected_pairs = list(c("KLA", "FTS"), c("NUU", "PRI"), c("KLA", "ZUR"))

df_pairs = do.call(rbind, lapply(selected_pairs, function(pair) {
  site_a = pair[1]; site_b = pair[2]
  idx = which(
    (pair_index$site1 == site_a & pair_index$site2 == site_b) |
      (pair_index$site1 == site_b & pair_index$site2 == site_a)
  )
  if (length(idx) == 0) stop(paste("Pair not found:", site_a, site_b))
  data.frame(value = betaDiv_post6[idx, ],
             pair  = paste(site_a, site_b, sep = "-"))
}))

pair_levels = sapply(selected_pairs, function(p) paste(p[1], p[2], sep = "-"))
df_pairs$pair = factor(df_pairs$pair, levels = pair_levels)

# plot
pdens = ggplot(df_pairs, aes(x = value)) +
  geom_density(fill = NA, color = "black", adjust = 2, linewidth = 0.25) +
  facet_grid(. ~ pair, scales = "free") +
  labs(x = expression(beta[pq]), y = "Density") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3.5)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 6.5), 
        axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),
        strip.background.x = element_rect(fill = "grey70", color = NA),
        strip.text.x  = element_text(size = 7.5))
pdens

# ggsave(plot = pdens, filename = "fungi_analysis/figures/density_beta.png",
#        width = 6, height = 2)

################################################################################
#------# group-specific rarefaction curves for different values of sigma #-----#

Ylist = convertY_tolist(Y, Xgroup$n)
Cpq = shared_species_pairs(Ylist, labels = Xgroup$site)
Kq = compute_Kq(Y_list = Ylist, labels = Xgroup$site)

mseed = 343
# empirical
accum_groups = accum_by_group(Y = Y, n = Xgroup$n, lab_group = Xgroup$site,
                              method = "random", permutations = 500,
                              mseed = mseed)

sites_info = list(list(site = "KLA", q = 15), list(site = "TOO", q = 26), 
                  list(site = "FTS", q = 10))

sigma_vals = c(0.05, 0.5, 0.95)
sigma_colors = setNames(
  colorRampPalette(c("steelblue3", "darkgoldenrod", "firebrick"))(length(sigma_vals)),
  as.character(sigma_vals)
)

# helper
make_site_plot = function(site_name, q, accum_groups, post_param, Xgroup, sigma_vals,
                          sigma_colors, r0, nu0, mseed, ylim = NULL, 
                          show_legend = FALSE, show_y_axis = TRUE) {
  acc_site = accum_groups %>% filter(group == site_name)
  pred_list = lapply(sigma_vals, function(sig) {
    pred = rar_pred_Kgroup(alphaq = post_param$alpha_out[, q], 
                           thetaq = post_param$theta_out[, q], r0 = r0, nu0 = nu0,
                           nq = 0, kq = 0, sq = Xgroup$n[q], seed = mseed,
                           sigmaq = sig)
    df = data.frame(s = 1:nrow(pred), mean = rowMeans(pred),
                    lo = apply(pred, 1, quantile, probs = 0.025),
                    hi = apply(pred, 1, quantile, probs = 0.975),
                    sigma = as.factor(sig))
    df = rbind(data.frame(s = 0, mean = 0, lo = 0, hi = 0, sigma = as.factor(sig)), 
               df)
    df
  })
  df_pred = do.call(rbind, pred_list)
  df_pred$facet = site_name  # dummy variable for facet strip
  
  p = ggplot() +
    geom_ribbon(data = df_pred, aes(x = s, ymin = lo, ymax = hi, fill = sigma),
                alpha = 0.2) +
    geom_line(data = df_pred, aes(x = s, y = mean, color = sigma), linewidth = 0.9) +
    scale_color_manual(values = sigma_colors, name = expression(sigma)) +
    scale_fill_manual(values  = sigma_colors, name = expression(sigma)) +
    coord_cartesian(ylim = ylim) +
    labs(x = "Sampling effort", y = TeX("$L_{n_p,n_q}$")) +
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

global_ylim = local({
  y_pred_hi = unlist(lapply(sites_info, function(s) {
    unlist(lapply(sigma_vals, function(sig) {
      pred = rar_pred_Kgroup(
        alphaq = post_fungi$alpha_out[, s$q], thetaq = post_fungi$theta_out[, s$q], 
        r0 = r0, nu0 = nu0, nq = 0, kq = 0, sq = Xgroup$n[s$q], seed = mseed, 
        sigmaq = sig
      )
      apply(pred, 1, quantile, probs = 0.975)
    }))
  }))
  y_emp = accum_groups$n_feat
  c(0, max(c(y_emp, y_pred_hi), na.rm = TRUE) * 1.05)
})

# build plots
plot_list = lapply(seq_along(sites_info), function(i) {
  s = sites_info[[i]]
  make_site_plot(site_name = s$site, q = s$q, accum_groups = accum_groups,
                 post_param = post_fungi, Xgroup = Xgroup, sigma_vals = sigma_vals,
                 sigma_colors = sigma_colors, r0 = r0, nu0 = nu0, mseed = mseed,
                 ylim = global_ylim,
                 show_legend  = (i == length(sites_info)),
                 show_y_axis  = (i == 1))           
})

shared_legend = get_legend(plot_list[[length(plot_list)]])
plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] +
  theme(legend.position = "none")

# combine
psigma_group = grid.arrange( 
  arrangeGrob(grobs = lapply(plot_list, ggplotGrob), nrow = 1),
  shared_legend,
  ncol = 2,
  widths = c(0.9, 0.1)
)

# ggsave("fungi_analysis/figures/fungi_sigma_group.png", plot = psigma_group,
#        device = 'png',width = 11, height = 4.3)


################################################################################
#-------------# shared species prediction curves with semi-test #--------------#

occ = gibbs_mosaic_parallel(Y = Y, alpha_post = post_fungi$alpha_out, 
                          theta_post = post_fungi$theta_out, Npost = Npost, n = n, 
                          burnin = burnin, thin = 1, sigma = sigma3, mseed = mseed,
                          n_cores = 34)

pi_post = Map(function(g_list, e_list) {
  Map(function(g_vec, e_vec) g_vec * e_vec, g_list, e_list)
}, occ$gamma, occ$eta)

# helper 
build_shared_data = function(q, p, pair_name, s_aug, B) {
  n_obs_q = Xgroup$n[q] ; n_obs_p = Xgroup$n[p]
  u = n_obs_p - n_obs_q
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
  
  curves = mclapply(1:B, function(b) {
    idx_q = ((b-1)*s_aug_q + 1):(b*s_aug_q)
    Yaugq = matrix(0, nrow = s_aug_q, ncol = max_len)
    for (i in seq_along(idx_q)) {
      pii = pi_post[[q]][[idx_q[i]]]
      Yaugq[i, 1:length(pii)] = rbinom(length(pii), 1, pii)
    }
    Yaugq_all = rbind(Yq_base, Yaugq)
    idx_p = ((b-1)*s_ext + 1):(b*s_ext)
    Yaugp = matrix(0, nrow = s_ext, ncol = max_len)
    for (i in seq_along(idx_p)) {
      piip = pi_post[[p]][[idx_p[i]]]
      Yaugp[i, 1:length(piip)] = rbinom(length(piip), 1, piip)
    }
    Yaugp_all = rbind(Yp_base, Yaugp)
    simplified_accum(n_pq = c(nrow(Yaugp_all), nrow(Yaugq_all)),
                     Yp = Yaugp_all, Yq = Yaugq_all, nperm = 100)
  }, mc.cores = 2)
  
  if (any(sapply(curves, is.null))) stop("simplified_accum returned NULL")
  
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


pairs_list = list(list(q = 15, p = 22, name = "KLA-FTS"),
                  list(q = 16, p = 18, name = "NUU-PRI"),
                  list(q = 15, p = 34, name = "KLA-ZUR"))
s_aug = 30
B = 30

results = lapply(pairs_list, function(pr) {
  cat("Processing pair:", pr$name, "\n")
  build_shared_data(pr$q, pr$p, pr$name, s_aug = s_aug, B = B)
})

# load("/Users/stolffederica/Library/CloudStorage/Dropbox/ActiveResearch/Pex_IBP/code/code_pBB/fungi_analysis/resDCC/res_semitrain4.Rdata")

all_emp = do.call(rbind, lapply(results, `[[`, "emp"))
all_semi = do.call(rbind, lapply(results, `[[`, "semi"))
all_pred = do.call(rbind, lapply(results, `[[`, "pred"))
all_vline = do.call(rbind, lapply(results, `[[`, "vlines"))

pair_levels = sapply(pairs_list, `[[`, "name")
all_emp$pair = factor(all_emp$pair,   pair_levels)
all_semi$pair = factor(all_semi$pair,  pair_levels)

all_pred$pair  = factor(all_pred$pair,  pair_levels)
all_vline$pair = factor(all_vline$pair, pair_levels)

# plot
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
  labs(x = "Sampling effort", y = TeX("$C_{n_p, n_q}$")) +
  theme_minimal() +
  theme(
    strip.background = element_rect(fill = "grey70", color = NA),
    strip.text = element_text(size = 19), axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 17), axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
    )
pall_semi

# ggsave(plot = pall_semi, filename = "fungi_analysis/figures/semitest.png",
#        width = 15, height = 5)
