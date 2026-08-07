
library(vegan)

################################################################################
#----------------# group-specific rarefaction curve #--------------------------#
################################################################################

#-----------# empirical #---------------#

#' Computes empirical rarefaction curves for each group using vegan::specaccum
accum_by_group = function(Y, n, method = "random", permutations, lab_group, 
                          mseed) {
  set.seed(mseed)
  Q = length(n)
  out = vector("list", Q)
  ends_n = cumsum(n)
  starts_n = c(1, head(ends_n, -1) + 1)
  for (q in 1:Q) {
    Yq = Y[starts_n[q]:ends_n[q], , drop = FALSE]
    acc = specaccum(Yq, method = method, permutations = permutations)
    richness = acc$richness
    out[[q]] = data.frame(x = acc$sites, n_feat = richness, group = lab_group[q])
  }
  do.call(rbind, out)
}

#-----------# model-prediction #---------------#

# building block for model-based rarefaction predictions.
d_n = function(n_q, a_q, b_q) {
  exp(lgamma(b_q + n_q) - lgamma(b_q) - (lgamma(b_q + a_q + n_q) - lgamma(b_q + a_q)))
}


# model-based prediction for group-specific rarefaction curves

rar_pred_Kgroup = function(alphaq, thetaq, r0, nu0, nq, kq, sq, seed, 
                           p0_all, k, sigmaq = NULL) {
  set.seed(seed)
  B = length(alphaq)
  s_vec = 1:sq
  ksn_chain = matrix(0, nrow = length(s_vec), ncol = B)
  for (b in 1:B) {
    alphab = alphaq[b]
    thetab = thetaq[b]
    dn = d_n(n = nq, a = -alphab, b = alphab + thetab)
    if (is.null(sigmaq)){
      ds = d_n(n = s_vec, a = -alphab, b = alphab + thetab + nq)
    }
    else{
      ds = d_n(n = s_vec, a = -alphab, b =  sigmaq*(alphab + thetab) + nq)
    }
    # Species never recorded at any site
    size1 = r0 + k
    mu = size1 * (1 - ds) * p0_all[b] / (r0 / nu0 + 1 - p0_all[b])
    u1 = rnbinom(length(s_vec), size = size1, mu = mu)
    # Species recorded at some other site but not at site q
    u2 = rbinom(length(s_vec), size = k - kq, prob = 1 - ds)
    ksn_chain[, b] = u1 + u2
  }
  return(ksn_chain + kq)
}


# rar_pred_Kgroup = function(alphaq, thetaq, r0, nu0, nq, kq, sq, seed, 
#                            sigmaq = NULL) {
#   set.seed(seed)
#   B = length(alphaq)
#   s_vec = 1:sq
#   ksn_chain = matrix(0, nrow = length(s_vec), ncol = B)
#   for (b in 1:B) {
#     alphab = alphaq[b]
#     thetab = thetaq[b]
#     dn = d_n(n = nq, a = -alphab, b = alphab + thetab)
#     if (is.null(sigmaq)){
#       ds = d_n(n = s_vec, a = -alphab, b = alphab + thetab + nq)
#     }
#     else{
#       ds = d_n(n = s_vec, a = -alphab, b =  sigmaq*(alphab + thetab) + nq)
#     }
#     size = r0 + kq
#     mu = size * (1 - ds) * dn / (r0 / nu0 + 1 - dn)
#     ksn_chain[, b] = rnbinom(length(s_vec), size = size, mu = mu)
#   }
#   return(ksn_chain + kq)
# }


################################################################################
#----------------------# shared-species rarefaction #--------------------------#
################################################################################

#------------------------------# empirical #-----------------------------------#

# Computes the empirical simplified shared-species accumulation curve 
# (diagonal of the 2D surface) 
simplified_accum = function(n_pq, Yp, Yq, nperm = 100){
  a = min(n_pq[1], n_pq[2])
  acc_simpl = matrix(0, nrow = a, ncol = nperm)
  for (b in 1:nperm) {
    ord_p = sample(n_pq[1])
    ord_q = sample(n_pq[2])
    acc_simpl[, b] = shared_accum_simplified_direct(Yp[ord_p, , drop = FALSE], 
                                                     Yq[ord_q, , drop = FALSE])[1:a]
  }
  return(rowMeans(acc_simpl))
}

# Computes the empirical simplified shared-species accumulation curve on 
# test data, anchored to species already observed in the training set.
simplified_accum_anchor = function(ntest_pq, Yp_test, Yq_test, Yp_train, 
                                   Yq_train, nperm = 100){
  a = min(ntest_pq[1], ntest_pq[2])
  acc_simpl = matrix(0, nrow = a, ncol = nperm)
  for (b in 1:nperm) {
    ord_p = sample(ntest_pq[1])
    ord_q = sample(ntest_pq[2])
    acc_simpl[, b] = shared_accum_simplified_train_anchor_direct(
      Yp_train, Yq_train,Yp_test[ord_p, , drop = FALSE], 
      Yq_test[ord_q, , drop = FALSE]
    )[1:a]
  }
  return(rowMeans(acc_simpl))
}


# Computes the full 2D shared-species accumulation surface 
surface_accum = function(n_pq, Yp, Yq, nperm = 100) {
  acc_full = array(0, dim = c(n_pq[1], n_pq[2], nperm))
  for (b in 1:nperm) {
    ord_p = sample(n_pq[1])
    ord_q = sample(n_pq[2])
    acc_full[, , b] = shared_accum_2D(Yp[ord_p, , drop = FALSE], 
                                      Yq[ord_q, , drop = FALSE])
  }
  return(apply(acc_full, c(1, 2), mean))
}

# Computes the full 2D shared-species accumulation surface on test data,
# anchored to species already observed in the training set.
surface_accum_anchor = function(ntest_pq, Yp_test, Yq_test, Yp_train, Yq_train,
                                nperm = 100) {
  acc_full = array(0, dim = c(ntest_pq[1], ntest_pq[2], nperm))
  for (b in 1:nperm) {
    ord_p = sample(ntest_pq[1])
    ord_q = sample(ntest_pq[2])
    acc_full[, , b] = shared_accum_2D_train_anchor(Yp_train, Yq_train,
                                                    Yp_test[ord_p, , drop = FALSE],
                                                    Yq_test[ord_q, , drop = FALSE])
  }
  return(apply(acc_full, c(1, 2), mean))
}


# helper for surface_accum
shared_accum_2D = function(Yp, Yq) {
  np = nrow(Yp); nq = nrow(Yq)
  res = matrix(0, np, nq)
  Ap_seen = rep(FALSE, ncol(Yp))
  for (sp in 1:np) {
    Ap_seen = Ap_seen | (Yp[sp, ] > 0)
    Aq_seen = rep(FALSE, ncol(Yq))
    for (sq in 1:nq) {
      Aq_seen = Aq_seen | (Yq[sq, ] > 0)
      res[sp, sq] = sum(Ap_seen & Aq_seen)
    }
  }
  res
}

# helper for simplified_accum
shared_accum_simplified_direct = function(Yp, Yq) {
  np = nrow(Yp); nq = nrow(Yq)
  a = max(np, nq)
  p_idx = pmin(seq_len(a), np) ; q_idx = pmin(seq_len(a), nq)
  out = numeric(a)
  Ap_seen = rep(FALSE, ncol(Yp)) ; Aq_seen = rep(FALSE, ncol(Yq))
  p_prev = 0 ; q_prev = 0
  for (u in seq_len(a)) {
    p_now = p_idx[u] ; q_now = q_idx[u]
    if (p_now > p_prev) {
      for (sp in (p_prev + 1):p_now) {
        Ap_seen = Ap_seen | (Yp[sp, ] > 0)
      }
    }
    if (q_now > q_prev) {
      for (sq in (q_prev + 1):q_now) {
        Aq_seen = Aq_seen | (Yq[sq, ] > 0)
      }
    }
    out[u] = sum(Ap_seen & Aq_seen)
    p_prev = p_now ; q_prev = q_now
  }
  out
}

# helper for surface_accum_anchor
shared_accum_2D_train_anchor = function(Yp_train, Yq_train, Yp_test, Yq_test) {
  np_te = nrow(Yp_test) ; nq_te = nrow(Yq_test)
  Ap_seen = if (nrow(Yp_train) > 0) colSums(Yp_train > 0) > 0 else 
    rep(FALSE, ncol(Yp_test))
  Aq_seen = if (nrow(Yq_train) > 0) colSums(Yq_train > 0) > 0 else 
    rep(FALSE, ncol(Yq_test))
  res = matrix(0, np_te, nq_te)
  Ap_curr = Ap_seen
  for (sp in 1:np_te) {
    Ap_curr = Ap_curr | (Yp_test[sp, ] > 0)
    Aq_curr = Aq_seen
    for (sq in 1:nq_te) {
      Aq_curr = Aq_curr | (Yq_test[sq, ] > 0)
      res[sp, sq] = sum(Ap_curr & Aq_curr)
    }
  }
  res
}

# helper for simplified_accum_anchor
shared_accum_simplified_train_anchor_direct = function(Yp_train, Yq_train, 
                                                        Yp_test, Yq_test) {
  np_te = nrow(Yp_test) ; nq_te = nrow(Yq_test)
  a = max(np_te, nq_te)
  p_idx = pmin(seq_len(a), np_te) ; q_idx = pmin(seq_len(a), nq_te)
  Ap_seen <- if (nrow(Yp_train) > 0) colSums(Yp_train > 0) > 0 else 
    rep(FALSE, ncol(Yp_test))
  Aq_seen <- if (nrow(Yq_train) > 0) colSums(Yq_train > 0) > 0 else 
    rep(FALSE, ncol(Yq_test))
  out = numeric(a)
  Ap_curr = Ap_seen ; Aq_curr = Aq_seen
  p_prev = 0 ; q_prev = 0
  for (u in seq_len(a)) {
    p_now = p_idx[u] ; q_now = q_idx[u]
    if (p_now > p_prev) {
      for (sp in (p_prev + 1):p_now) {
        Ap_curr <- Ap_curr | (Yp_test[sp, ] > 0)
      }
    }
    if (q_now > q_prev) {
      for (sq in (q_prev + 1):q_now) {
        Aq_curr <- Aq_curr | (Yq_test[sq, ] > 0)
      }
    }
    out[u] = sum(Ap_curr & Aq_curr)
    p_prev = p_now ; q_prev = q_now
  }
  out
}


#' Counts the number of species shared between each pair of groups in Y_list,
#' returning a named vector with keys of the form "p_q".
shared_species_pairs = function(Y_list, labels = NULL) {
  Q = length(Y_list)
  if (!is.list(Y_list) || Q < 2) {
    stop("Y_list must be a list with at least two matrices.")
  }
  if (!is.null(labels) && length(labels) != Q) {
    stop("labels must have length Q.")
  }
  A_list = lapply(Y_list, function(Y) colSums(Y > 0) > 0)
  out = numeric(Q * (Q - 1) / 2)
  names(out) = character(length(out))
  idx = 1
  for (p in 1:(Q - 1)) {
    for (q in (p + 1):Q) {
      out[idx] = sum(A_list[[p]] & A_list[[q]])
      names(out)[idx] = if (!is.null(labels)) paste0(labels[p], "_", 
                                                     labels[q]) else paste0(p, "_", q)
      idx = idx + 1
    }
  }
  out
}

#-------------------------------# model based #--------------------------------#

# Model-base predictive shared-species 2D surface (full=T) or diagonal (full=F) 
# for a pair (p,q), returning an array of draws.
rar_pred_shared = function(k, k_p, k_q, c_pq, r0, nu0, alpha, theta, n, s_pq, p, q,
                           full = FALSE, sigma_pq = NULL) {
  # alpha, theta: B x Q matrix
  # s_pq: vector of length 2, n: vector of length Q
  stopifnot(ncol(alpha) == ncol(theta), ncol(alpha) == length(n))
  stopifnot(k >= k_p + k_q - c_pq, k_p >= c_pq, k_q >= c_pq)
  B = nrow(alpha)
  k_un = k_p + k_q - c_pq
  n_sp = s_pq[1]
  n_sq = s_pq[2]
  if (full) {
    result = array(0L, dim = c(n_sp, n_sq, B))
  } else {
    smax = min(n_sp, n_sq)
    sp_vec = sq_vec = 1:smax
    result = matrix(0L, nrow = smax, ncol = B)
  }
  for (b in 1:B) {
    alpha_b = alpha[b, ]
    theta_b = theta[b, ]
    c0 = 1 - prod(d_n(n, -alpha_b, alpha_b + theta_b))
    if (full) {
      if (is.null(sigma_pq)){
        d0p_vec = d_n(1:n_sp, -alpha_b[p], theta_b[p] + alpha_b[p] + n[p])  # length n_sp
        d0q_vec = d_n(1:n_sq, -alpha_b[q], theta_b[q] + alpha_b[q] + n[q])  # length n_sq
      }
      else{
        d0p_vec = d_n(1:n_sp, -alpha_b[p], sigma_pq[1]*(theta_b[p] + alpha_b[p]) 
                      + n[p])  # length n_sp
        d0q_vec = d_n(1:n_sq, -alpha_b[q], sigma_pq[2]*(theta_b[q] + alpha_b[q]) 
                      + n[q])  # length n_sq
      }
      omp = 1 - d0p_vec
      omq = 1 - d0q_vec
      prob_both = outer(omp, omq)
      nb_mu = (r0 + k) * (1 - c0) * prob_both / (r0 / nu0 + c0)
      result[,, b] = array(
        rnbinom(n_sp*n_sq, size = r0 + k, mu = as.vector(nb_mu)) +
          rbinom(n_sp*n_sq, size = k - k_un,   prob = as.vector(prob_both)) +
          rbinom(n_sp*n_sq, size = k_p - c_pq, prob = as.vector(outer(rep(1, n_sp), omq))) +
          rbinom(n_sp*n_sq, size = k_q - c_pq, prob = as.vector(outer(omp, rep(1, n_sq)))),
        dim = c(n_sp, n_sq)
      )
      
    } else {
      if (is.null(sigma_pq)){
        d0p_vec = d_n(sp_vec, -alpha_b[p], theta_b[p] + alpha_b[p] + n[p])  # length smax
        d0q_vec = d_n(sq_vec, -alpha_b[q], theta_b[q] + alpha_b[q] + n[q])  # length smax
      }
      else{
        d0p_vec = d_n(sp_vec, -alpha_b[p], sigma_pq[1]*(theta_b[p] + alpha_b[p])
                      + n[p])  # length smax
        d0q_vec = d_n(sq_vec, -alpha_b[q], sigma_pq[2]*(theta_b[q] + alpha_b[q])
                      + n[q])  # length smax
      }
      omp = 1 - d0p_vec
      omq = 1 - d0q_vec
      nb_mu = (r0 + k) * (1 - c0) * omp * omq / (r0 / nu0 + c0)
      result[, b] = rnbinom(smax, size = r0 + k, mu = nb_mu) + 
        rbinom (smax, size = k - k_un, prob = omp * omq) +
        rbinom (smax, size = k_p - c_pq, prob = omq) +
        rbinom (smax, size = k_q - c_pq, prob = omp)
    }
  }
  result + c_pq   # full: n_sp x n_sq x B,  diagonal: smax x B
}


#' Computes both the empirical anchored surface and the posterior predictive 
#' surface for a single pair (p, q), returning tidy dataframes for each.
compute_surface_pair = function(p, q, tt_split, c_train, post_param, r0, nu0,
                                pair_key, nperm = 100) {
  # empirical
  emp_mat = surface_accum_anchor(ntest_pq = tt_split$n_test[c(p, q)],
                                 Yp_test = tt_split$Y_test[[p]], 
                                 Yq_test = tt_split$Y_test[[q]], 
                                 Yp_train = tt_split$Y_train[[p]], 
                                 Yq_train = tt_split$Y_train[[q]], nperm = nperm)
  dfemp = expand.grid(s = seq_len(nrow(emp_mat)), t = seq_len(ncol(emp_mat)))
  dfemp$shared = c(emp_mat)
  
  # predicted
  pred_arr = rar_pred_shared(k = tt_split$k_train, k_p = tt_split$kq_train[p], 
                             k_q = tt_split$kq_train[q], c_pq = c_train[pair_key], 
                             r0 = r0, nu0 = nu0, alpha = post_param$alpha_out, 
                             theta = post_param$theta_out, n = tt_split$n_train, 
                             s_pq = tt_split$n_test[c(p, q)], p = p, q = q, 
                             full = TRUE)
  dfpred = expand.grid(s = seq_len(nrow(emp_mat)), t = seq_len(ncol(emp_mat)))
  dfpred$mean  = c(apply(pred_arr, c(1,2), mean))
  dfpred$lower = c(apply(pred_arr, c(1,2), quantile, probs = 0.025))
  dfpred$upper = c(apply(pred_arr, c(1,2), quantile, probs = 0.975))
  list(dfemp = dfemp, dfpred = dfpred)
}
