split_train_test_new <- function(data, n_train = NULL, n_test = NULL, 
                             train_prop = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  Y <- data$Y
  n <- data$n
  Q <- length(n)
  
  ends_n <- cumsum(n)
  starts_n <- c(1, head(ends_n, -1) + 1)
  
  # recycle scalar inputs
  if (!is.null(n_train) && length(n_train) == 1) n_train <- rep(n_train, Q)
  if (!is.null(n_test)  && length(n_test)  == 1) n_test  <- rep(n_test, Q)
  
  Y_train <- vector("list", Q)
  Y_test <- vector("list", Q)
  Y_semitrain <- vector("list", Q)
  
  n_train_out <- numeric(Q)
  n_test_out <- numeric(Q)
  n_semitrain <- numeric(Q)
  
  for (q in seq_len(Q)) {
    rows_q <- starts_n[q]:ends_n[q]
    Yq <- Y[rows_q, , drop = FALSE]
    
    # define train/test sizes
    if (!is.null(n_train) && !is.null(n_test)) {
      ntrain_q <- n_train[q]
      ntest_q  <- n_test[q]
    } else if (!is.null(n_train) && is.null(n_test)) {
      ntrain_q <- n_train[q]
      ntest_q  <- n[q] - ntrain_q
    } else if (is.null(n_train) && !is.null(n_test)) {
      ntest_q  <- n_test[q]
      ntrain_q <- floor(train_prop * (n[q] - ntest_q))
    } else {
      ntrain_q <- floor(train_prop * n[q])
      ntest_q  <- n[q] - ntrain_q
    }
    
    nsemi_q <- n[q] - ntrain_q - ntest_q
    
    if (nsemi_q < 0) {
      stop(sprintf("For group %d: n_train + n_test > n[q]", q))
    }
    
    idx <- sample(seq_len(n[q]))
    
    idx_tr <- idx[seq_len(ntrain_q)]
    idx_te <- idx[seq(from = ntrain_q + 1, length.out = ntest_q)]
    idx_st <- idx[seq(from = ntrain_q + ntest_q + 1, length.out = nsemi_q)]
    
    Y_train[[q]] <- Yq[idx_tr, , drop = FALSE]
    Y_test[[q]] <- Yq[idx_te, , drop = FALSE]
    Y_semitrain[[q]] <- Yq[idx_st, , drop = FALSE]
    
    n_train_out[q] <- ntrain_q
    n_test_out[q] <- ntest_q
    n_semitrain[q] <- nsemi_q
  }
  
  # stack training data WITHOUT removing all-zero columns
  Y_train_all <- do.call(rbind, Y_train)
  Y_train_all1 <- Y_train_all[, colSums(Y_train_all) > 0, drop = FALSE]
  k <- ncol(Y_train_all1)
  # k = sum(colSums(Y_train_all)>0)
  
  ends_train <- cumsum(n_train_out)
  starts_train <- c(1, head(ends_train, -1) + 1)
  
  m_train <- matrix(0, nrow = k, ncol = Q)
  kq_train <- numeric(Q)
  
  for (q in seq_len(Q)) {
    if (n_train_out[q] > 0) {
      Yq <- Y_train_all1[starts_train[q]:ends_train[q], , drop = FALSE]
      cs <- colSums(Yq)
      m_train[, q] <- cs
      kq_train[q] <- sum(cs > 0)
    }
  }
  
  list(
    Y_train = Y_train,
    Y_test = Y_test,
    Y_semitrain = Y_semitrain,
    n_train = n_train_out,
    n_test = n_test_out,
    n_semitrain = n_semitrain,
    Y_train_all = Y_train_all,
    k_train = k,
    m_train = m_train,
    kq_train = kq_train
  )
}

#------------------------------------------------------------------------------#
#---------------# prediction model based shared-species #----------------------#


d_n = function(n_q, a_q, b_q) {
  exp(lgamma(b_q + n_q) - lgamma(b_q) - (lgamma(b_q + a_q + n_q) - lgamma(b_q + a_q)))
}

# d_n defined above

rar_pred_shared = function(k, k_p, k_q, c_pq, r0, nu0, alpha, theta, n, s_pq, p, q,
                           full = FALSE) {
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
    smax   = max(n_sp, n_sq)
    sp_vec = pmin(1:smax, n_sp)
    sq_vec = pmin(1:smax, n_sq)
    result = matrix(0L, nrow = smax, ncol = B)
  }
  
  for (b in 1:B) {
    alpha_b = alpha[b, ]
    theta_b = theta[b, ]
    
    c0 = 1 - prod(d_n(n, -alpha_b, alpha_b + theta_b))
    
    if (full) {
      d0p_vec = d_n(1:n_sp, -alpha_b[p], theta_b[p] + alpha_b[p] + n[p])  # length n_sp
      d0q_vec = d_n(1:n_sq, -alpha_b[q], theta_b[q] + alpha_b[q] + n[q])  # length n_sq
      omp = 1 - d0p_vec
      omq = 1 - d0q_vec
      prob_both = outer(omp, omq)
      nb_mu     = (r0 + k) * (1 - c0) * prob_both / (r0 / nu0 + c0)
      result[,, b] = array(
        rnbinom(n_sp * n_sq, size = r0 + k,     mu = as.vector(nb_mu))        +
          rbinom (n_sp * n_sq, size = k - k_un,   prob = as.vector(prob_both))     +
          rbinom (n_sp * n_sq, size = k_p - c_pq, prob = as.vector(outer(rep(1, n_sp), omq))) +
          rbinom (n_sp * n_sq, size = k_q - c_pq, prob = as.vector(outer(omp, rep(1, n_sq)))),
        dim = c(n_sp, n_sq)
      )
      
    } else {
      d0p_vec = d_n(sp_vec, -alpha_b[p], theta_b[p] + alpha_b[p] + n[p])  # length smax
      d0q_vec = d_n(sq_vec, -alpha_b[q], theta_b[q] + alpha_b[q] + n[q])  # length smax
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

# volendo si può anche fare una cosa più efficiente riscrivendo la funzione
# sopra e riutilizzando dq e dp
rar_pred_shared_all_pairs = function(k, kq, c_pairs, r0, nu0, alpha, theta, n, s,
                                     full = FALSE) {
  Q       = length(kq)
  pairs   = combn(Q, 2)
  n_pairs = ncol(pairs)
  result        = vector("list", n_pairs)
  names(result) = apply(pairs, 2, function(x) paste0(x[1], "_", x[2]))
  for (i in 1:n_pairs) {
    p = pairs[1, i]
    q = pairs[2, i]
    result[[i]] = rar_pred_shared(k = k, k_p = kq[p], k_q = kq[q],
                                  c_pq = c_pairs[i], r0 = r0, nu0 = nu0,
                                  alpha = alpha, theta = theta, n = n,
                                  s_pq = s[c(p, q)], p = p, q = q,
                                  full = full)
  }
  result   # full=FALSE: list of smax x B;  full=TRUE: list of n_sp x n_sq x B
}


#------------------------------------------------------------------------------#
#---------------- shared-species accumulation: clean version ------------------#
#------------------------------------------------------------------------------#

#==============================================================================#
# 1) FULL 2D accumulation for one ordering
#==============================================================================#
shared_accum_2D <- function(Yp, Yq) {
  np <- nrow(Yp)
  nq <- nrow(Yq)
  
  res <- matrix(0L, nrow = np, ncol = nq)
  
  Ap <- rep(FALSE, ncol(Yp))
  for (i in seq_len(np)) {
    Ap <- Ap | (Yp[i, ] > 0)
    
    Aq <- rep(FALSE, ncol(Yq))
    for (j in seq_len(nq)) {
      Aq <- Aq | (Yq[j, ] > 0)
      res[i, j] <- sum(Ap & Aq)
    }
  }
  
  res
}

#==============================================================================#
# 2) FULL 2D accumulation averaged over permutations
#==============================================================================#
shared_accum_perm_pair_full <- function(Yp, Yq, nperm = 100) {
  np <- nrow(Yp)
  nq <- nrow(Yq)
  
  acc_full <- array(0, dim = c(np, nq, nperm))
  
  for (b in seq_len(nperm)) {
    Ypb <- Yp[sample(np), , drop = FALSE]
    Yqb <- Yq[sample(nq), , drop = FALSE]
    acc_full[, , b] <- shared_accum_2D(Ypb, Yqb)
  }
  
  apply(acc_full, c(1, 2), mean)
}

#==============================================================================#
# 3) SIMPLIFIED diag accumulation for one ordering
#    Requires same number of rows in the two groups
#==============================================================================#
shared_accum_pair_diag <- function(Yp, Yq) {
  np <- nrow(Yp)
  nq <- nrow(Yq)
  
  if (np != nq) {
    stop("shared_accum_pair_diag requires nrow(Yp) = nrow(Yq).")
  }
  if (ncol(Yp) != ncol(Yq)) {
    stop("Yp and Yq must have the same number of columns.")
  }
  
  Ap <- rep(FALSE, ncol(Yp))
  Aq <- rep(FALSE, ncol(Yq))
  out <- integer(np)
  
  for (t in seq_len(np)) {
    Ap <- Ap | (Yp[t, ] > 0)
    Aq <- Aq | (Yq[t, ] > 0)
    out[t] <- sum(Ap & Aq)
  }
  
  out
}

#==============================================================================#
# 4) SIMPLIFIED diag accumulation averaged over permutations
#==============================================================================#
shared_accum_perm_pair_diag <- function(Yp, Yq, nperm = 100) {
  np <- nrow(Yp)
  nq <- nrow(Yq)
  
  if (np != nq) {
    stop("shared_accum_perm_pair_diag requires nrow(Yp) = nrow(Yq).")
  }
  
  acc <- matrix(0, nrow = np, ncol = nperm)
  
  for (b in seq_len(nperm)) {
    Ypb <- Yp[sample(np), , drop = FALSE]
    Yqb <- Yq[sample(nq), , drop = FALSE]
    acc[, b] <- shared_accum_pair_diag(Ypb, Yqb)
  }
  
  rowMeans(acc)
}

#==============================================================================#
# 5) SIMPLIFIED diag accumulation for anchored continuation
#    Yp0, Yq0 = anchor already observed
#    Yp1, Yq1 = future block to be accumulated
#==============================================================================#
shared_accum_pair_anchor_diag <- function(Yp0, Yq0, Yp1, Yq1,
                                          include_anchor = TRUE) {
  K <- ncol(Yp0)
  
  if (ncol(Yq0) != K || ncol(Yp1) != K || ncol(Yq1) != K) {
    stop("All matrices must have the same number of columns.")
  }
  
  np1 <- nrow(Yp1)
  nq1 <- nrow(Yq1)
  
  if (np1 != nq1) {
    stop("shared_accum_pair_anchor_diag requires nrow(Yp1) = nrow(Yq1).")
  }
  
  Ap <- if (nrow(Yp0) > 0) colSums(Yp0 > 0) > 0 else rep(FALSE, K)
  Aq <- if (nrow(Yq0) > 0) colSums(Yq0 > 0) > 0 else rep(FALSE, K)
  
  out <- integer(np1 + 1)
  out[1] <- sum(Ap & Aq)
  
  for (t in seq_len(np1)) {
    Ap <- Ap | (Yp1[t, ] > 0)
    Aq <- Aq | (Yq1[t, ] > 0)
    out[t + 1] <- sum(Ap & Aq)
  }
  
  if (!include_anchor) out <- out[-1]
  out
}

#==============================================================================#
# 6) SIMPLIFIED anchored continuation averaged over permutations
#==============================================================================#
shared_accum_perm_pair_anchor_diag <- function(Yp0, Yq0, Yp1, Yq1,
                                               nperm = 100,
                                               include_anchor = TRUE) {
  np1 <- nrow(Yp1)
  nq1 <- nrow(Yq1)
  
  if (np1 != nq1) {
    stop("shared_accum_perm_pair_anchor_diag requires nrow(Yp1) = nrow(Yq1).")
  }
  
  len <- np1 + as.integer(include_anchor)
  acc <- matrix(0, nrow = len, ncol = nperm)
  
  for (b in seq_len(nperm)) {
    ord_p <- if (np1 > 0) sample(np1) else integer(0)
    ord_q <- if (nq1 > 0) sample(nq1) else integer(0)
    
    Yp1b <- if (np1 > 0) Yp1[ord_p, , drop = FALSE] else Yp1
    Yq1b <- if (nq1 > 0) Yq1[ord_q, , drop = FALSE] else Yq1
    
    acc[, b] <- shared_accum_pair_anchor_diag(
      Yp0 = Yp0, Yq0 = Yq0,
      Yp1 = Yp1b, Yq1 = Yq1b,
      include_anchor = include_anchor
    )
  }
  
  rowMeans(acc)
}

#==============================================================================#
# 7) All-pairs wrappers
#==============================================================================#
shared_accum_perm_allpairs_full <- function(Y_list, nperm = 100) {
  Q <- length(Y_list)
  pairs <- combn(Q, 2, simplify = FALSE)
  
  out <- list()
  for (pq in pairs) {
    p <- pq[1]
    q <- pq[2]
    key <- paste0(p, "_", q)
    out[[key]] <- shared_accum_perm_pair_full(Y_list[[p]], Y_list[[q]], nperm = nperm)
  }
  
  out
}

shared_accum_perm_allpairs_diag <- function(Y_list, nperm = 100) {
  Q <- length(Y_list)
  pairs <- combn(Q, 2, simplify = FALSE)
  
  out <- list()
  for (pq in pairs) {
    p <- pq[1]
    q <- pq[2]
    key <- paste0(p, "_", q)
    out[[key]] <- shared_accum_perm_pair_diag(Y_list[[p]], Y_list[[q]], nperm = nperm)
  }
  
  out
}

#==============================================================================#
# 8) Observed shared species for all pairs
#==============================================================================#
shared_species_pairs <- function(Y_list) {
  Q <- length(Y_list)
  A <- lapply(Y_list, function(Y) colSums(Y > 0) > 0)
  pairs <- combn(Q, 2)
  
  out <- apply(pairs, 2, function(pq) sum(A[[pq[1]]] & A[[pq[2]]]))
  names(out) <- apply(pairs, 2, function(pq) paste0(pq[1], "_", pq[2]))
  out
}

