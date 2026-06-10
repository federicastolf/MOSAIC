
#---# generate synthetic data #---#

synt_dataH = function(n, Q, beta, a_phi, b_phi, r0, nu0 = 3000, mseed = 123,
                      tau   = NULL,    # if not NULL: misspecified model with random intercept sd=tau
                      N     = NULL,    # if not NULL: use fixed N instead of sampling it
                      phi   = NULL) {  # if not NULL: use fixed phi instead of sampling it
  set.seed(mseed)
  H = length(beta)
  X = cbind(1, matrix(rnorm(Q*(H - 1)), nrow = Q))
  
  if (is.null(N))   N   = rnbinom(1, size = r0, mu = nu0)
  if (is.null(phi)) phi = rgamma(Q, a_phi, b_phi)
  
  eta = X %*% beta
  if (!is.null(tau)) eta = eta + rnorm(Q, sd = tau)
  mu = plogis(eta)
  
  Y = vector("list", Q)
  pi = matrix(0, Q, N)
  for (q in 1:Q) {
    Y[[q]] = matrix(0, n[q], N)
    pi[q, ] = rbeta(N, mu[q]*phi[q], (1 - mu[q])*phi[q])
    for (i in 1:n[q])
      Y[[q]][i, ] = rbinom(N, 1, pi[q, ])
  }
  
  Yall = do.call(rbind, Y)
  Yall1 = Yall[, colSums(Yall) > 0, drop = FALSE]
  k = ncol(Yall1)
  
  ends_n = cumsum(n)
  starts_n = c(1L, head(ends_n, -1L) + 1L)
  m = matrix(0, nrow = k, ncol = Q)
  for (q in 1:Q) {
    m[, q] = colSums(Yall1[starts_n[q]:ends_n[q], ])
  }
  list(Y = Yall, X = X, beta = beta, N = N, phi = phi, pi = pi,
       n = n, k = k, m = m, Ylist = Y)
}


#---# Splits a dataset into train/test subsets per group #----#
split_train_test = function(data, train_prop = 0.7,  test_n = NULL, mseed) {
  set.seed(mseed)
  Y = data$Y
  n = data$n
  Q = length(n)
  
  ends_n = cumsum(n)
  starts_n = c(1, head(ends_n, -1) + 1)
  Y_train = list()
  Y_test = list()
  n_train = n_test = numeric(Q)
  
  for (q in 1:Q) {
    rows_q = starts_n[q]:ends_n[q]
    Yq = Y[rows_q, , drop = FALSE]
    
    if (!is.null(test_n)) {
      ntest = test_n
      ntrain = n[q] - ntest
    } else {
      ntrain = floor(train_prop*n[q])
      ntest = n[q] - ntrain
    }
    idx = sample(seq_len(n[q]))
    idx_tr = idx[1:ntrain]
    idx_te = idx[(ntrain + 1):n[q]]
    Y_train[[q]] = Yq[idx_tr, , drop = FALSE]
    Y_test[[q]] = Yq[idx_te, , drop = FALSE]
    n_train[q] = ntrain
    n_test[q] = ntest
  }
  
  # stack training data
  Y_train_all = do.call(rbind, Y_train)
  Y_train_all = Y_train_all[, colSums(Y_train_all != 0) > 0, drop = FALSE]
  k = ncol(Y_train_all)
  ends_train = cumsum(n_train)
  starts_train = c(1, head(ends_train, -1) + 1)
  m_train = matrix(0, nrow = k, ncol = Q) # frequencies matrix
  kq_train = rep(0, Q)
  for (q in 1:Q) {
    Yq = Y_train_all[starts_train[q]:ends_train[q], , drop = FALSE]
    cs = colSums(Yq)
    m_train[, q] = cs
    kq_train[q] = sum(cs>0)
  }
  
  list("Y_train" = Y_train, "Y_test" = Y_test, "n_train" = n_train, 
       "k_train" = k, "n_test" = n_test, "Y_train_all" = Y_train_all, 
       "m_train" = m_train, "kq_train" = kq_train)
}

#----# Converts Y into a list of per-group matrices #----#
convertY_tolist = function(Y, n){
  Q = length(n)
  ends_n = cumsum(n)
  starts_n = c(1, head(ends_n, -1) + 1)
  Y_list = list()
  for (q in 1:Q) {
    rows_q = starts_n[q]:ends_n[q]
    Y_list[[q]]= Y[rows_q, , drop = FALSE]
  }
  return(Y_list)
}

#---# Counts the number of columns with at least one 1 per group #---#
compute_Kq = function(Y_list, labels = NULL) {
  out = sapply(Y_list, function(Y) sum(colSums(Y) > 0))
  if (!is.null(labels)) names(out) = labels
  out
}
