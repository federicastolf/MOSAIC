library(Rcpp)
library(RcppArmadillo)
library(parallel)

sourceCpp("functs/Cwrapper.cpp")


#---------------------# gibbs-sampler for MOSAIC #----------------------------#

# alpha_post and theta_post from the output of stan

gibbs_mosaic = function(Y, alpha_post, theta_post, Npost, n, burnin,  thin, 
                        sigma, mseed){
  set.seed(mseed)
  Q = ncol(alpha_post)
  n_iter = length(Npost)
  
  thinning = seq(burnin+1, n_iter, by = thin)
  eta = gamma = vector("list", Q)
  
  ends_n = cumsum(n)
  starts_n = c(1, head(ends_n, -1) + 1)
  
  for(q in 1:Q){
    cat(paste0("group: ", q, "\n"))
    
    Yq = Y[starts_n[q]:ends_n[q], ]
    alphaq = alpha_post[,q]
    thetaq = theta_post[,q]
    sigmaq = sigma[q]
  
    resYq = GibbsDA(Yq, n_iter, Npost, alphaq, thetaq, sigmaq)
    
    eta[[q]] = resYq$eta[thinning]
    gamma[[q]] = resYq$gamma[thinning]
  }
  return(list("gamma" = gamma, "eta" = eta, "sigma" = sigma))
}


# Data-augmentation steps

GibbsDA = function(Yq, n_iter, Npost, alphaq, thetaq, sigmaq){
  nq = nrow(Yq)
  k = ncol(Yq)
  
  # Initialization 
  eta_samples = gamma_samples = vector("list", n_iter)
  N1 = Npost[1]
  eta = rep(0.2, N1)
  eta_samples[[1]] = eta
  gamma = rep(0.3, N1)
  gamma_samples[[1]] = gamma
  
  # Pre-allocate W and Ta at maximum size
  N_max = max(Npost)
  W = Ta = matrix(0, nrow = nq, ncol = N_max)
  
  for (iter in 2:n_iter){
    
    alpha = alphaq[iter]
    theta = thetaq[iter]
    N = Npost[iter]
    N_prev = Npost[iter - 1]
    
    # columns for which Yq==1
    W[, 1:k] = Yq
    Ta[, 1:k] = Yq
    
    W[, (k+1):N_max] = 0
    Ta[, (k+1):N_max] = 0
    
    # update W and Ta for columns 1:min(N, N_prev) for which Yq==0
    # for columns beyond N_prev W=T=0
    update_WT_inplace(W, Ta, eta, gamma, N, N_prev)
    
    # update eta and gamma
    b1 = sigmaq*(alpha + theta)
    b2 = (1 - sigmaq)*(alpha + theta)
    wsum = colSums(W[, 1:N, drop=FALSE])
    tsum = colSums(Ta[, 1:N, drop=FALSE])
    eta = rbeta(N, -alpha + wsum, b1 + nq - wsum) 
    gamma = rbeta(N, -alpha + b1 + tsum, b2 + nq - tsum)
    
    eta_samples[[iter]] = eta
    gamma_samples[[iter]] = gamma
    
    if(iter %% 1000 == 0) {
      cat(paste0("iteration: ", iter, "\n"))
    }
  }
  return(list("eta" = eta_samples, "gamma" = gamma_samples))
}


#----------# posterior for N #----------------#

genNNegBinom = function(postfit, n_q, r0, nu0, k, mseed) {
  set.seed(mseed)
  alpha = postfit$alpha_out
  theta = postfit$theta_out
  B = nrow(alpha)
  Q = ncol(alpha)
  nq_mat = matrix(rep(n_q, each = B), nrow = B, ncol = Q)
  log_ratio_mat = (lgamma(alpha + theta + nq_mat) - lgamma(alpha + theta)) -
    (lgamma(theta + nq_mat) - lgamma(theta))
  log_Pi = rowSums(log_ratio_mat)
  Pi = exp(log_Pi)
  r = r0 + k
  denom = (r0/nu0) + 1 - Pi
  mu = r*Pi/denom
  Npost = k + rnbinom(B, size = r, mu = mu)
  return(Npost)
}


# parallel version 

gibbs_mosaic_parallel = function(Y, alpha_post, theta_post, Npost, n, burnin,
                               thin, sigma, mseed, n_cores = NULL){
  
  Q = ncol(alpha_post)
  n_iter = length(Npost)
  thinning = seq(burnin + 1, n_iter, by = thin)
  ends_n = cumsum(n)
  starts_n = c(1, head(ends_n, -1) + 1)
  
  if (is.null(n_cores)) n_cores = min(Q, detectCores() - 1L)
  n_cores = max(1L, n_cores)
  
  results = mclapply(
    seq_len(Q),
    FUN = function(q) {
      set.seed(mseed + q)
      Yq = Y[starts_n[q]:ends_n[q], , drop = FALSE]
      alphaq = alpha_post[, q]
      thetaq = theta_post[, q]
      sigmaq = sigma[q]
      resYq = GibbsDA(Yq, n_iter, Npost, alphaq, thetaq, sigmaq)
      list(eta = resYq$eta[thinning], gamma = resYq$gamma[thinning])
    },
    mc.cores = n_cores, mc.preschedule = FALSE   
  )

  failed = vapply(results, inherits, logical(1), "try-error")
  if (any(failed)) stop("Workers failed for groups: ", paste(which(failed), collapse = ", "))
  eta = lapply(results, `[[`, "eta")
  gamma = lapply(results, `[[`, "gamma")
  return(list("gamma" = gamma, "eta" = eta, "sigma" = sigma))
}