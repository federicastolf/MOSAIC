functions {
  // log rising factorial: log (a)_n = log Gamma(a+n) - log Gamma(a)
  real log_rising(real a, int n) {
    if (n == 0) return 0;
    return lgamma(a + n) - lgamma(a);
  }
}

data {
  int<lower=1> Q;                 // number of groups/sites
  int<lower=1> H;                 // predictors dimension
  int<lower=1> k;                 // number of observed species
  int<lower=1> n_q[Q];            // sample sizes for site q
  int<lower=0> m[k, Q];           // species frequencies
  matrix[Q, H] X;                 
  vector[H] b0;
  matrix[H, H] B0;
  real<lower=0> a_phi;            // gamma shape for phi
  real<lower=0> b_phi;            // gamma rate for phi
  real<lower=0> r0;               
  real<lower=0> nu0;              
}

parameters {
  vector[H] beta;
  vector<lower=0>[Q] phi;         
}

transformed parameters {
  vector[Q] mu;
  vector[Q] alpha;
  vector[Q] theta;

  for (q in 1:Q) {
    mu[q]    = inv_logit( X[q] * beta );
    theta[q] = phi[q];
    alpha[q] = -mu[q] * phi[q];            
  }
}

model {
  // Priors
  beta ~ multi_normal(b0, B0);
  phi  ~ gamma(a_phi, b_phi);              // shape, rate

  // log marginal likelihood 
  {
    real c = r0 / (nu0 + r0);
    real p = nu0 / (nu0 + r0);

    // log R = sum_q [ log (theta+alpha)_{nq} - log (theta)_{nq} ]
    real logR = 0;
    for (q in 1:Q) {
      real ta = theta[q] + alpha[q];       // = (1-mu_q)*phi_q > 0
      logR += log_rising(ta, n_q[q]) - log_rising(theta[q], n_q[q]);
    }
    real R = exp(logR);

    // log choose(k+r0-1, k) = lgamma(k+r0) - lgamma(k+1) - lgamma(r0)
    target += lgamma(k + r0) - lgamma(k + 1) - lgamma(r0);

    // r0*log c + k*log c
    target += (r0 + k) * log(c);

    // -(k+r0)*log(1 - p*R)
    target += -(k + r0) * log1m(p * R);

    // + sum_q k[ log(-alpha_q) - log (theta_q)_{nq} ]
    for (q in 1:Q) {
      target += k * ( log(-alpha[q]) - log_rising(theta[q], n_q[q]) );
    }

    // + sum_{j=1}^k sum_{q=1}^Q [ log(1-alpha_q)_{mjq-1} + log(theta_q+alpha_q)_{nq-mjq} ]
    if (k > 0) {
      for (j in 1:k) {
        for (q in 1:Q) {
          int mjq = m[j, q];
          target += log_rising(1 - alpha[q], mjq - 1);
          target += log_rising(theta[q] + alpha[q], n_q[q] - mjq);
        }
      }
    }
  }
}

generated quantities {
  vector[Q] alpha_out = alpha;
  vector[Q] theta_out = theta;
}
