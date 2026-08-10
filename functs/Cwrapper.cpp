#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
void update_WT_inplace(arma::mat& W, arma::mat& T, 
                        const arma::vec& eta, const arma::vec& gamma,
                        int N, int N_prev) {
  
  int nq = W.n_rows;
  int N_update = std::min(N, N_prev);
  
  for(int j = 0; j < N_update; j++) {
    
    double prob1 = (1.0 - eta(j)) * (1.0 - gamma(j));  // (0,0)
    double prob2 = (1.0 - eta(j)) * gamma(j);          // (0,1)
    double prob3 = eta(j) * (1.0 - gamma(j));          // (1,0)
    double prob_sum = prob1 + prob2 + prob3;
    
    Rcpp::NumericVector probs = Rcpp::NumericVector::create(
      prob1 / prob_sum, 
      prob2 / prob_sum, 
      prob3 / prob_sum
    );
    
    std::vector<int> zero_rows;
    for(int i = 0; i < nq; i++) {
      if(W(i, j) == 0.0 && T(i, j) == 0.0) {
        zero_rows.push_back(i);
      }
    }
    
    if(zero_rows.size() == 0) continue;
    
    // Sample from multinomial
    Rcpp::IntegerVector draws = Rcpp::sample(3, zero_rows.size(), true, probs);
    
    
    for(size_t idx = 0; idx < zero_rows.size(); idx++) {
      int i = zero_rows[idx];
      
      if(draws[idx] == 1) {
        W(i, j) = 0.0;
        T(i, j) = 0.0;
      } else if(draws[idx] == 2) {
        W(i, j) = 0.0;
        T(i, j) = 1.0;
      } else {
        W(i, j) = 1.0;
        T(i, j) = 0.0;
      }
    }
  }
  
  if(N > N_prev) {
    for(int j = N_prev; j < N; j++) {
      W.col(j).zeros();
      T.col(j).zeros();
    }
  }
}


// [[Rcpp::export]]
double BetaDiversity_vec(const arma::vec& etap, const arma::vec& etaq) {
  // Both vectors should have the same length for a given iteration
  if(etap.n_elem != etaq.n_elem) {
    Rcpp::stop("etap and etaq must have the same length");
  }
  
  double num = arma::dot(etap, etaq);
  double d1 = arma::dot(etap, etap);
  double d2 = arma::dot(etaq, etaq);
  
  if(d1 < 1e-10 || d2 < 1e-10) {
    return 1.0; // Maximum dissimilarity if either vector is all zeros
  }
  
  double corr = num / std::sqrt(d1 * d2);
  return 1.0 - corr;
}

// [[Rcpp::export]]
arma::mat BetaDiversity_Mat(List eta) {
  int ng = eta.size();
  
  List eta0 = as<List>(eta[0]);
  int Nsampl = eta0.size();
  
  int n_pairs = ng * (ng - 1) / 2;
  arma::mat betaAll(n_pairs, Nsampl);
  
  int id = 0;
  for (int p = 0; p < ng - 1; ++p) {
    List eta_p = as<List>(eta[p]);
    
    for (int q = p + 1; q < ng; ++q) {
      List eta_q = as<List>(eta[q]);
      
      for(int s = 0; s < Nsampl; ++s) {
        arma::vec vec_p = as<arma::vec>(eta_p[s]);
        arma::vec vec_q = as<arma::vec>(eta_q[s]);
        
        betaAll(id, s) = BetaDiversity_vec(vec_p, vec_q);
      }
      id++;
    }
  }
  
  return betaAll;
}

