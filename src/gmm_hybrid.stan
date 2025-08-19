data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points for the main group
  array[N] real y;         // observations for the main group
  int<lower=0> N_neg;      // number of data points for the separate group
  array[N_neg] real y_neg; // observations for the negative group
  real y_max;
  real y_min;
}
 
parameters {
  simplex[K] theta;          // mixing proportions
  positive_ordered[K] mu;    // locations of mixture components
  vector<lower=0>[K] sigma;  // scales of mixture components
}

model {
  
  vector[K] log_theta = log(theta);  // cache log calculation

  // Priors
  sigma ~ uniform(0, 10);
  mu ~ uniform(y_min, y_max);
  

  // Main group data
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    target += log_sum_exp(lps);
  }

  // Separate group data modeled by the lowest Gaussian component only
  for (n_neg in 1:N_neg) {
    target += normal_lpdf(y_neg[n_neg] | mu[1], sigma[1]);
  }
}


generated quantities {
  vector[N + N_neg] log_likelihood_all;
  real total_log_likelihood;
  int idx = 1;
  matrix[N, K] prob_y;       // Probabilities for main group
  matrix[N_neg, K] prob_y_neg; // Probabilities for negative group

  total_log_likelihood = 0;

  // Log likelihood for the main group data
  for (n in 1:N) {
    vector[K] lps = log(theta);
    for (k in 1:K) {
      lps[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    log_likelihood_all[idx] = log_sum_exp(lps);
    total_log_likelihood += log_likelihood_all[idx];

    // Calculate probabilities of each observation belonging to each component
    for (k in 1:K) {
      prob_y[n, k] = exp(lps[k] - log_sum_exp(lps));
    }
    idx += 1;
  }

  // Log likelihood for the separate group data
  for (n_neg in 1:N_neg) {
    log_likelihood_all[idx] = normal_lpdf(y_neg[n_neg] | mu[1], sigma[1]);
    total_log_likelihood += log_likelihood_all[idx];

    // Calculate probabilities of each observation belonging to each component
    for (k in 1:K) {
      if (k == 1) {
        prob_y_neg[n_neg, k] = 1;
      } else {
        prob_y_neg[n_neg, k] = 0;
      }
    }
    idx += 1;
  }
}




