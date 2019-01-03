data {
  int<lower=0> N;                     // number of genes
  int<lower=0> M;                     // number of samples
  int<lower=0> K;                     // number of variants
  int Z[N, M];                        // matrix of read counts
  matrix[K, M] X;                     // matrix of genotypes
  row_vector<lower=0>[M] s;           // vector of library sizes
}

transformed data {
  row_vector[M] ls = log(s); // log library sizes
}

parameters {
  vector[N] b0;                      // baseline gene expression (log-scale)
  matrix[N, K] B;                    // matrix of regression coefficients
  vector[N] lphi;                    // log-dispersion
  real<lower=0> eta;                 // global scale parameter
  matrix<lower=0>[N, K] zeta;        // local scale parameter
}

transformed parameters {
  vector<lower=0>[N] phi = exp(lphi);        // dispersion
  vector<lower=0>[N] alpha = 1.0 ./ phi;     // inverse-dispersion
}

model {
  matrix[N, M] coefs = B * X;
  eta ~ cauchy(0, 1);
  for(i in 1:N) {
    zeta[i] ~ cauchy(0, 1);
    B[i] ~ normal(0, eta * zeta[i]);
    Z[i] ~ neg_binomial_2_log(ls + b0[i] + coefs[i], alpha[i]);
  }
}
