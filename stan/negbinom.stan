data {
  int<lower=0> N;                     // number of genes
  int<lower=0> M;                     // number of samples
  int<lower=0> K;                     // number of variants
  matrix[M, K] X;                     // matrix of genotypes
  int<lower=0> Z[N, M];               // array of read counts
  vector<lower=0>[M] c;               // vector of normalisation factors
  vector<lower=0>[M] s;               // vector of library sizes
}

transformed data {
  vector[M] lc = log(c);              // log normalisation factors
  vector[M] ls = log(s);              // log library sizes
}

parameters {
  vector[N] b0;                      // baseline gene expression (log-scale)
  vector[N] lphi;                    // log-dispersion
  real<lower=0> eta;                 // global scale parameter
  vector<lower=0>[K] zeta[N];        // local scale parameters
  vector[K] B[N];                    // regression coefficients
}

transformed parameters {
  vector<lower=0>[N] phi = exp(lphi);         // dispersion
  vector<lower=0>[N] alpha = 1.0 ./ phi;      // inverse-dispersion
}

model {
  real sc = 1.0 / sqrt(N*K);
  eta ~ cauchy(0, 1);
  for(i in 1:N) {
    zeta[i] ~ cauchy(0, 1);
    B[i] ~ normal(0, eta * zeta[i] * sc);
    Z[i] ~ neg_binomial_2_log(ls + lc + b0[i] + X * B[i], alpha[i]);
  }
}
