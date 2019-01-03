data {
  int<lower=0> N;                     // number of genes
  int<lower=0> M;                     // number of samples
  int<lower=0> K;                     // number of variants
  matrix[N, M] Y;                     // matrix of transformed read counts
  matrix[K, M] X;                     // matrix of genotypes
}

parameters {
  vector[N] b0;                       // baseline gene expression
  matrix[N, K] B;                     // matrix of regression coefficients
  vector[N] ls2;                      // log-variance of gene expression
  real<lower=0> eta;                  // global scale parameter
  matrix<lower=0>[N, K] zeta;         // local scale parameter
}

transformed parameters {
  vector<lower=0>[N] sigma = sqrt(exp(ls2));     // standard deviation of gene expression
}

model {
  matrix[N, M] coefs = B * X;
  real sc = mean(sigma) / sqrt(N*K);
  eta ~ cauchy(0, 1);
  for(i in 1:N) {
    zeta[i] ~ cauchy(0, 1);
    B[i] ~ normal(0, eta * zeta[i] * sc);
    Y[i] ~ normal(b0[i] + coefs[i], sigma[i]);
  }
}
