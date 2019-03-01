data {
  int<lower=0> N;                     // number of genes
  int<lower=0> M;                     // number of samples
  int<lower=0> K;                     // number of variants
  matrix[M, K] X;                     // matrix of genotypes
  vector[M] Y[N];                     // matrix of transformed read counts
}

parameters {
  vector[N] b0;                       // baseline gene expression
  vector[N] ls2;                      // log-variance of gene expression
  real<lower=0> eta;                  // global scale parameter
  vector<lower=0>[K] zeta[N];         // local scale parameters
  vector[K] B[N];                     // regression coefficients
}

transformed parameters {
  vector<lower=0>[N] sigma = sqrt(exp(ls2));     // standard deviation of gene expression
}

model {
  real sc = mean(sigma) / sqrt(N*K);
  eta ~ cauchy(0, 1);
  for(i in 1:N) {
    zeta[i] ~ cauchy(0, 1);
    B[i] ~ normal(0, eta * zeta[i] * sc);
    Y[i] ~ normal(b0[i] + X * B[i], sigma[i]);
  }
}
