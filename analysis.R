library(tidyverse)
source('R/utils.R')
source('R/viz.R')

doMC::registerDoMC(8)

obs_geno = load_genotypes()
plot_genotypes(obs_geno)

obs_counts = load_counts()
obs_counts_stats = calculate_gene_stats(obs_counts)
plot_mean_variance(obs_counts, obs_counts_stats)

normal = rstan::stan_model(file.path('stan', 'normal.stan'), auto_write = T)
poissonln = rstan::stan_model(file.path('stan', 'poissonln.stan'), auto_write = T)
negbinom = rstan::stan_model(file.path('stan', 'negbinom.stan'), auto_write = T)

cpp_object_initializer = rstan::cpp_object_initializer

set.seed(42)

fitted_models =
  plyr::ldply(lst(500, 1500, 2500), function(N) {
    sim = simulate_data(obs_counts_stats, obs_geno, nsamples = N)
    plyr::ldply(lst(normal, poissonln, negbinom), function(mdl) {
      fit_model(sim, mdl)
    }, .parallel = T, .id = 'MODEL')
  }, .progress = 'text', .id = 'NSAMPLES') %>%
  as_data_frame()

plot_fitted_models(fitted_models)
