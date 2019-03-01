library(tidyverse)

load_genotypes = function() {
  # read sample names
  samples =
    file.path('data', 'chr7.012.indv') %>%
    read_tsv(col_names = 'SAMPLE') %>%
    pull(SAMPLE)

  # read loci
  pos =
    file.path('data', 'chr7.012.pos') %>%
    read_tsv(col_names = c('CHROM', 'POS')) %>%
    unite('POS', CHROM, POS, sep = ':') %>%
    pull(POS)

  # read genotypes
  geno =
    file.path('data', 'chr7.012') %>%
    read_tsv(na = '-1', col_names = c('SAMPLE', pos)) %>%
    select(-1) %>%
    as.matrix()
  dimnames(geno) = list(samples = samples, variants = pos)

  # remove loci with the same genotype across all samples
  sds = apply(geno, 2, sd, na.rm = T)
  geno[,sds > 0]
}

load_counts = function(pop = 'CEU') {
  # load samples
  samples =
    read_delim(file.path('data', 'montpick_phenodata.txt'), delim = ' ') %>%
    filter(population == pop) %>%
    pull(sample.id)

  # load read counts
  counts =
    read_tsv(file.path('data', 'montpick_count_table.txt')) %>%
    select(GENE = gene, samples) %>%
    as.data.frame() %>%
    column_to_rownames('GENE') %>%
    as.matrix()
  dimnames(counts) = list(genes = rownames(counts), samples = samples)

  # remove genes with the same number of counts across all samples
  sds = apply(counts, 1, sd, na.rm = T)
  counts[sds > 0,]
}

calculate_norm_factors = function(counts) {
  lcounts = log(counts)
  lgm = rowSums(lcounts) / ncol(lcounts)
  idxs = is.finite(lgm)
  lratios = sweep(lcounts[idxs,], 1, lgm[idxs], '-')
  apply(exp(lratios), 2, median)
}

lognbinom = function(pars, z, cc) {
  mu = pars[1]
  phi = pars[2]

  m = cc * mu
  alpha = 1 / phi
  ll = lgamma(z + alpha) - lgamma(alpha) - lgamma(z + 1) +
    z * log(m) + alpha * log(alpha) - (z + alpha) * log(m + alpha)

  sum(ll)
}

calculate_gene_stats = function(counts) {
  sizes = calculate_norm_factors(counts)
  optim = possibly(optim, otherwise = NULL)
  stats =
    counts %>%
    plyr::alply(1, function(cnts) {
      optim(par = c(1, 1),
            fn = lognbinom,
            z = cnts,
            cc = sizes,
            control = list(fnscale = -1),
            method = 'L-BFGS-B',
            hessian = T,
            lower = 1e-12)
    }, .parallel = T, .dims = T) %>%
    compact() %>%
    discard(~.x$convergence > 0) %>%
    map('par') %>%
    enframe() %>%
    mutate(value = map(value, str_c, collapse = ',')) %>%
    separate(value, into = c('MEAN', 'PHI'), sep = ',', convert = T) %>%
    mutate(VAR = MEAN + PHI * MEAN^2) %>%
    rename(genes = name) %>%
    as.data.frame() %>%
    column_to_rownames('genes') %>%
    as.matrix()
  dimnames(stats) = list(genes = rownames(stats), statistics = colnames(stats))
  stats
}

simulate_data = function(count_stats, genotypes, nsamples = 1000, ngenes = 100, nvars = 50, nhits = 10, rate = 4) {
  # fetch genotypes
  X = genotypes[sample(1:nrow(genotypes), nsamples), sample(1:ncol(genotypes), nvars)]
  X = X[,apply(X, 2, sd, na.rm = T) > 0]
  dimnames(X) = list(samples = str_c('S', 1:nrow(X)), variants = str_c('V', 1:ncol(X)))

  # simulate matrix of coefficients
  B = matrix(0, nrow = ngenes, ncol = ncol(X), dimnames = list(genes = str_c('G', 1:ngenes),
                                                               variants = str_c('V', 1:ncol(X))))
  hits = c(1 + rexp(0.5 * nhits, rate = rate), -1 - rexp(0.5 * nhits, rate = rate))
  B[sample(length(B), length(hits))] = hits


  # simulate counts
  # df = count_stats[sample(1:nrow(count_stats), ngenes),]
  df = count_stats[order(count_stats[,'MEAN'], decreasing = T),][1:ngenes,]

  X0 = scale(X, center = T, scale = T)
  m = sweep(exp(B %*% t(X0)), 1, df[,'MEAN'], '*')
  alpha = 1 / df[,'PHI']
  Z = matrix(rnbinom(length(m), mu = m, size = alpha), nrow = nrow(m))
  dimnames(Z) = list(genes = str_c('G', 1:nrow(Z)), samples = str_c('S', 1:ncol(Z)))

  # output
  lst(B, X, Z, stats = df)
}

fit_model = function(data, model, fcts = calculate_norm_factors(data$Z), ...) {
  Z_tilde = sweep(data$Z, 2, fcts, '/')             # normalise count data
  Y = log(Z_tilde + 1)                              # transform normalised count data
  X0 = scale(data$X, center = T, scale = T)         # standardise genotypes

  # MAP estimation
  fit = rstan::optimizing(object = model,
                          data = list(Z = data$Z,
                                      Y = Y,
                                      X = X0,
                                      c = fcts,
                                      s = colSums(data$Z),
                                      N = nrow(Y),
                                      M = ncol(Y),
                                      K = ncol(X0)),
                          seed = 42,
                          ...)

  # extract estimated matrix of regression coefficients B in vector format
  tibble(EST =
           fit %>%
           pluck('par') %>%
           enframe() %>%
           filter(str_detect(name, '^B\\[')) %>%
           pull(value),
         TRU = as.vector(data$B),
         IDX = 1:length(TRU))
}
