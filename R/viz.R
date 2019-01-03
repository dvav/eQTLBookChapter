library(tidyverse)

plot_genotypes = function(geno) {
  geno %>%
    t() %>%
    cor() %>%
    reshape2::melt() %>%
    ggplot() +
    geom_raster(aes(x = Var1, y = Var2, fill = value)) +
    theme(legend.position = 'none',
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_fill_gradient(low = 'black', high = 'white') +
    labs(x = 'genomic position', y = 'genomic position')
}

plot_mean_variance = function(counts, stats) {
  plot_fn = function(df, xlabel, ylabel) {
    df %>%
      ggplot() +
      geom_point(aes(x = MEAN, y = VAR), size = 0.1) +
      geom_smooth(aes(x = MEAN, y = VAR), size = 0.5,
                  color = 'red', method = 'lm', formula = y ~ poly(x, 2)) +
      scale_x_continuous(trans = 'log10') +
      scale_y_continuous(trans = 'log10') +
      labs(x = xlabel, y = ylabel)
  }

  # normalise count data
  sizes = calculate_norm_factors(counts)
  counts = sweep(counts, 2, sizes, '/')

  # plot of estimated stats
  gg1 =
    stats %>%
    as.data.frame() %>%
    plot_fn(xlabel = 'estimated mean', ylabel = 'estimated variance')

  # plot of observed stats
  gg2 =
    data.frame(
      MEAN = rowMeans(counts),
      VAR = apply(counts, 1, var)
    ) %>%
    plot_fn(xlabel = 'observed mean', ylabel = 'observed variance')

  # combine plots
  cowplot::plot_grid(gg1, gg2, align = 'vh', labels = 'AUTO')
}

plot_fitted_models = function(df, thr = 0.1, ylim = c(-2, 2)) {
  df %>%
    mutate(TRU = na_if(TRU, 0),
           EST = if_else(abs(EST) < thr, NA_real_, EST)) %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', size = 0.2) +
    geom_point(aes(x = IDX, y = TRU), color = 'red') +
    geom_point(aes(x = IDX, y = EST), shape = 4) +
    facet_grid(MODEL~NSAMPLES) +
    scale_y_continuous(limits = ylim) +
    labs(x = 'index of regression coefficients (genes x variants)',
         y = 'value of regression coefficients')
}

