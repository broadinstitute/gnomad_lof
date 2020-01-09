sum_stats = load_sum_stats()
tau_star = function(save_plot=F, normalize=F) {
  
  filt_sum_stats = sum_stats %>%
    filter(cond == 'all100' &
             grepl('oe_lof_upper_quantile_', name) &
             (is.na(cor_rm0.2) | cor_rm0.2 == 0)  # Select uncorrelated variables
    )
  tau_data = ddply(filt_sum_stats, 'name', function(x) {
    taustar = x$taustar
    se = x$se
    meta_enrichment = metagen(taustar, se)$TE.random
    meta_sd = metagen(taustar, se)$seTE.random
    return(data.frame(meta_enrichment, meta_sd))
  }) %>%
    mutate(decile = as.integer(str_sub(name, -1)),
           enrichment = ifelse(normalize, meta_enrichment / mean(meta_enrichment), meta_enrichment)
    )
  
  p = tau_data %>%
    ggplot + aes(x = decile, y = meta_enrichment,
                 ymin = meta_enrichment - 1.96 * meta_sd, ymax = meta_enrichment + 1.96 * meta_sd) + 
    geom_pointrange() +
    theme_classic() + ylab(expression(paste(tau^"*"))) + oe_x_axis +
    geom_hline(yintercept = 0, col = "darkgray", linetype = 'dashed')
  
  if (save_plot) {
    pdf('e9a_tau_star.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

enrichment_corrections = function(save_plot=F) {

  filt_sum_stats = sum_stats %>%
    filter(cond %in% c('all100', 'condbrain', 'condgenesize') &
             grepl('oe_lof_upper_quantile_', name) &
             (is.na(cor_rm0.2) | cor_rm0.2 == 0)  # Select uncorrelated variables
    )
  enrichment_correction_data = ddply(filt_sum_stats, c('name', 'cond'), function(x) {
    taustar = x$taustar
    se = x$se
    meta_enrichment = metagen(taustar, se)$TE.random
    meta_sd = metagen(taustar, se)$seTE.random
    return(data.frame(meta_enrichment, meta_sd))
  }) %>%
    mutate(decile = as.integer(str_sub(name, -1)),
           cond=fct_recode(cond, 'Baseline' = 'all100', 
                           'Baseline + brain\nexpression' = 'condbrain', 
                           'Baseline + gene size +\nexon count' = 'condgenesize')
    )
  
  p = enrichment_correction_data %>%
    ggplot + aes(x = decile, y = meta_enrichment, 
                 ymin = meta_enrichment - 1.96 * meta_sd, ymax = meta_enrichment + 1.96 * meta_sd) + 
    geom_point(stat="identity", col="black", position=position_dodge()) + 
    geom_pointrange() + 
    theme_classic() + ylab(expression(paste(tau^"*"))) + 
    geom_hline(yintercept = 0, col = "darkgray", linetype = 'dashed') +
    oe_x_axis + facet_wrap(~cond) #+ theme(panel.spacing = unit(0.1, "lines"))
  
  if (save_plot) {
    pdf('e10b_enrichment_corrections.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

window_assessment = function(save_plot=F) {
  filt_data = sum_stats %>%
    filter(cond %in% c('all100', 'all50', 'all10') &
             grepl('oe_lof_upper_quantile_', name) &
             (is.na(cor_rm0.2) | cor_rm0.2 == 0)  # Select uncorrelated variables
    )
  
  window_data = ddply(filt_data, c('name', 'cond'), function(x) {
    enrichment = x$enrichment
    enrichment_SE = x$enrichment_SE
    meta_enrichment = metagen(enrichment, enrichment_SE)$TE.random
    meta_sd = metagen(enrichment, enrichment_SE)$seTE.random
    return(data.frame(meta_enrichment, meta_sd))
  }) %>%
    mutate(decile = as.integer(str_sub(name, -1)),
           cond=fct_recode(cond, '100 kb' = 'all100',
                           '50 kb' = 'all50', 
                           '10 kb' = 'all10'),
           cond=fct_relevel(cond, '100 kb', '50 kb', '10 kb')
    )
  
  p = window_data %>%
    ggplot + aes(x = decile, y = meta_enrichment, col = cond,
                 ymin = meta_enrichment - 1.96 * meta_sd, ymax = meta_enrichment + 1.96 * meta_sd) + 
    geom_pointrange(position=position_dodge(width=0.4)) +
    theme_classic() + ylab('Enrichment') + scale_color_discrete(name="Window size") + 
    oe_x_axis +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1))
  
  if (save_plot) {
    pdf('e9c_window_assessment.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

efigure10 = function() {
  e10a = tau_star()
  e10b = window_assessment()
  e10c = enrichment_corrections()
  e10top = ggarrange(e10a, e10b, nrow = 1, ncol = 2, labels = 'auto', vjust = 1)
  pdf('extended_data_figure10.pdf', height=6, width=7.5)
  print(ggarrange(e10top, e10c, nrow = 2, ncol = 1, labels = c('', 'c'), vjust = 1))
  dev.off()
  png('extended_data_figure10.png', height=6, width=7.5, units = 'in', res=300)
  print(ggarrange(e10top, e10c, nrow = 2, ncol = 1, labels = c('', 'c'), vjust = 1))
  dev.off()
}


