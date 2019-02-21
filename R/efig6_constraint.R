
methylation_hist = function(save_plot=F) {
  methylation_data = read_delim(gzfile('data/methylation_hist.txt.bgz'), delim = '\t')
  p = ggplot(methylation_data) + aes(x = edge, y = freq) + geom_bar(stat='identity') +
    theme_classic() + 
    geom_vline(xintercept = 0.2, linetype='dashed', color='darkgray') +
    geom_vline(xintercept = 0.6, linetype='dashed', color='darkgray') +
    xlab('Mean methylation') + ylab('Number of bases')
  
  if (save_plot) {
    pdf('e5a_methylation.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

shapes = c('0' = 16, '1' = 15, '2' = 17)
compare_mutation_rates = function(save_plot=F, log=F, legend=T) {
  new_mu_data = read_delim(gzfile('data/mutation_rate_methylation_bins.txt.bgz'), delim = '\t') %>%
    annotate_variant_types()
  old_mu_data = read_delim('data/fordist_1KG_mutation_rate_table.txt', delim = ' ') %>%
    transmute(context = from, ref = substr(from, 2, 2), alt = substr(to, 2, 2), old_mu = mu_snp)
  
  new_mu_data %>%
    left_join(old_mu_data) %>%
    filter(variant_type != 'CpG transition') %$%
    cor.test(mu_snp, old_mu) %>% print
  
  p = new_mu_data %>%
    left_join(old_mu_data) %>%
    ggplot + aes(x = old_mu, y = mu_snp, color = variant_type, shape = as.character(methylation_level)) + 
    geom_point() + scale_shape_manual(values=shapes) +
    theme_classic() + scale_color_manual(values=variant_type_colors, name = '') +
    xlab('Old mu') + ylab('New mu') + guides(color=F, shape=F) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  if (legend) {
    p = p + annotate('text', 0, 1.4e-7, hjust=0, label='CpG transition', color = color_cpg, size = 3) +
      annotate('text', 0, 1.25e-7, hjust=0, label='non-CpG transition', color = color_ti, size = 3) +
      annotate('text', 0, 1.1e-7, hjust=0, label='transversion', color = color_tv, size = 3) +
      annotate('text', 7e-9, 8e-8, hjust=0, label='high', color = color_cpg, size = 3) +
      annotate('text', 7e-9, 6.5e-8, hjust=0, label='medium', color = color_cpg, size = 3) +
      annotate('text', 7e-9, 5e-8, hjust=0, label='low', color = color_cpg, size = 3)+
      annotate('point', 0, 8e-8, color = color_cpg, shape = shapes['2'], size = 2) +
      annotate('point', 0, 6.5e-8, color = color_cpg, shape = shapes['1'], size = 2) +
      annotate('point', 0, 5e-8, color = color_cpg, shape = shapes['0'], size = 2)
  }
  
  if (log) {
    p = p + scale_x_log10() + scale_y_log10()
  }
  if (save_plot) {
    pdf('e5a_mu_comparison.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

load_calibration_data = function(high_coverage_cutoff=40) {
  coverage_data = read_delim(gzfile('data/prop_observed_by_coverage_no_common_pass_filtered_bins.txt.bgz'), delim = '\t')
   
  coverage_data %>%
    filter(exome_coverage >= high_coverage_cutoff) %>%
    group_by(context, ref, alt, methylation_level, mu_snp, variant_type, cpg) %>%
    summarize(obs = sum(variant_count, na.rm=T), poss = sum(possible_variants, na.rm=T), prop_observed = obs / poss) %>%
    group_by(cpg) %>%
    do(tidy(lm(prop_observed ~ mu_snp, data=.))) %>% ungroup -> high_coverage_models
  
  coverage_data %>%
    left_join(high_coverage_models %>% select(cpg, term, estimate)) %>%
    group_by_at(vars(-term, -estimate)) %>%
    summarize(pred_prop_observed = sum((term == 'mu_snp') * mu_snp * estimate + (term == '(Intercept)') * estimate)) %>% 
    ungroup %>% return
}

high_coverage_prop_observed = function(save_plot=F, x_axis='mu_snp', x_label='Mu', legend=F) {
  high_coverage_data = load_calibration_data() %>% 
    filter(exome_coverage >= 40) %>%
    group_by(context, ref, alt, methylation_level, mu_snp, variant_type, cpg, pred_prop_observed) %>%
    summarize(obs = sum(variant_count, na.rm=T), poss = sum(possible_variants, na.rm=T), prop_observed = obs / poss)
  
  p = high_coverage_data %>%
    ggplot + aes_string(x = x_axis) +
    aes(y = prop_observed, color = variant_type, shape = as.character(methylation_level)) + 
    geom_point() + theme_classic() + scale_shape_manual(values=shapes, guide=F) +
    scale_color_manual(values=variant_type_colors) + xlab(x_label) + ylab('Proportion observed')
  
  if (x_axis == 'mu_snp') {
    lms = high_coverage_data %>%
      group_by(cpg) %>%
      do(l = lm(prop_observed ~ mu_snp, data = .)) %>%
      tidy(l) %>%
      select(cpg, term, estimate) %>%
      spread('term', 'estimate')
    p = p + geom_abline(aes(slope = mu_snp, intercept= `(Intercept)`), 
                        data=lms, linetype='dashed', color='darkgray')
  }
  
  if (!legend) {
    p = p + guides(color=FALSE)
  }
  
  if (save_plot) {
    pdf('prop_observed_vs_mu.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

low_coverage_obs_exp = function(save_plot=F) {
  low_coverage_data = load_calibration_data() %>% 
    group_by(exome_coverage) %>%
    summarize(obs = sum(variant_count, na.rm=T), poss = sum(possible_variants, na.rm=T),
              mu_sum=sum(possible_variants * mu_snp, na.rm=T),
              exp=sum(possible_variants * pred_prop_observed), oe = obs / exp)
  
  low_cov_lm = low_coverage_data %>%
    filter(exome_coverage < 40) %>%
    do(l = lm(oe ~ exome_coverage, data = .)) %>% tidy(l) %>%
    select(term, estimate) %>% spread('term', 'estimate')
  
  p = low_coverage_data %>%
    ggplot + aes(x = exome_coverage, y = oe) + geom_point() + theme_classic() +
    xlab('Median coverage') + ylab('Observed / Expected') +
    geom_smooth(aes(exome_coverage, oe),
                data = low_coverage_data %>% filter(exome_coverage < 40 & exome_coverage > 0),
                method = lm, formula = y ~ log10(x), 
                linetype = 'dashed', se=F, color='darkgray'
                )
  
  if (save_plot) {
    pdf('obs_exp_vs_coverage.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

correlation_plots = function(save_plot=F, data_type='syn') {
  plot_gene_data = load_constraint_data() %>%
    filter(gene != 'TTN' & gene != 'MUC16') %>%
    rename(obs=paste0('obs_', data_type),
           exp=paste0('exp_', data_type)) %>%
    filter(!is.na(obs) & !is.na(exp))
  var_type = var_type_aliases[data_type]
  
  r = cor(plot_gene_data$obs, plot_gene_data$exp)
  print(r^2)
  varmax = max(max(plot_gene_data$obs), max(plot_gene_data$exp))
  var_color = colors[[var_type]]

  p = plot_gene_data %>%
    ggplot + aes(x = exp, y = obs) + 
    geom_point_rast(color=var_color, raster.width = 3, raster.height = 9 / 4) +
    theme_classic() + xlab('Expected variants') + ylab('Observed variants') + 
    geom_abline(slope = 1, intercept = 0) +
    annotate('text', x = 0, y = varmax, hjust = 0, color = var_color, label = var_type) +
    annotate('text', x = 0, y = varmax * 0.9, hjust = 0, label = paste('r =', round(r, 4))) +
    xlim(0, varmax) + ylim(0, varmax)
  
  if (save_plot) {
    pdf(paste0('obs_exp_', data_type, '.pdf'), height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

efigure6 = function() {
  e6a = methylation_hist()
  e6b = compare_mutation_rates()
  e6c = high_coverage_prop_observed()
  e6d = high_coverage_prop_observed(x_axis = 'pred_prop_observed', x_label = 'Predicted proportion observed')
  e6e = low_coverage_obs_exp()
  e6f = correlation_plots()
  e6g = correlation_plots(data_type = 'mis')
  e6h = correlation_plots(data_type = 'lof')
  pdf('extended_data_figure6.pdf', height=9, width=6)
  ggarrange(e6a, e6b, e6c, e6d, e6e, e6f, e6g, e6h, nrow = 4, ncol = 2, labels='auto', vjust = 1)
  dev.off()
  png('extended_data_figure6.png', height=9, width=6, units = 'in', res=300)
  ggarrange(e6a, e6b, e6c, e6d, e6e, e6f, e6g, e6h, nrow = 4, ncol = 2, labels='auto', vjust = 1)
  dev.off()
  # print(ef5)
}

