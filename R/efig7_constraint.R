gene_data = load_constraint_data()

oe_2d_density = function(save_plot=F, recompute_ridges=F) {
  upper_threshold = 100
  breaks = c(0, 1, 3, 10, 30, 100)
  
  if (recompute_ridges) {
    # all_combos = data.frame(obs = rep(0:100, times = 1001), exp = rep(0:1000, each = 101)/10)
    all_combos = data.frame(obs = rep(0:100, times = 3001), exp = rep(0:3000, each = 101)/100) %>%
      union(data.frame(obs = rep(0:100, times = 701), exp = rep(300:1000, each = 101)/10))
    
    # Confirm coverage:
    # all_combos %>% ggplot + aes(x = exp, y = obs) + geom_bin2d(bins=10)
    
    all_combos %>%
      rowwise %>%
      do(cbind(., generate_ci(.$obs, .$exp))) %>%
      ungroup -> all_cis
  } else {
    fname = get_or_download_file('all_cis.RData', subfolder = 'misc_files/')
    all_cis = get(load(fname))
  }
  
  gene_data %>%
    filter(!is.na(oe_lof_upper) & oe_lof_upper_bin < 9) %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(max_limit=max(oe_lof_upper),
              min_limit=min(oe_lof_upper)) -> limits
  
  all_cis %>%
    crossing(limits) %>%
    filter(abs(upper - max_limit) <= 0.005) -> boundaries
  
  p = gene_data %>%
    mutate(exp_lof_int = ifelse(exp_lof > upper_threshold, upper_threshold, floor(exp_lof)),
           obs_lof_int = ifelse(obs_lof > upper_threshold, upper_threshold, obs_lof)) %>%
    # count(obs_lof_int, exp_lof_int) %>%
    ggplot + aes(x = exp_lof_int + 1, y = obs_lof_int + 1) + theme_classic() + 
    # geom_bin2d(aes(alpha=..count..), bins = 12, fill='darkblue') +
    geom_bin2d(bins = 12) +
    scale_x_log10(breaks = breaks + 1, labels = breaks) + 
    scale_y_log10(breaks = breaks + 1, labels = breaks) +
    xlab('pLoFs expected') + ylab('pLoFs observed') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') + 
    geom_line(aes(x = exp + 1, y = obs + 1, color = oe_lof_upper_bin + 1, group = oe_lof_upper_bin + 1), 
              data = boundaries) +
    scale_color_gradientn(colors = gradient_colors, values = rescale(gradient_values), 
                          name=paste(constraint_metric_name, 'decile', 'border', sep='\n'),
                          breaks=seq(1, 9, 2)) +
    scale_fill_gradient(name='\nNumber of\ngenes', low='white', high='springgreen4') +
    # scale_alpha_continuous(name='Number of\ngenes') +
    theme(legend.title = element_text(size = 9), legend.box = "horizontal")
  
  if (save_plot) {
    pdf('oe_2d_density.pdf', height=3, width=5)
    print(p)
    dev.off()
  }
  return(p)
}

oe_distribution = function(save_plot=F, metric='oe_lof') {
  # gene_data %>%
  #   filter(oe_lof < 1.5 & exp_lof > 0) %$%
  #   density(oe_lof, n = 2^12) %$%
  #   data.frame(x=x, y=y) -> data_dens
  # 
  # p = data_dens %>%
  #   ggplot + aes(x = x, y = y) + geom_line() + 
  #   geom_segment(aes(xend = x, yend = 0, color = x)) +
  #   scale_color_gradientn(colors = gradient_colors, values = rescale(gradient_values)) +
  #   xlim(c(0, 1.5)) + theme_classic() + xlab('Observed/Expected') + ylab('Density') + guides(color=F)
  
  gene_data %>%
    summarize(mean = mean(oe_lof, na.rm=T),
              median = median(oe_lof, na.rm=T),
              zeroes = sum(oe_lof == 0, na.rm=T),
              zeroes_few_expected = sum(oe_lof == 0 & exp_lof <= 5, na.rm=T),
              mean_upper = mean(oe_lof_upper, na.rm=T),
              median_upper = median(oe_lof_upper, na.rm=T))

  n_bins = 30
  plot_data = gene_data %>%
    filter(oe_lof < 1.5 & exp_lof > 0) %>%
    mutate_at(metric, funs(plot_metric = . + 0)) %>%
    mutate(plot_metric_bin = floor((plot_metric - min(plot_metric)) * n_bins / max(plot_metric))) %>%
    group_by(plot_metric_bin) %>%
    summarize(mid=min(plot_metric), count=n()) 
  
  xlabel = case_when(metric == 'oe_lof' ~ 'Observed/Expected',
                     metric == 'oe_lof_upper' ~ constraint_metric_name,
                     TRUE ~ metric)
  
  p = plot_data %>%
    ggplot + aes(x = mid, y = count, fill = mid, color = mid) + 
    geom_bar(stat='identity', width=max(plot_data$mid)/n_bins) + 
    scale_color_gradientn(colors = gradient_colors, values = rescale(gradient_values)) +
    scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values)) +
    theme_classic() + xlab(xlabel) + ylab('Number of genes') + guides(fill=F, color=F)
  
  if (save_plot) {
    pdf('oe_distribution.pdf', height=2, width=3)
    print(p)
    dev.off()
  }
  return(p)
}

proportion_high_pLI = function(save_plot=F) {
  p = gene_data %>%
    filter(!is.na(oe_lof_upper_bin)) %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(num_high_pLI = sum(exac_pLI > 0.9, na.rm=T), n=n(), 
              prop_high_pLI = num_high_pLI / n, sem = sqrt(prop_high_pLI * (1 - prop_high_pLI) / n),
              prop_upper = prop_high_pLI + 1.96 * sem, prop_lower = prop_high_pLI - 1.96 * sem) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_high_pLI, ymin = prop_lower, ymax = prop_upper) + 
    geom_pointrange() + scale_y_continuous(labels=percent_format(accuracy = 1)) +
    theme_classic() + oe_x_axis + ylab('Percent high\npLI (ExAC)')
  
  if (save_plot) {
    pdf('proportion_high_pLI.pdf', height=3, width=5)
    print(p)
    dev.off()
  }
  return(p)
}

caf_vs_constraint = function(save_plot=F) {
  gene_data %$%
    cor.test(oe_lof_upper, p, method='spearman')
  gene_data %>%
    filter(p > 0) %$%
    cor.test(oe_lof_upper, log10(p))
  
  p = gene_data %>%
    filter(exp_lof > 10) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = p, group = oe_lof_upper_bin, fill = oe_lof_upper_bin) + 
    geom_boxplot() + scale_y_log10() + theme_classic() + oe_x_axis +
    scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values), guide = FALSE) +
    ylab('pLoF frequency')
  
  if (save_plot) {
    pdf('caf_vs_constraint.pdf', height=3, width=5)
    print(p)
    dev.off()
  }
  return(p)
}

fname = get_or_download_file('constraint_standard.txt.bgz', subfolder = 'standard/')
full_constraint_data = read_delim(gzfile(fname), delim = '\t') %>%
  mutate(obs_os = obs_lof_with_os - obs_lof,
         exp_os = exp_lof_with_os - exp_lof)

variant_depletions = function(save_plot=F, legend=T) {
  # full_constraint_data %>% ggplot + aes(x = exp_os) + geom_density() + scale_x_log10()
  # p = full_constraint_data %>%
  #   filter(canonical == 'true') %>%
  #   mutate(oe_lof_upper_bin=ntile(oe_lof_upper, 10)) %>%
  #   group_by(oe_lof_upper_bin) %>%
  #   summarize(obs_os=sum(obs_os, na.rm=T), exp_os=sum(exp_os, na.rm=T),
  #             oe_os=obs_os/exp_os, n=n()) %>%
  #   ggplot + aes(x = oe_lof_upper_bin, y = oe_os) + 
  #   geom_point(color=color_lof, size=2) + theme_classic() +
  #   oe_x_axis + geom_hline(yintercept = 1, linetype='dashed', color='darkgray') +
  #   scale_y_continuous(label=percent_format(accuracy = 1)) +
  #   ylab('Observed/Expected\nOther splice variants')
  
  p = full_constraint_data %>%
    filter(canonical & exp_lof > 10) %>%
    mutate(oe_lof_upper_bin=ntile(oe_lof_upper, 10) - 1) %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(obs_os=sum(obs_os, na.rm=T), 
              exp_os=sum(exp_os, na.rm=T),
              oe_os=obs_os/exp_os,
              obs_lof=sum(obs_lof, na.rm=T), 
              exp_lof=sum(exp_lof, na.rm=T),
              oe_lof=obs_lof/exp_lof, n=n(),
              obs_mis=sum(obs_mis, na.rm=T), 
              exp_mis=sum(exp_mis, na.rm=T),
              oe_mis=obs_mis/exp_mis,
              obs_syn=sum(obs_syn, na.rm=T), 
              exp_syn=sum(exp_syn, na.rm=T),
              oe_syn=obs_syn/exp_syn, n=n()) %>%
    gather('variant_type', 'value', oe_lof, oe_os, oe_mis, oe_syn) %>%
    separate('variant_type', c('oe', 'variant_type')) %>%
    mutate(variant_type=fct_relevel(var_type_aliases[variant_type],
                                    'Synonymous', 'Missense')) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = value, color=variant_type) + 
    geom_point(size=2) + theme_classic() +
    oe_x_axis + geom_hline(yintercept = 1, linetype='dashed', color='darkgray') +
    scale_color_manual(values=colors, name=NULL) + 
    scale_y_continuous(label=percent_format(accuracy = 1),
                       breaks=seq(2, 12, 2)/10) +
    ylab('Observed/Expected') + guides(color=F)
  
  if (legend) {
    # p = p + theme(legend.justification = c(1, 0), legend.position = c(1, 0.02),
    #               legend.background = element_rect(fill=alpha('white', 0)))
    top_legend = 0.5
    p = p + annotate('text', 9, top_legend, hjust=1, vjust=1, label='synonymous', color = color_syn, size = 3) +
      annotate('text', 9, 0.8*top_legend, hjust=1, vjust=1, label='missense', color = color_mis, size = 3) +
      annotate('text', 9, 0.6*top_legend, hjust=1, vjust=1, label='other splice', color = color_os, size = 3) +
      annotate('text', 9, 0.4*top_legend, hjust=1, vjust=1, label='pLoF', color = color_lof, size = 3)
  }
  
  if (save_plot) {
    pdf('e5_os_depletion.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

prop_homozygous = function(save_plot=F) {
  (p = gene_data %>%
     filter(exp_lof >= 10) %>%
     mutate(obs_hom_lof = ifelse(obs_hom_lof == 0, NA, obs_hom_lof)) %>%
     group_by(oe_lof_upper_bin) %>%
     summarize(at_least_one_hom = sum(obs_hom_lof > 0, na.rm=T),
               total = n(), prop_hom = at_least_one_hom / total,
               sem = 1.96*sqrt(prop_hom * (1 - prop_hom) / total)) %>% 
               {print(.)} %>%
     ggplot + aes(x = oe_lof_upper_bin, y = prop_hom, color = oe_lof_upper_bin) + geom_point() +
     geom_pointrange(aes(ymin = prop_hom - sem, ymax = prop_hom + sem)) + 
     theme_classic() + oe_x_axis +
     scale_y_continuous(labels=percent_format(accuracy = 1), limits = c(0, NA)) +
     scale_color_gradientn(colors = gradient_colors, values = rescale(gradient_value_deciles), guide='none') +
     ylab('Percent with at least\n1 homozygous pLoF'))
  if (save_plot) {
    pdf('proportion_hom.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

efigure7 = function() {
  e7a = oe_distribution()
  e7b = oe_distribution(metric = 'oe_lof_upper')
  e7c = oe_2d_density()
  e7d = variant_depletions()
  e7e = proportion_high_pLI()
  e7f = prop_homozygous()
  e7g = caf_vs_constraint()
  e7row1 = egg::ggarrange(plots=list(e7a, e7b), labels=c('a', 'b'), align = 'v', ncol = 2, 
                        label.args=list(gp=grid::gpar(font=2, cex=1.2)))
  e7row34 = egg::ggarrange(plots=list(e7d, e7e, e7f, e7g), labels=c('d', 'e', 'f', 'g'), align = 'v', nrow = 2, ncol = 2, 
                         label.args=list(gp=grid::gpar(font=2, cex=1.2)))
  pdf('extended_data_figure7.pdf', height=8, width=7.5, onefile = F)
  # egg::ggarrange(plots=list(a, b, c, d, e, f, g),
  #                nrow = 4, ncol = 2, align = 'v', 
  #                labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g'), 
  #                label.args=list(gp=grid::gpar(font=2, cex=1.2)), vjust = 1)
  ggarrange(e7row1, egg::set_panel_size(e7c, width=unit(2, 'in')), e7row34, align='hv',
            nrow = 3, ncol = 1,
            labels = c('', 'c', ''),
            heights = c(1, 1, 2)
  )
  dev.off()
  png('extended_data_figure7.png', height=8, width=7.5, units = 'in', res=300)
  # egg::ggarrange(plots=list(a, b, c, d, e, f, g),
  #                nrow = 4, ncol = 2, align = 'v', 
  #                labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g'), 
  #                label.args=list(gp=grid::gpar(font=2, cex=1.2)), vjust = 1)
  ggarrange(e7row1, egg::set_panel_size(e7c, width=unit(2, 'in')), e7row34, align='hv',
            nrow = 3, ncol = 1,
            labels = c('', 'c', ''),
            heights = c(1, 1, 2)
  )
  dev.off()
 }




