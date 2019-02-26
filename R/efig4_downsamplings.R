data_type = 'exomes'
obs_poss_ds = load_observed_possible_data(data_type)
indel_summary = load_indel_data(data_type)

num_observed_by_class = function(save_plot=F) {
  pe4a = obs_poss_ds %>%
    filter(csq == 'synonymous' & downsampling >= 100) %>%
    group_by(variant_type, downsampling) %>%
    summarize(observed = sum(observed), possible = sum(as.numeric(possible))) %>%
    ggplot + aes(x = downsampling, y = observed, color = variant_type) %>%
    geom_line(lwd=2) + theme_classic() + 
    downsampling_x_axis() + scale_y_log10(labels=comma) +
    scale_color_manual(values=variant_type_colors, name='') + 
    xlab('Sample size') + ylab('Number of variants observed') +
    guides(color=F) + 
    annotate('text', 100, 1e6, hjust=0, label='CpG transition', color = color_cpg, size = 3) +
    annotate('text', 100, 0.7e6, hjust=0, label='non-CpG transition', color = color_ti, size = 3) +
    annotate('text', 100, 0.49e6, hjust=0, label='transversion', color = color_tv, size = 3)

  if (save_plot) {
    pdf('e4a_num_observed_by_class.pdf', height=3, width=4)
    print(pe4a)
    dev.off()
  }
  return(pe4a)
}

prop_observed_by_class = function(save_plot=F, plot_log=T) {
  pe4b = obs_poss_ds %>%
    filter(csq == 'synonymous' & downsampling >= 100 & coverage >= 30) %>%
    group_by(variant_type, downsampling) %>%
    summarize(observed = sum(observed), possible = sum(as.numeric(possible))) %>%
    ggplot + aes(x = downsampling, y = observed/possible, color = variant_type) %>%
    geom_line(lwd=2) + theme_classic() + 
    scale_color_manual(values=variant_type_colors, name='') + 
    downsampling_x_axis(log=plot_log) +
    ylab('Percent observed') + guides(color=F)
  
  if (plot_log) pe4b = pe4b + scale_y_log10(labels=percent)
  
  if (save_plot) {
    pdf('e4b_prop_observed_by_class.pdf', height=3, width=4)
    print(pe4b)
    dev.off()
  }
  return(pe4b)
}

titv = function(save_plot=F, plot_log=T, syn=T) {
  pe4c = obs_poss_ds %>%
    filter((syn | csq == 'synonymous') & downsampling >= 100) %>%
    group_by(variant_type, downsampling) %>%
    summarize(observed = sum(observed), possible = sum(as.numeric(possible))) %>%
    group_by(downsampling) %>%
    summarize(ti = sum(observed * (variant_type != 'transversion')),
              tv = sum(observed * (variant_type == 'transversion')),
              titv = ti / tv) %>%
    ggplot + aes(x = downsampling, y = titv) %>%
    geom_line(lwd=2) + theme_classic() +
    ylab('Ti/Tv ratio') + xlab('Sample size')
  
  if (plot_log) {
    pe4c = pe4c + scale_x_log10(labels=comma, breaks=ds_breaks_log)
  }
  
  if (save_plot) {
    pdf('e4c_titv.pdf', height=3, width=4)
    print(pe4c)
    dev.off()
  }
  return(pe4c)
}

obs_poss_supp_table = function() {
  data_type = 'exomes'
  obs_poss_ds = load_observed_possible_data(data_type)
  indel_summary = load_indel_data(data_type)
  
  obs_poss_ds %>%
    filter(csq %in% c('synonymous', 'nonsense', 'missense') & downsampling == 125748 & coverage >= 30) %>%
    group_by(csq, variant_type, methylation_level) %>%
    summarize(observed = sum(observed), possible = sum(as.numeric(possible)), prop = observed / possible, 
              singletons = sum(singletons), prop_singleton = singletons / observed) %>% ungroup -> snp_collapsed
  
  indel_summary %>%
    filter(csq %in% c('frameshift', 'inframe insertion', 'inframe deletion') & coverage >= 30) %>%
    group_by(csq, variant_type = 'indel') %>%
    summarize(observed = sum(observed), singletons=sum(singletons), 
              prop_singleton = singletons / observed) %>% ungroup -> indel_collapsed
  
  snp_collapsed %>%
    bind_rows(indel_collapsed) %>%
    mutate(prop = sprintf('%.2f', 100 * prop),
           prop_singleton = sprintf('%.2f', 100 * prop_singleton)) %>%
    write_delim('data/s_table_number_observed.tsv', delim='\t')
}

obs_poss_ds %>%
  filter(csq %in% c('synonymous', 'nonsense', 'missense') & downsampling == 125748 & coverage >= 30) %>%
  group_by(csq) %>%
  summarize(observed = sum(observed), possible = sum(as.numeric(possible)), 
            prop = observed / possible)

observed_by_function = function(save_plot=F, plot_log=T, plot_type='prop_observed', grouping=c('csq', 'downsampling'), plot_csqs=c('synonymous', 'nonsense', 'missense')) {
  plotting_csq = 'csq' %in% grouping
  plot_data = obs_poss_ds
  if ('pLoF' %in% plot_csqs) {
    plot_data = plot_data %>%
      mutate(csq = ifelse(csq %in% c('splice donor', 'splice acceptor', 'nonsense'), 'pLoF', csq))
  }
  plot_data = plot_data %>%
    filter((!plotting_csq | csq %in% plot_csqs) & downsampling >= 100 & coverage >= 30) %>%
    group_by_at(vars(grouping)) %>%
    summarize(observed = sum(observed), possible = sum(as.numeric(possible)),
              prop_observed = observed / possible)
  ytop = max(plot_data[,plot_type]) * ifelse(save_plot, 0.98, 0.9)
  y1 = ytop
  pe4d = plot_data %>%
    ggplot + aes(x = downsampling)
  if (plotting_csq) {
    pe4d = pe4d + aes(color = csq)
    for (csq in plot_csqs) {
      pe4d = pe4d + 
        annotate('text', 100, y1, hjust=0, label=csq, color = colors[[csq]], size = 3)
      y1 = ifelse(plot_log, 0.7 * y1, y1 - ytop / 10)
    }
  }
  pe4d = pe4d + aes_string(y = plot_type) +
    geom_line(lwd=ifelse(save_plot, 1.5, 2)) + theme_classic() +
    scale_color_manual(values=c(variant_category_colors, colors), name='') + 
    guides(color=F) +
    xlab('Sample size') + 
    ylab(if_else(plot_type=='prop_observed', 'Percent observed', 'Number observed')) 
  axis_format = ifelse(plot_type == 'prop_observed', percent, comma)
  if (data_type == 'genomes') {
    ds_breaks = c(0, 4e3, 8e3, 12e3, 15708)
    ds_breaks_log = c(1e2, 1e3, 15708)
  } else {
    ds_breaks = c(0, 4e4, 8e4, 125748)
    ds_breaks_log = c(1e2, 1e3, 1e4, 125748)
  }
  if (plot_log) {
    pe4d = pe4d + scale_x_log10(labels=comma, breaks=ds_breaks_log) + 
      scale_y_log10(labels=axis_format)
  } else {
    pe4d = pe4d + scale_x_continuous(labels=comma, breaks=ds_breaks) + 
      scale_y_continuous(labels=axis_format)#, breaks=seq(0, 200e6, 20e6))
  }
  
  if (save_plot) {
    pdf(paste0('e4b_', plot_type, '_by_', 
               paste(grouping, collapse = '_'), '.pdf'), height=3, width=4)
    print(pe4d)
    dev.off()
  }
  return(pe4d)
}

efigure4 = function() {
  e4a = num_observed_by_class()
  e4b = prop_observed_by_class()
  e4c = titv()
  e4d = observed_by_function()
  
  pdf('extended_data_figure4.pdf', height=6, width=8)
  ggarrange(e4a, e4b, e4c, e4d, ncol = 2, nrow = 2, labels = 'auto', align = 'v')
  dev.off()
  png('extended_data_figure4.png', height=6, width=8, units = 'in', res=300)
  ggarrange(e4a, e4b, e4c, e4d, ncol = 2, nrow = 2, labels = 'auto', align = 'v')
  dev.off()
}
  


