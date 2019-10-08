source('fig1_summary.R')
data_type = 'exomes'
obs_poss_ds = load_observed_possible_data(data_type)
indel_summary = load_indel_data(data_type)
indel_ds_file = get_or_download_file('indels_summary_exomes.downsampling.txt.bgz', subfolder = 'summary_results/')
indel_ds = read_delim(gzfile(indel_ds_file), delim='\t')

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

prop_observed_by_class = function(save_plot=F, plot_log=T, split_methylation=F) {
  regroupings = c('variant_type', 'downsampling')
  if (split_methylation) {
    regroupings = c(regroupings, 'methylation_level')
  }
  pe4b = obs_poss_ds %>%
    mutate(methylation_level = as.factor(methylation_level)) %>%
    filter(csq == 'synonymous' & downsampling >= 100 & coverage >= 30) %>%
    group_by_at(regroupings) %>%
    summarize(observed = sum(observed), possible = sum(as.numeric(possible))) %>%
    ggplot + aes(x = downsampling, y = observed/possible, color = variant_type) +
    geom_line(lwd=2) + theme_classic() + 
    scale_color_manual(values=variant_type_colors, name='') + 
    downsampling_x_axis(log=plot_log) +
    ylab('Percent observed') + guides(color=F)
  
  if (split_methylation) {
    pe4b = pe4b + aes(linetype = methylation_level) + 
      scale_linetype_manual(name='Methylation\nLevel', limits=c(2, 1, 0),
                            labels=c('High', 'Medium', 'Unmethylated'),
                            values=c('dotted', '22', 'solid'))
  }

  if (plot_log) pe4b = pe4b + scale_y_log10(labels=percent)
  
  if (save_plot) {
    pdf(paste0('e4b_prop_observed_by_class', ifelse(split_methylation, '_methylation_split', ''), '.pdf'), height=3, width=4)
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

observed_by_function = function(save_plot=F, plot_log=T, plot_type='prop_observed', grouping=c('csq', 'downsampling'), plot_csqs=c('synonymous', 'nonsense', 'missense'), coverage_cutoff = 30) {
  plotting_csq = 'csq' %in% grouping
  plot_data = obs_poss_ds
  if ('pLoF' %in% plot_csqs) {
    plot_data = plot_data %>%
      mutate(csq = ifelse(csq %in% c('splice donor', 'splice acceptor', 'nonsense'), 'pLoF', csq))
    plot_data = indel_ds %>%
      filter(worst_csq == 'frameshift_variant' & downsampling >= 100 & coverage >= coverage_cutoff) %>%
      mutate(csq = 'pLoF') %>%
      union_all(plot_data)
  }
  plot_data = plot_data %>%
    filter((!plotting_csq | csq %in% plot_csqs) & downsampling >= 100 & coverage >= coverage_cutoff) %>%
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
  e4e = downsampling_by_pop()
  e4f = downsampling_by_pop(plot_log=F)
  extra_margin = 6
  pdf('extended_data_figure4.pdf', height=9, width=8.5)
  print(ggarrange(e4a + theme(plot.margin = unit(c(5.5, 5.5, 5.5 + extra_margin, 5.5), "pt")), 
                  e4b + theme(plot.margin = unit(c(5.5, 5.5, 5.5 + extra_margin, 5.5), "pt")),
                  e4c + theme(plot.margin = unit(c(5.5 + extra_margin, 5.5, 5.5 + extra_margin, 5.5), "pt")),
                  e4d + theme(plot.margin = unit(c(5.5 + extra_margin, 5.5, 5.5 + extra_margin, 5.5), "pt")),
                  e4e + theme(plot.margin = unit(c(5.5 + extra_margin, 5.5, 5.5, 5.5), "pt")),
                  e4f + theme(plot.margin = unit(c(5.5 + extra_margin, 5.5, 5.5, 5.5), "pt")),
                  ncol = 2, nrow = 3, labels = 'auto', align = 'v'))
  dev.off()
  png('extended_data_figure4.png', height=9, width=8.5, units = 'in', res=300)
  print(ggarrange(e4a, e4b, e4c, e4d, e4e, e4f, ncol = 2, nrow = 3, labels = 'auto', align = 'v'))
  dev.off()
}

sfigure5 = function() {
  fig1split = summary_figure('exomes', group_splice = T, intergenic = F,
                             maps_limits=c(NA, 0.18), po_limits=c(NA, 1), split_methylation = T, keep_x_labels = T)
  s5a = fig1split[[2]]
  s5b = prop_observed_by_class(split_methylation=T)
  pdf('supplementary_figure5.pdf', height=7, width=8)
  print(ggarrange(s5a, s5b, ncol = 1, nrow = 2, labels = 'auto', align = 'v'))
  dev.off()
  png('supplementary_figure5.png', height=7, width=8, units = 'in', res=300)
  print(ggarrange(s5a, s5b, ncol = 1, nrow = 2, labels = 'auto', align = 'v'))
  dev.off()
}

ds_data = load_downsampled_gene_data()

downsampling_by_pop = function(save_plot=F, plot_log=T) {
  collapsed_ds = ds_data %>% filter(canonical) %>%
    group_by(downsampling, pop) %>%
    summarize_if(is.numeric, sum, na.rm=T) %>% ungroup
  
  collapsed_ds = ds_data %>% filter(canonical) %>%
    group_by(downsampling, pop) %>%
    summarize_at(vars(exp_syn:caf), .funs = funs(sum = sum(., na.rm=T),
                                                 over5raw = sum(. >= 5, na.rm=T)/n(),
                                                 over10raw = sum(. >= 10, na.rm=T)/n(),
                                                 over20raw = sum(. >= 20, na.rm=T)/n(),
                                                 over5 = sum(. >= 5, na.rm=T)/sum(. > 0, na.rm=T),
                                                 over10 = sum(. >= 10, na.rm=T)/sum(. > 0, na.rm=T),
                                                 over20 = sum(. >= 20, na.rm=T)/sum(. > 0, na.rm=T))
    ) %>% ungroup
  
  all_pops = collapsed_ds %>%
    gather(key='metric', value='count', -downsampling, -pop) %>%
    separate(metric, into = c('obs_exp', 'variant_type', 'func')) %>%
    mutate(variant_type = fct_recode(variant_type, 'Synonymous' = 'syn', 'Missense' = 'mis', 'LoF' = 'lof'),
           obs_exp = fct_recode(obs_exp, 'Expected' = 'exp', 'Observed' = 'obs'))
  
  all_pops %>%
    filter(func == 'sum' & pop != 'global' & downsampling <= 8128 &
             # obs_exp == 'n'
             obs_exp == 'Observed' & variant_type == 'LoF'
    ) %>%
    ggplot + aes(x = downsampling, y = count, color = pop) +
    geom_line() + scale_color_manual(values=pop_colors, labels=pop_names)
  
  pop_ds_data = all_pops %>%
    filter(func == 'sum' & 
             # pop != 'global' & 
             # downsampling <= 8128 &
             obs_exp == 'n'  # includes indels
           # obs_exp == 'Observed' & variant_type == 'LoF'  # constraint-variants only
    ) %>%
    group_by(pop) %>%
    arrange(downsampling) %>%
    mutate(variant_diff=count-lag(count),
           sample_diff=downsampling-lag(downsampling),
           variants_per_sample = variant_diff / sample_diff) %>%
    filter(sample_diff != 40)
  
  top = pop_ds_data %$% max(variants_per_sample)
  right = pop_ds_data %$% max(downsampling)
  
  pop_colors[['global']] = 'black'
  pop_names[['global']] = 'gnomAD'
  p = pop_ds_data %>%
    ggplot + aes(x = downsampling, y = variants_per_sample, color = pop) +
    geom_line(size=1) + xlab('Sample size') + theme_classic() + 
    ylab('New pLoF variants\nper additional sample') + coord_cartesian(ylim=c(0, top)) +
    scale_color_manual(values=pop_colors, guide=F)
  
  current = 1
  for(i in pop_ds_data %$% unique(pop)) {
    p = p + annotate('text', x = right, hjust = 1, y = top * current,
                     size = 3,
                     label = pop_names[[i]], color=pop_colors[[i]])
    current = current - 0.07
  }
  
  if (plot_log) {
    p = p + scale_x_log10(labels=comma, breaks=ds_breaks_log)
  } else {
    p = p + scale_x_continuous(labels=comma, breaks=ds_breaks)
  }
  
  if (save_plot) {
    pdf('e4d_variants_per_new_sample.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

indel_summary = function() {
  genome_indel_ds_file = get_or_download_file('indels_summary_genomes.downsampling.txt.bgz', subfolder = 'summary_results/')
  genome_indel_ds_file = read_delim(gzfile(genome_indel_ds_file), delim='\t')
  genome_indel_ds_file %>% 
    filter(worst_csq == 'intergenic_variant' & abs(indel_length) <= 10) %>% 
    group_by(indel_length, downsampling) %>% 
    summarize(singletons=sum(singletons), observed=sum(observed)) %>% 
    filter(downsampling == 15708) %>% ungroup %>%
    ggplot + aes(x = indel_length, y = singletons/observed) +
    geom_bar(stat='identity')
}

site_statistics = function(data_type) {
  obs_poss_sites = load_observed_possible_sites_data(data_type)
  
  obs_poss_sites %>%
    filter(coverage >= 0) %>%
    summarize(observed = sum(n_sites_observed), possible = sum(n_sites_possible),
              prop = observed / possible) %>%
    print
  
  obs_poss_sites %>%
    filter(coverage >= 30) %>%
    summarize(observed = sum(n_sites_observed), possible = sum(n_sites_possible),
              prop = observed / possible) %>%
    print
  
  obs_poss_sites %>%
    filter(coverage >= 0) %>%
    arrange(desc(coverage)) %>%
    mutate(n_sites_observed_c = cumsum(n_sites_observed),
           n_sites_possible_c = cumsum(n_sites_possible),
           prop_observed_c = n_sites_observed_c / n_sites_possible_c) %>%
    ggplot + aes(x = coverage, y = prop_observed_c) + 
    geom_line(size = 2) + 
    xlab('Coverage >= X') + ylab('Cumulative proportion observed')
}
