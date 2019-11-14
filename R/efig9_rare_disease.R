source('fig5_disease.R')
gene_data = load_constraint_data()

constraint_by_omim_discovery = function(save_plot=F, histograms=T) {
  omim_data = load_omim_by_year_data()
  
  gene_data %>%
    left_join(omim_data, by=c('gene_id' = 'Ensembl.Gene.ID')) %>% 
    filter(!is.na(year)) %>%
    mutate(discoverybyNGS=fct_recode(as.character(discoverybyNGS), 
                                     `WES/WGS` = 'TRUE', conventional = 'FALSE'),
           discoverybyNGS=fct_relevel(discoverybyNGS, 'WES/WGS')) %>%
    skim_tee() -> constraint_with_year  # 3343 genes discovered
  
  constraint_with_year %>%
    mutate(NGS=discoverybyNGS == 'WES/WGS') %>%
    summarize(ttest = list(t.test(oe_lof_upper ~ NGS))) %>%
    tidy(ttest) %>%
    print
  
  constraint_with_year %>%
    mutate(NGS=discoverybyNGS == 'WES/WGS') %>%
    summarize(ttest = list(glm(NGS ~ oe_lof_upper + cds_length, family='binomial'))) %>%
    tidy(ttest) %>%
    print
  
  constraint_with_year %>%
    ggplot + aes(x = year, y = oe_lof_upper) + geom_point() + geom_smooth()
  
  constraint_with_year %>%
    ggplot + aes(x = year, y = oe_lof_upper, group = year) + 
    geom_boxplot(fill='lightblue') + xlab('Year') + ylab('LOEUF') + theme_classic()
  
  if (histograms) {
    p = constraint_with_year %>%
      filter(exp_lof > 10) %>%
      ggplot + aes(x = oe_lof_upper, y = discoverybyNGS, fill = ..x..) + 
      geom_density_ridges_gradient(scale = 2, rel_min_height = 0.03, gradient_lwd = 1.) +
      theme_ridges(font_size = 11, center_axis_labels=T) +
      scale_fill_gradientn(colors = gradient_colors2, values = rescale(gradient_values2), guide='none') +
      xlab(constraint_metric_name) + ylab(NULL)
  } else {
    p = constraint_with_year %>%
      filter(exp_lof > 10) %>%
      count(oe_lof_upper_bin, discoverybyNGS) %>%
      add_count(discoverybyNGS, wt=n, name='nn') %>%
      mutate(prop_in_bin = n / nn,
             sem = 1.96 * sqrt(prop_in_bin * (1 - prop_in_bin) / nn)) %>%
      ggplot + aes(x = oe_lof_upper_bin, y = prop_in_bin, ymin = prop_in_bin - sem, ymax = prop_in_bin + sem, color = discoverybyNGS) + 
      geom_pointrange(position = position_dodge(width = 0.6)) + theme_classic() +
      oe_x_axis + scale_y_continuous(labels=percent_format(accuracy=1), name='Percent of OMIM\ngenes by technology') +
      scale_color_manual(values=c('conventional' = 'darkorchid', 'WES/WGS' = color_lof),
                         name='Technology') +
      theme(legend.justification = c(1, 1), legend.position = c(1, 1))
  }
  
  if (save_plot) {
    pdf('constraint_omim_by_tech.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

omim_over_time = function(save_plot=F) {
  omim_data = load_omim_by_year_data()
  
  rounding = 5
  gene_data %>%
    left_join(omim_data, by=c('gene_id' = 'Ensembl.Gene.ID')) %>% 
    filter(!is.na(year)) %>%
    mutate(rounded_year = round((year) / rounding) * rounding) %>%
    skim_tee() -> constraint_with_year  # 3343 genes discovered
  
  # constraint_with_year %>%
  #   filter(exp_lof > 10) %>%
  #   ggplot + aes(x = oe_lof_upper, y = year_group, fill = ..x..) + 
  #   geom_density_ridges_gradient(scale = 10 / rounding, rel_min_height = 0.03, gradient_lwd = 1.) + theme_ridges() +
  #   scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values), guide = 'none') +
  #   xlab('Observed/expected (upper)') + ylab(NULL)
  # 
  # constraint_with_year %>%
  #   filter(exp_lof > 10) %>%
  #   ggplot + aes(x = oe_lof_upper, group = year_group, fill = year_group) + 
  #   geom_density(alpha=0.3) + theme_classic()
  # 
  # constraint_with_year %>%
  #   filter(exp_lof > 10) %>%
  #   ggplot + aes(x = oe_lof_upper, group = year_group, fill = year_group) + 
  #   geom_histogram(position='dodge') + theme_classic()
  # 
  # constraint_with_year %>%
  #   ggplot + aes(x = oe_lof_upper, y = year, group = year) +
  #   geom_density_ridges(scale = 10, rel_min_height = 0.03) + theme_ridges() +
  #   scale_y_reverse(breaks=c(2015, 2010, 2005, 2000, 1995, 1990), expand = c(0.01, 0))
  # 
  # constraint_with_year %>%
  #   filter(exp_lof > 10) %>%
  #   group_by(rounded_year, oe_lof_upper_bin) %>%
  #   summarize(n = n()) %>%
  #   ggplot + aes(x = rounded_year, y = n, fill = oe_lof_upper_bin) + geom_bar(stat='identity') + theme_classic() +
  #   scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values))
  
  p = constraint_with_year %>%
    filter(exp_lof > 10) %>%
    group_by(rounded_year, oe_lof_upper_bin) %>%
    summarize(n = n()) %>%
    group_by(rounded_year) %>%
    mutate(norm_n = n / sum(n)) %>%
    ggplot + aes(x = rounded_year, y = norm_n, fill = oe_lof_upper_bin + 1) + geom_bar(stat='identity') + theme_classic() +
    scale_x_continuous(breaks=seq(1990, 2015, 5)) +
    scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values), breaks=seq(1, 10, 2), name='Constraint\ndecile\n(lower = more\nconstrained)') +
    scale_y_continuous(labels=percent, name='Percent of genes') + xlab('Year (rounded down to nearest 5)') +
    theme(legend.title = element_text(size = 9),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  if (save_plot) {
    pdf('e8a_lof_freq.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

scz_rare_variants = function(save_plot=F, csqs_to_plot=c('Synonymous', 'pLoF')) {
  fname = get_or_download_file('scz.summary.constraint_bins.110919.tsv', 'misc_files/')
  
  scz_data = read_tsv('data/scz.summary.constraint_bins.110919.tsv') %>%
    mutate(csq=fct_relevel(fct_recode(csq, Missense = 'missense', Synonymous = 'synonymous', pLoF = 'LOF')))
  
  p = scz_data %>%
    filter(csq %in% csqs_to_plot) %>%
    ggplot + aes(x = constraint_bin, y = OR, colour = csq, ymax = top_limit, ymin = bottom_limit) +
    geom_point(position = position_dodge(width = 0.6), size = 3) + theme_classic() +
    geom_errorbar(position = position_dodge(width = 0.6), width = 0) +
    scale_color_manual(values=colors, guide=F) +
    ylab('Odds Ratio for\nschizophrenia cases \ncompared to controls') +
    oe_x_axis + geom_hline(yintercept = 1, linetype = "dashed")
  
  if (save_plot) {
    pdf('e9d_scz.pdf', height=3, width=4)
    print(p)
    dev.off()
    png('e9d_scz.png', height=3*300, width=4*300, res=300)
    print(p)
    dev.off()
  }
  return(p)
}

proportion_in_omim = function(save_plot=F) {
  omim_data = load_omim_by_year_data()
  gene_omim_data = gene_data %>%
    mutate(in_omim = gene_id %in% omim_data$Ensembl.Gene.ID)
  
  gene_omim_data %>%
    summarize(ttest = list(t.test(oe_lof_upper ~ in_omim))) %>%
    tidy(ttest)
  
  gene_omim_data %>%
    summarize(ttest = list(glm(in_omim ~ oe_lof_upper + cds_length, family='binomial'))) %>%
    tidy(ttest)
  
  p = gene_omim_data %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(num_in_omim = sum(in_omim, na.rm=T), n=n(), 
              prop_in_omim = num_in_omim / n, sem = sqrt(prop_in_omim * (1 - prop_in_omim) / n),
              prop_upper = prop_in_omim + 1.96 * sem, prop_lower = prop_in_omim - 1.96 * sem) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_in_omim, ymin = prop_lower, ymax = prop_upper) + 
    geom_pointrange() + scale_y_continuous(labels=percent_format(accuracy=1)) +
    theme_classic() + oe_x_axis + ylab('Proportion in OMIM')
  
  if (save_plot) {
    pdf('e9_proportion_in_omim.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

proportion_in_omim_powered = function(save_plot=F) {
  omim_data = load_omim_by_year_data()
  gene_omim_data = gene_data %>%
    filter(!is.na(exp_lof)) %>%
    mutate(in_omim = gene_id %in% omim_data$Ensembl.Gene.ID,
           powered=exp_lof > 10)
  
  # gene_omim_data %>%
  #   ggplot + aes(x = in_omim, y = exp_lof) + 
  #   geom_boxplot() + scale_y_log10()
  
  gene_omim_data %>% count(powered, in_omim) %>% print
  
  gene_omim_data %>% count(powered, in_omim) %$% matrix(n, nrow = 2, ncol = 2) %>% fisher.test

  gene_omim_data %>% left_join(load_all_gene_list_data()) %>%
    filter(gene_list == 'ACMG_2_0' & !powered) %>%
    select(gene, obs_lof, exp_lof, oe_lof, oe_lof_upper, pLI) %>% print
}

efigure9 = function() {
  e9a = proportion_in_omim()
  e9b = constraint_by_omim_discovery(histograms = F)
  e9c = plot_rare_disease(phenotype = 'ASD')
  e9d = scz_rare_variants()
  pdf('extended_data_figure9.pdf', height=6, width=7.5)
  print(ggarrange(e9a, e9b, e9c, e9d, nrow = 2, ncol = 2, align = 'hv', labels = 'auto', vjust = 1))
  dev.off()
  png('extended_data_figure9.png', height=6, width=7.5, units = 'in', res=300)
  print(ggarrange(e9a, e9b, e9c, e9d, nrow = 2, ncol = 2, align = 'hv', labels = 'auto', vjust = 1))
  dev.off()
}


