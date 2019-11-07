gene_data = load_constraint_data()
transcript_data = load_constraint_data(level = 'transcript')
sv_data = load_sv_data()

gene_list_spectrum = function(save_plot=F,
                              gene_lists_to_plot=c('Haploinsufficient', 'Autosomal Recessive', 'Olfactory Genes'),
                              gene_lists_to_skip=c('')) {
  or_genes = gene_data %>%
    filter(grepl('^OR', gene)) %>%
    transmute(gene = gene, gene_list = 'Olfactory Genes')

  gene_lists = load_all_gene_list_data() %>%
    bind_rows(or_genes)
  
  # gene_list_spectrum_data = gene_data %>%
  #   mutate(gene_list=strsplit(gene_lists, '\\|')) %>%
  #   unnest(gene_list) %>%
  #   mutate(gene_list = if_else(grepl('HI', gene_list), 'Haploinsufficient', gene_list),
  #          presence = 1) %>%
  #   complete(gene_list, oe_lof_upper_bin, fill = list(presence = 0)) %>%
  #   count(gene_list, oe_lof_upper_bin, wt = presence) %>%
  #   group_by(gene_list) %>%
  #   mutate(prop_in_bin = n / sum(n)) %>% 
  #   ungroup %>%
  #   mutate(gene_list = fct_relevel(gene_list, 'Haploinsufficient')) 
  
  test_data = gene_data %>%
    left_join(gene_lists) %>%
    mutate(gene_list = if_else(grepl('haploinsufficiency', gene_list), 'Haploinsufficient', gene_list))
  
  all_sig_data = map_df(test_data %$% na.omit(unique(gene_list)), function(genel) {
    test_data %>%
      mutate(presence = gene_list == genel) %>%
      do(tidy(glm(presence ~ oe_lof_upper + cds_length, ., family = 'binomial'))) %>%
      mutate(gene_list=genel)
  })
  all_sig_data %>% filter(gene_list %in% c('Olfactory Genes', 'Haploinsufficient')) # , 'all_ar'))
  
  gene_list_spectrum_data = gene_data %>%
    left_join(gene_lists) %>%
    mutate(gene_list = if_else(grepl('haploinsufficiency', gene_list), 'Haploinsufficient', gene_list),
           presence = 1) %>%
    complete(gene_list, oe_lof_upper_bin, fill = list(presence = 0)) %>%
    count(gene_list, oe_lof_upper_bin, wt = presence) %>%
    group_by(gene_list) %>%
    mutate(prop_in_bin = n / sum(n)) %>% 
    ungroup %>%
    mutate(
      gene_list = fct_recode(gene_list, 'Autosomal Recessive' = "all_ar", 'Autosomal Dominant' = "all_ad"),
      gene_list = fct_relevel(gene_list, 'Haploinsufficient'))
  
  gene_list_spectrum_data %>%
    filter(gene_list %in% c('Haploinsufficient', 'Autosomal Dominant', 'Autosomal Recessive', 'Olfactory Genes')) %>%
    filter(!is.na(gene_list) & !is.na(oe_lof_upper_bin)) %>%
    select(-prop_in_bin) %>%
    spread(gene_list, n) %>%
    write_delim('data/s_table_gene_lists.tsv', delim='\t')
  
  top_legend = max(gene_list_spectrum_data$prop_in_bin)# + ko_data$sem)
  gene_list_colors_plot = gene_list_colors
  for (gene_list in gene_lists_to_skip) {
    gene_list_colors_plot[gene_list] = NA
  }
  p3a = gene_list_spectrum_data %>%
    filter(gene_list %in% gene_lists_to_plot) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_in_bin, color = gene_list, fill = gene_list) + 
    # geom_line(lwd=1.5) + 
    geom_bar(position='dodge', stat='identity', width=0.9) + 
    theme_classic() + 
    scale_color_manual(values=gene_list_colors_plot, guide=F) + # name=NULL) +
    scale_fill_manual(values=gene_list_colors_plot, guide=F) + # name=NULL) +
    annotate('text', 4.5, top_legend, hjust=0.5, vjust=1, 
             label='Haploinsufficient', color=gene_list_colors['Haploinsufficient']) +
    annotate('text', 4.5, top_legend*0.88, hjust=0.5, vjust=1, 
             label='Autosomal Recessive', color=gene_list_colors['Autosomal Recessive']) +
    annotate('text', 4.5, top_legend*0.76, hjust=0.5, vjust=1, 
             label='Olfactory Genes', color=gene_list_colors['Olfactory Genes']) +
    ylab('Percent of gene list') + oe_x_axis + scale_y_continuous(labels=percent_format(accuracy = 1)) +
    theme(legend.justification = c(0.5, 0.8), legend.position = c(0.5, 1))
  
  if (save_plot) {
    pdf('3a_gene_list.pdf', height = 3, width = 5)
    print(p3a)
    dev.off()
  }
  return(p3a)
}

acmg_spectrum = function(save_plot=F) {
  acmg_spectrum_data = gene_data %>%
    left_join(load_all_gene_list_data()) %>%
    mutate(presence = 1) %>%
    complete(gene_list, oe_lof_upper_bin, fill = list(presence = 0)) %>%
    count(gene_list, oe_lof_upper_bin, wt = presence) %>%
    group_by(gene_list) %>%
    mutate(prop_in_bin = n / sum(n)) %>% 
    ungroup %>%
    filter(gene_list == 'ACMG_2_0')
  
  p = acmg_spectrum_data %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_in_bin) + 
    geom_bar(position='dodge', stat='identity', width=0.9, fill='darkblue') + 
    theme_classic() + 
    ylab('Percent of ACMG59') + oe_x_axis + 
    scale_y_continuous(labels=percent_format(accuracy = 1))
  
  if (save_plot) {
    pdf('sX_acmg_gene_list.pdf', height = 3, width = 5)
    print(p)
    dev.off()
  }
  return(p)
}

sv_spectrum = function(save_plot=F) {
  gene_plus_sv = gene_data %>%
    left_join(sv_data)
  
  lm(observed_rare_biallelic_LoF_deletions/expected_rare_biallelic_LoF_deletions ~ oe_lof_upper, gene_plus_sv) %>%
    summary
  
  lm(observed_rare_biallelic_LoF_deletions/expected_rare_biallelic_LoF_deletions ~ oe_lof_upper + cds_length, gene_plus_sv) %>%
    summary
  
  gene_plus_sv %>%
    # filter(!is.na(expected_rare_biallelic_LoF_deletions)) %>%
    # mutate(oe_lof_upper_bin=ntile(oe_lof_upper_bin, 10)) %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(n=sum(!is.na(expected_rare_biallelic_LoF_deletions)),
              observed_rare_biallelic_LoF_deletions=sum(observed_rare_biallelic_LoF_deletions, na.rm=T),
              expected_rare_biallelic_LoF_deletions=sum(expected_rare_biallelic_LoF_deletions, na.rm=T),
              oe_sv=observed_rare_biallelic_LoF_deletions/expected_rare_biallelic_LoF_deletions) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = oe_sv) + geom_point(size=3, color=color_lof) + theme_classic() +
    oe_x_axis + geom_hline(yintercept = 1, linetype='dashed', color='darkgray') + 
    ylab('Observed / Expected (SV)')
  
  compute_permutation = function() {
    gene_plus_sv %>%
      group_by(oe_lof_upper_bin) %>%
      sample_n(1000, replace = TRUE) %>%
      summarize(oe_sv=sum(observed_rare_biallelic_LoF_deletions, na.rm=T)/
                  sum(expected_rare_biallelic_LoF_deletions, na.rm=T))
  }
  permuted_data = do.call(rbind, rerun(1000, compute_permutation()))
  
  permutation_summary = permuted_data %>%
    group_by(oe_lof_upper_bin) %>%
    summarise(oe_sv = list(enframe(quantile(oe_sv, probs=c(0.05, 0.5, 0.95))))) %>%
    unnest %>% spread('name', 'value') %>% rename(lower = '5%', mean = '50%', upper = '95%')
  
  p = permutation_summary %>%
    ggplot + aes(x = oe_lof_upper_bin, y = mean, ymin = lower, ymax = upper) +
    geom_pointrange(size=0.5, color=color_lof) + theme_classic() +
    oe_x_axis + geom_hline(yintercept = 1, linetype='dashed', color='darkgray') + 
    ylab('Aggregate deletion\nSV observed/expected') + 
    scale_y_continuous(breaks=seq(0, 1.4, 0.2)) +
    coord_cartesian(ylim=c(0, 1.4))
  
  if (save_plot) {
    pdf('3b_sv.pdf', height=3, width=5)
    print(p)
    dev.off()
  }
  return(p)
}

get_proportion_in_gene_lists = function(gene_data) {
  gene_lists = get_ko_gene_lists()
  gene_lists %>%
    inner_join(gene_data) %>%
    count(gene_list, oe_lof_upper_bin) %>%
    add_count(gene_list, wt = n, name = 'nn') %>%
    mutate(prop_in_bin = n / nn,
           sem = 1.96 * sqrt(prop_in_bin * (1 - prop_in_bin) / nn)) %>%
    return
}

t_test_gene_list = function(gene_data, specific_gene_list) {
  gene_lists = get_ko_gene_lists()
  gene_lists %>%
    filter(gene_list == specific_gene_list) %>%
    right_join(gene_data) %>%
    mutate(in_gene_list = !is.na(gene_list)) %T>%
    {count(., in_gene_list) %>% print} %>%
    summarize(ttest = list(t.test(oe_lof_upper ~ in_gene_list))) %>%
    tidy(ttest)
}

mouse_ko_comparison = function(save_plot=F) {
  gene_data %>%
    left_join(get_ko_gene_lists() %>% filter(gene_list == 'mouse_het_lethal_genes')) %>%
    mutate(mouse_ko = !is.na(gene_list)) %>%
    glm(mouse_ko ~ oe_lof_upper + cds_length, ., family='binomial') %>%
    summary
  
  p = get_proportion_in_gene_lists(gene_data) %>%
    filter(gene_list == 'mouse_het_lethal_genes') %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_in_bin, 
                 ymin = prop_in_bin - sem, ymax = prop_in_bin + sem) + 
    geom_bar(stat='identity', fill = color_lof, width=0.6) +
    theme_classic() + 
    ylab('Percent of mouse\nhet lethal knockout genes') + oe_x_axis + 
    scale_y_continuous(labels=percent_format(accuracy = 1))
  
  t_test_gene_list(gene_data, 'mouse_het_lethal_genes') %>% print
  
  if (save_plot) {
    pdf('3c_mouse_ko.pdf', height = 3, width = 5)
    print(p)
    dev.off()
  }
  return(p)
}

cell_ko_comparison = function(save_plot=F) {
  gene_data %>%
    left_join(get_ko_gene_lists() %>% filter(gene_list == 'CEGv2')) %>%
    mutate(cell_essential = !is.na(gene_list)) %>%
    glm(cell_essential ~ oe_lof_upper + cds_length, ., family='binomial') %>%
    summary
  gene_data %>%
    left_join(get_ko_gene_lists() %>% filter(gene_list == 'NEGv1')) %>%
    mutate(cell_non_essential = !is.na(gene_list)) %>%
    glm(cell_non_essential ~ oe_lof_upper + cds_length, ., family='binomial') %>%
    summary
  
  lists_to_plot = c('CEGv2' = color_lof,
                    'NEGv1' = color_benign)
  lists_labels = c('CEGv2' = 'Cell essential',
                    'NEGv1' = 'Cell non-essential')
  ko_data = get_proportion_in_gene_lists(gene_data) %>%
    filter(gene_list %in% names(lists_to_plot))
  
  t_test_gene_list(gene_data, 'CEGv2') %>% print
  t_test_gene_list(gene_data, 'NEGv1') %>% print
  
  top_legend = max(ko_data$prop_in_bin)# + ko_data$sem)
  p = ko_data %>%
    ggplot + aes(x = oe_lof_upper_bin, y = prop_in_bin, fill = gene_list,
                 ymin = prop_in_bin - sem, ymax = prop_in_bin + sem) + 
    geom_bar(stat = 'identity', position = 'dodge') + 
    theme_classic() + 
    ylab('Percent of essential /\n non-essential genes') + oe_x_axis + 
    scale_y_continuous(labels=percent_format(accuracy = 1)) +
    scale_fill_manual(values=lists_to_plot, guide=F) + # labels=lists_labels, name=NULL) +
    annotate('text', 4.5, top_legend, hjust=0.5, vjust=1, 
             label='Cell essential', color=gene_list_colors['Haploinsufficient']) +
    annotate('text', 4.5, top_legend*0.88, hjust=0.5, vjust=1, 
             label='Cell non-essential', color=gene_list_colors['Olfactory Genes'])
    # theme(legend.justification = c(0.5, 0.8), legend.position = c(0.5, 1))
  
  if (save_plot) {
    pdf('3d_cell_ko.pdf', height = 3, width = 5)
    print(p)
    dev.off()
  }
  return(p)
}


figure3 = function() {
  p3a = gene_list_spectrum()
  p3b = sv_spectrum()
  p3c = mouse_ko_comparison()
  p3d = cell_ko_comparison()
  pdf('figure3.pdf', height=5, width=7.5)
  print(ggarrange(p3a, p3b, p3c, p3d, nrow = 2, ncol = 2, labels = 'auto', align='hv', vjust = 1))
  dev.off()
  png('figure3.png', height=5, width=7.5, units = 'in', res=300)
  print(ggarrange(p3a, p3b, p3c, p3d, nrow = 2, ncol = 2, labels = 'auto', align='hv', vjust = 1))
  dev.off()
}
