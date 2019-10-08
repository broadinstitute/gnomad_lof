gene_data = load_constraint_data()
transcript_data = load_constraint_data(level = 'transcript')
tx_summary = load_tx_summary()
tx_disease_gene = load_gene_tissue_mapping()
ds_data = load_downsampled_gene_data()

canonical_most_often_constrained = function(save_plot=F) {
  transcript_data = transcript_data %>%
    filter(!is.na(oe_lof_upper)) %>%
    group_by(gene) %>%
    mutate(most_constrained_transcript = min(oe_lof_upper) == oe_lof_upper) %>% ungroup
  
  odds_ratio = transcript_data %$% 
    table(canonical, most_constrained_transcript) %>%
    fisher.test %$% estimate[[1]]
  
  permute_or = function() {
    transcript_data %>% 
      group_by(gene) %>%
      mutate(most_constrained_transcript = sample(most_constrained_transcript)) %>%
      ungroup %>%
      summarize(both=sum(canonical & most_constrained_transcript, na.rm=T),
                constrained_not_canonical=sum(!canonical & most_constrained_transcript, na.rm=T),
                canonical_not_constrained=sum(canonical & !most_constrained_transcript, na.rm=T),
                neither=sum(!canonical & !most_constrained_transcript, na.rm=T),
                or=both * neither/(constrained_not_canonical * canonical_not_constrained)) %$% or
  }
  or_background = pbreplicate(100, permute_or())
  
  background_color = 'steelblue'
  data.frame(x = or_background) %>%
    ggplot + aes(x = x) +
    geom_density(fill = background_color, alpha = 0.2, color = background_color) +
    theme_classic() + scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(name='Odds ratio of canonical\nand most constrained') +
    geom_vline(xintercept = odds_ratio, linetype='dashed', color = color_lof) +
    annotate('text', 5, 7, label = 'background', color = 'steelblue') +
    annotate('text', 5, 6.5, label = 'actual', color = color_lof)
  # Likely confounded by transcript length
 }

number_of_transcripts_by_constraint = function(save_plot=F) {
  data = transcript_data %>%
    count(gene_id) %>%
    rename(n_transcripts = n) %>%
    inner_join(gene_data) %>%
    left_join(tx_summary)
  
  # p = data %>%
  #   mutate(n_transcripts = ifelse(n_transcripts > 10, 10, n_transcripts)) %>%
  #   ggplot + aes(x = n_transcripts, y = oe_lof_upper, group = n_transcripts) + 
  #   geom_boxplot() + theme_classic() +
  #   scale_x_continuous(breaks=1:10, name='Number of transcripts') +
  #   ylab('Observed/expected ratio (upper)')
  
  p = data %>%
    ggplot + aes(x = oe_lof_upper_bin, y = n_transcripts, group = oe_lof_upper_bin, fill = oe_lof_upper_bin) + 
    geom_boxplot() + theme_classic() + oe_x_axis + coord_cartesian(ylim=c(1, 10)) +
    scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values), guide = FALSE) +
    ylab('Number of transcripts')
  
  data %$%
    glm(oe_lof_upper ~ max_expression + cds_length + n_transcripts) %>%
    summary
  data %$%
    glm(oe_lof_upper ~ cds_length + n_transcripts + mean_expression) %>%
    summary
  
  if (save_plot) {
    pdf('efig8_number_of_transcripts_by_constraint.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

expression_by_constraint = function(save_plot=F, only_canonical=F) {
  all_tx = transcript_data %>%
    left_join(tx_summary)
  
  plot_data = all_tx %>%
    mutate(mean_expression=ifelse(mean_expression > 5, 5, mean_expression)) %>%
    filter(!only_canonical | canonical) 
  
  collapsed_tx = plot_data %>% 
    group_by(oe_lof_upper_bin) %>% 
    summarize(mean_mean_expression = mean(mean_expression, na.rm=T),
              sd_mean_expression = sd(mean_expression, na.rm=T),
              n = n(),
              sem_mean_expression = 1.96 * sd_mean_expression / sqrt(n))
  
  p = ggplot() +
    geom_violin(aes(x = oe_lof_upper_bin, y = mean_expression, group = oe_lof_upper_bin),
                data = plot_data,
                fill = color_lof, alpha = 0.2, color = F) +
    theme_classic() + oe_x_axis +
    geom_pointrange(aes(x = oe_lof_upper_bin, y = mean_mean_expression, 
                        ymin = mean_mean_expression - sem_mean_expression, 
                        ymax = mean_mean_expression + sem_mean_expression), 
                    data = collapsed_tx, color = color_lof) +
    ylab('Mean expression\n(capped at 5)')
  
  if (save_plot) {
    pdf('efig8_expression_by_constraint.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

number_of_go_terms_by_constraint = function(save_plot=F) {
  go_data = load_go_data()
  go_constraint = gene_data %>%
    left_join(go_data, by=c('gene' = 'DB_Object_Symbol')) %>%
    count(oe_lof_upper_bin, gene, gene_length)
  
  mean_go_data = go_constraint %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(mean_go = mean(n, na.rm=T))
  
  # summary(lm(n ~ oe_lof_upper_bin + gene_length, go_constraint))
  
  p = go_constraint %>%
    mutate(n = ifelse(n > 50, 50, n)) %>%
    ggplot + aes(x = oe_lof_upper_bin, y = n, group = oe_lof_upper_bin) + 
    geom_violin(alpha=0.5, fill=color_lof, color=F) + theme_classic() +
    geom_point(aes(y = mean_go), data = mean_go_data, color=color_lof) +
    oe_x_axis + ylab('Number of GO terms')
  
  if (save_plot) {
    pdf('efig8_mean_number_of_go_terms.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

get_pharos_categories = function() {
  fname = '../misc_files/pharos_data.tsv'
  if (file.exists(fname)) {
    pharos_data = read.delim(fname, sep = '\t')
  } else {
    pharos.ids <- gene_data$gene %>% unique
    pharos.ids <- paste0("'", pharos.ids, "'")
    pharos.db = dbConnect(dbDriver("MySQL"), user='tcrd', dbname='tcrd540',
                          host='tcrd.kmc.io')
    query <- paste0('SELECT id, sym FROM protein WHERE sym IN (', 
                    paste0(pharos.ids, collapse = ', '), ');')
    ids <- dbGetQuery(pharos.db, query) %>%
      distinct()
    
    query2 <-  paste0('SELECT id, tdl FROM target WHERE id IN (', 
                      paste0(ids$id, collapse = ', '), ');')
    
    tdl <- dbGetQuery(pharos.db, query2)
    
    dbDisconnect(pharos.db)
    
    pharos_data = tdl %>% left_join(ids, by='id')
    write.table(pharos_data, fname, quote=F, row.names=F, sep='\t')
  }
  return(pharos_data)
}
tdl = get_pharos_categories()
  
functional_categorization = function(save_plot=F) {
  tdl_gene_data = tdl %>%
    left_join(dplyr::select(gene_data, oe_lof_upper_bin, gene), by=c('sym' = 'gene'))
  
  rename_list = c(
    'approved drugs' = 'Tclin',
    'ChEMBL drug activities' = 'Tchem',
    'weaker drug activites' = 'Tbio',
    'no known ligands' = 'Tdark'
  )
  tdl_name_colors = c(
    'approved drugs' = 'firebrick',  # color_lof,
    'ChEMBL drug activities' = 'tomato',  # color_dominant,
    'weaker drug activites' = 'tan1',  # color_recessive,
    'no known ligands' = color_syn # color_benign
  )
  tdl.freqs <- tdl_gene_data %>%
    mutate(tdl=fct_recode(tdl, !!!rename_list),
           tdl=fct_relevel(tdl, names(rename_list))) %>%
    count(oe_lof_upper_bin, tdl) %>% 
    na.omit() %>%
    group_by(tdl) %>%
    mutate(prop = n/sum(n))
  
  # Plotting frequencies
  # p = ggplot(tdl.freqs, aes(oe_lof_upper_bin, n)) +
  #   geom_bar(aes(fill = tdl), position = "dodge", stat="identity") +
  #   # geom_line(aes(color = tdl), size = 2) +
  #   theme_classic() +
  #   oe_x_axis
  
  # Plotting as percentage of group
  p = ggplot(tdl.freqs, aes(oe_lof_upper_bin, prop)) +   
    geom_bar(aes(fill = tdl), position = "dodge", stat="identity") +
    # geom_line(aes(color = tdl), size = 2) +
    # labs(fill='Targets with:') +
    scale_fill_manual(values=tdl_name_colors, name='Targets with:') +
    scale_y_continuous(name = 'Percentage of category', limits = c(0, 0.2), labels=percent_format(accuracy = 1)) +
    theme_classic() +
    oe_x_axis
  
  if (save_plot) {
    pdf('efig8_functional_categories.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

mean_number_of_diseases = function(save_plot=F) {
  omim = load_omim_data()
  
  collapsed_disease_gene = gene_data %>% 
    left_join(omim) %>%
    mutate(associated=!is.na(phenotype)) %>%
    group_by(oe_lof_upper_bin, gene) %>%
    summarize(num_diseases=sum(associated))
  
  collapsed_disease_bins = collapsed_disease_gene %>%
    group_by(oe_lof_upper_bin) %>%
    summarize(mean_num_diseases = mean(num_diseases, na.rm=T),
              sd_num_diseases = sd(num_diseases, na.rm=T),
              n=n())
  
  p = collapsed_disease_bins %>%
    ggplot() + aes(x = oe_lof_upper_bin, y = mean_num_diseases, 
                   ymin = mean_num_diseases - 1.96 * sd_num_diseases / sqrt(n), 
                 ymax = mean_num_diseases + 1.96 * sd_num_diseases / sqrt(n)) +
    theme_classic() + oe_x_axis +
    geom_pointrange(color = color_lof) +
    ylab('Mean number of\ndiseases associated')
  
  if (save_plot) {
    pdf('efig8_mean_number_of_diseases.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

permute_most_constrained = function(tx_disease_gene) {
  tx_disease_gene %>%
    mutate(most_constrained_transcript = sample(most_constrained_transcript)) %>%
    filter(most_expressed) %$%
    mean(most_constrained_transcript)
}

percent_most_expressed_also_most_constrained = function(save_plot=F, force_permutation=F) {
  fname = 'most_constrained_transcript.RData'
  if (!force_permutation) {
    fname = get_or_download_file(fname, subfolder = 'misc_files/')
  } else {
    set.seed(42)
    most_constrained_background_distribution = replicate(10000, permute_most_constrained(tx_disease_gene))
    save(most_constrained_background_distribution, file=paste0('data/', fname))
  }
  most_constrained_background_distribution = get(load(paste0('data/', fname)))
  tx_disease_gene %>%
    # group_by(gene) %>%
    filter(!is.na(most_expressed)) %$% # & n() > 8) %>% ungroup %$% 
    table(most_constrained_transcript, most_expressed) %T>%
    {fisher.test(.) %>% print}
  
  prop_most_constrained = tx_disease_gene %>%
    filter(most_expressed) %$%
    mean(most_constrained_transcript)
  
  background_color = 'steelblue'
  p = data.frame(x = most_constrained_background_distribution) %>%
    ggplot + aes(x = x) +
    geom_density(fill = background_color, alpha = 0.2, color = background_color) +
    theme_classic() + scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(labels=percent_format(accuracy = 1), 
                       name='Percent of most expressed transcripts\nthat are also most constrained') +
    geom_vline(xintercept = prop_most_constrained, linetype='dashed', color = color_lof) +
    annotate('text', 0.4, 80, label = 'background', color = 'steelblue') +
    annotate('text', 0.4, 72, label = 'actual', color = color_lof)
  
  # In cases where a gene-tissue mapping is known for at least one disease,
  # the transcript with the highest expression is also the most constrained
  # 55.4% of the time (OR = 5.4; Fisher's exact p < 1e-100).
  
  # No transcript/tissue pairs found more than once
  # tx_disease_gene %>% filter(!is.na(most_expressed)) %>% count(gene_id, transcript, tissue) %>% filter(n > 1) %>% nrow
  
  if (save_plot) {
    pdf('e8_percent_most_expressed_also_most_constrained.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

load_pop_specific_constraint = function(pops_to_use=c('afr', 'nfe')) {
  pop_counts = ds_data %>% 
    group_by(pop) %>%
    summarize(pop_size=max(downsampling))
  comp_size = pop_counts %>% filter(pop %in% pops_to_use) %$% min(pop_size)
  print(paste('Running', paste(pops_to_use, collapse=', '), 'at n =', comp_size))
  
  pop_data = ds_data %>%
    filter(canonical & downsampling == comp_size & pop %in% pops_to_use) %>%
    transmute(gene, transcript, pop, obs_lof, exp_lof, oe_lof = obs_lof / exp_lof, caf)
  
  min_exp_by_transcript = pop_data %>%
    group_by(transcript) %>%
    summarize(min_exp=min(exp_lof),
              max_oe_lof=max(oe_lof)) %>% ungroup
  
  pop_data %>% 
    inner_join(min_exp_by_transcript) %>%
    return
}
 
plot_pop_specific_constraint_global = function(save_plot=F, pops_to_use=c('afr', 'amr', 'eas', 'nfe', 'sas')) {
  exp_cutoff = 10
  pop_data = load_pop_specific_constraint(pops_to_use=pops_to_use) %>%
    filter(min_exp >= exp_cutoff & max_oe_lof < 1.5) %>%
    rowwise %>%
    do(cbind(., generate_ci(.$obs_lof, .$exp_lof))) %>%
    ungroup
  
  print(paste(nrow(count(pop_data, transcript)), 'genes with >=', exp_cutoff, 'expected LoFs'))
  
  p = pop_data %>%
    group_by(pop) %>%
    summarize(mean_oe_lof = mean(upper, na.rm=T),
              sd_oe_lof = sd(upper, na.rm=T),
              n = n(), sem_oe_lof = 1.96 * sd_oe_lof / sqrt(n)) %>%
    ggplot + aes(x = pop, y = mean_oe_lof, color = pop,
                 ymin = mean_oe_lof - sem_oe_lof, ymax = mean_oe_lof + sem_oe_lof) + 
    geom_pointrange() + theme_classic() + 
    scale_x_discrete(labels=gsub('/', '/\n', pop_names), name=NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(0, 5.5, 0, 20.5)) +
    scale_color_manual(values=pop_colors, guide=F) +
    ylab(constraint_metric_name)
  
  # pop_data %>%
  #   ggplot + aes(x = upper, fill = pop) + 
  #   geom_density(alpha=0.2) + theme_classic() + 
  #   scale_fill_manual(values=pop_colors)
  
  if (save_plot) {
    pdf('efig8_pop_comparison_global.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

plot_pop_specific_constraint_comparison = function(save_plot=F, pops_to_use=c('afr', 'nfe'),
                                                   col_to_use='upper') {
  exp_cutoff = 10
  pop_data = load_pop_specific_constraint(pops_to_use=pops_to_use) %>%
    filter(min_exp >= exp_cutoff & max_oe_lof < 1.5) %>%
    rowwise %>%
    do(cbind(., generate_ci(.$obs_lof, .$exp_lof))) %>%
    ungroup
  
  # pop_data %>%
  #   group_by(gene, transcript) %>%
  #   filter(sum((pop == 'afr')*upper) < sum((pop == 'nfe')*lower) & 
  #            pop %in% c('afr', 'nfe')) %>%
  #   data.frame
  
  spread_pop_data = pop_data %>%
    select(gene, transcript, pop, col_to_use) %>%
    spread(pop, col_to_use)

  print(paste('Got', nrow(spread_pop_data), 'genes'))
  
  print(cor.test(spread_pop_data[[pops_to_use[1]]], spread_pop_data[[pops_to_use[2]]]))
  print(t.test(spread_pop_data[[pops_to_use[1]]], spread_pop_data[[pops_to_use[2]]]))
  print(wilcox.test(spread_pop_data[[pops_to_use[1]]], spread_pop_data[[pops_to_use[2]]]))

  p = spread_pop_data %>%
    ggplot +
    aes_string(x = pops_to_use[1], y = pops_to_use[2]) +
    # geom_point() + aes(gene = gene) +
    geom_hex() +
    theme_classic() + geom_abline(slope = 1, intercept = 0) +
    scale_x_log10() + scale_y_log10() +
    scale_fill_continuous(breaks=seq(5, 20, 5), name='Number\nof genes') +
    xlab(paste0(constraint_metric_name, '\n(', pop_names[pops_to_use[1]], ')')) +
    ylab(paste0(constraint_metric_name, '\n(', pop_names[pops_to_use[2]], ')')) + 
    theme(legend.justification = c(1, 0), legend.position = c(1, 0.01),
          legend.key.size = unit(0.8, "lines"))
  
  if (save_plot) {
    pdf(paste('efig8_pop_comparison', pops_to_use, '.pdf', collapse='_'), height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

efigure8 = function() {
  e8a = functional_categorization()
  e8b = plot_expression_metrics_by_constraint(canonical_all_comparison = T)
  e8c = percent_most_expressed_also_most_constrained()
  e8d = plot_pop_specific_constraint_comparison()
  e8e = plot_pop_specific_constraint_global()
  # e8bottom = ggarrange(e8b, e8c, e8d, e8e, nrow = 2, ncol = 2, align = 'v', labels = c('b', 'c', 'd', 'e'), vjust = 1)
  # ggarrange(e8a, e8bottom, nrow = 2, ncol = 1, align = 'v', labels = c('a', ''), vjust = 1, heights=c(1, 2))
  pdf('extended_data_figure8.pdf', height=8, width=7.5)
  print(ggarrange(e8a, 
            ggarrange(e8b, e8c, e8d, e8e, nrow = 2, ncol = 2, align = 'v', 
                      labels = c('b', 'c', 'd', 'e'), vjust = 1), 
            nrow = 2, ncol = 1, align = 'v', labels = c('a', ''), vjust = 1, heights=c(1, 2)))
  dev.off()
  png('extended_data_figure8.png', height=8, width=7.5, units = 'in', res=300)
  print(ggarrange(e8a, 
            ggarrange(e8b, e8c, e8d, e8e, nrow = 2, ncol = 2, align = 'v', 
                      labels = c('b', 'c', 'd', 'e'), vjust = 1), 
            nrow = 2, ncol = 1, align = 'v', labels = c('a', ''), vjust = 1, heights=c(1, 2)))
  dev.off()
}


