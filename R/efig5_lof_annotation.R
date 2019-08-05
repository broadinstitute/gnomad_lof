gene_data = load_constraint_data()

lof_frequency_spectrum = function(save_plot=F, cumulative=F) {
  if (cumulative) {
    p = gene_data %>%
      mutate(p_rank = percent_rank(p)) %>%
      ggplot + aes(x = p, y = p_rank) + geom_line(lwd=1.5, color=color_lof) + theme_classic() + 
      xlab('LoF frequency') + scale_x_log10() +
      ylab('Percentile') + scale_y_continuous(labels=percent)
  } else {
    p = gene_data %>%
      ggplot + aes(x = p) + geom_density(fill = color_lof, color = color_lof) + 
      theme_classic() + scale_x_log10() + xlab('LoF frequency') + ylab('Density')
  }
  
  sum(gene_data$p > 0.001, na.rm=T)
  sum(gene_data$p_afr > 0.001 | gene_data$p_amr > 0.001 | gene_data$p_eas > 0.001 | 
        gene_data$p_nfe > 0.001 | gene_data$p_sas > 0.001, na.rm=T)
  
  if (save_plot) {
    pdf(paste0('e5a_lof_freq', 
               ifelse(cumulative, '_cumulative', ''),
               '.pdf'), height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

selection_coefficient = function(save_plot=F, cumulative=F) {
  if (cumulative) {
    p = gene_data %>%
      mutate(s_rank = percent_rank(s)) %>%
      ggplot + aes(x = s, y = s_rank) + geom_line(lwd=1.5, color=color_lof) + theme_classic() + 
      xlab('Selection coefficient') + scale_x_log10() +
      ylab('Percentile') + scale_y_continuous(labels=percent)
  } else {
    p = gene_data %>%
      ggplot + aes(x = s) + 
      geom_density(fill = color_lof, color = color_lof) + 
      scale_x_log10() + theme_classic() + 
      xlab('Selection coefficient') + ylab('Density')
  }
  
  if (save_plot) {
    pdf(paste0('e5_selection_coeff', 
               ifelse(cumulative, '_cumulative', ''),
               '.pdf'), height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}
s_vs_lof_frequency = function() {
  p = gene_data %>%
    ggplot + aes(x = s, y = p) + 
    geom_point(fill = color_lof, color = color_lof, alpha=0.2) + 
    scale_x_log10() + scale_y_log10() + theme_classic() + 
    xlab('Selection coefficient') + ylab('LoF Frequency')
  library(cowplot)
  xp = axis_canvas(p, axis = "x") +
    geom_density(data=gene_data, aes(x = log10(s)), 
                 fill=color_lof, color=color_lof, alpha=0.5)
  yp = axis_canvas(p, axis = "y", coord_flip = T) +
    geom_density(data=gene_data, aes(x = log10(p)),
                 fill=color_lof, color=color_lof, alpha=0.5) +
    coord_flip()
  p1 = insert_xaxis_grob(p, xp, grid::unit(.2, "null"), position = "top")
  p2 = insert_yaxis_grob(p1, yp, grid::unit(.2, "null"), position = "right")
  ggdraw(p2) %>%
    return
  # ggarrange(xp, NULL, p, yp, align='hv', widths = c(2, 1), heights = c(1, 2))
}

expected_lof_spectrum = function() {
  gene_data %>%
    ggplot + aes(x = exp_lof) + geom_density(fill = color_lof, color = color_lof) + 
    theme_classic() + scale_x_log10(labels=comma_format(accuracy=0.1)) +
    xlab('Expected LoFs') + ylab('Density')
  
  gene_data %>%
    mutate(proportion = percent_rank(exp_lof)) %>%
    ggplot + aes(x = exp_lof, y = proportion) + 
    geom_line(color = color_lof, lwd=2) + theme_classic() +
    scale_x_log10(labels=comma_format(accuracy=0.1)) +
    scale_y_continuous(labels=percent) +
    xlab('Expected LoFs') + ylab('Cumulative proportion of genes')
}

caf_vs_sample_size = function() {
  ds_melt = load_downsampled_data()
  ds_melt %>%
    filter(obs_exp == 'caf' & variant_type == 'sum') %>%
    ggplot + aes(x = downsampling, y = count) + geom_line() + 
    theme_classic() + downsampling_x_axis()
  
  ds_data = load_downsampled_gene_data()
  ds_data %>%
    filter(canonical) %>%
    head(2e5) %>%
    ggplot + aes(x = downsampling, y = caf, group = transcript) + 
    geom_line(alpha=0.2) + theme_classic() + 
    scale_y_log10() + downsampling_x_axis(log=F) + ylab('CAF')
}

end_trunc_assessment = function(save_plot=F) {
  # load_maps_data(cut = 'end_trunc_50bp') %>% filter(!is.na(fifty_bp_rule))
  
  fifty_bp = load_maps_data(cut = 'end_trunc_50bp') %>%
    filter(fifty_bp_rule == 'FAIL')
  gerp = load_maps_data(cut = 'end_trunc_gerp')
  
  gerp_plot = gerp %>%
    filter(worst_csq == 'stop_gained' & !is.na(gerp)) %>%
    cumulative_maps(prefix='forward_maps') %>%
    arrange(desc(gerp)) %>%
    cumulative_maps(prefix = 'reverse_maps') %>%
    mutate(proportion_filtered = forward_maps_variant_count/sum(variant_count)) %>%
    filter(maps > 0) %>%
    select(gerp, proportion_filtered, forward_maps, forward_maps_upper, forward_maps_lower,
           reverse_maps, reverse_maps_upper, reverse_maps_lower)
  
  gerp_plot %>%
    arrange(abs(proportion_filtered - 0.5)) %>%
    head(2) %$%
    lm(gerp ~ proportion_filtered) %>% 
    predict(data.frame(proportion_filtered=0.5))  # 180
  
  gerp_plot %>%
    arrange(abs(proportion_filtered - 0.5)) %>%
    head(2)
  
  p = gerp_plot %>%
    ggplot + aes(x = proportion_filtered) + 
    geom_pointrange(aes(y = forward_maps, ymin = forward_maps_lower, ymax = forward_maps_upper), color = '#1F77B4') + 
    geom_pointrange(aes(y = reverse_maps, ymin = reverse_maps_lower, ymax = reverse_maps_upper), color = 'darkred') + 
    theme_classic() +
    geom_hline(yintercept = fifty_bp$maps, linetype = 'dashed', color = 'darkgray') + 
    scale_x_continuous(labels=percent_format(accuracy=1)) +
    xlab('Percent filtered') + ylab('MAPS')
  
  if (save_plot) {
    pdf('e5_end_trunc.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

loftee_comparisons = function() {
  pre_loftee_data = load_constraint_data(loftee=F) %>%
    transmute(gene, transcript, 
              p_no_loftee=p, classic_caf_no_loftee=classic_caf,
              n_sites_no_loftee=n_sites, max_af_no_loftee=max_af)
  
  gene_data %>%
    left_join(pre_loftee_data) -> loftee_compare_data
  loftee_compare_data %>%
    ggplot + aes(x = p_no_loftee, y = p) + geom_point(alpha=0.2) + 
    theme_classic() + scale_x_log10() + scale_y_log10()
  loftee_compare_data %>%
    ggplot + aes(x = p_no_loftee, y = p, color = oe_lof_upper_bin) + geom_point(alpha=0.2) + 
    theme_classic() + scale_x_log10() + scale_y_log10() +
    scale_color_gradientn(colors = gradient_colors, values = rescale(gradient_value_deciles), 
                          guide='none')
  
  loftee_compare_data %>%
    ggplot + aes(x = oe_lof_upper, y = p_no_loftee) + 
    geom_hex() + theme_classic() + scale_y_log10()
  
  loftee_compare_data %>%
    ggplot + aes(x = oe_lof_upper, y = p) + 
    geom_hex() + theme_classic() + scale_y_log10()
  
  loftee_compare_data %>%
    ggplot + aes(x = oe_lof_upper, y = p_no_loftee) + 
    geom_point() + theme_classic() + scale_y_log10()
  
  loftee_compare_data %>%
    ggplot + aes(x = oe_lof_upper, y = p) + 
    geom_point() + theme_classic() + scale_y_log10()
  
  loftee_compare_data %>%
    ggplot + aes(x = oe_lof_upper, y = max_af_no_loftee) + geom_point() +
    theme_classic() + scale_y_log10() 
  
  loftee_compare_data %>%
    ggplot + aes(x = oe_lof_upper, y = max_af) + geom_point() +
    theme_classic() + scale_y_log10() 
  
  loftee_compare_data %>%
    mutate(p = if_else(is.na(p) & !is.na(lof_z), 0, p) + 1e-10,
           p_no_loftee = if_else(is.na(p_no_loftee) & !is.na(lof_z), 0, p_no_loftee) + 1e-10
    ) %>%
    filter(!is.na(p) & !is.na(p_no_loftee)) %>%
    mutate(predicted_p = exp(lm(log(p_no_loftee) ~ oe_lof, .)$fitted.values),
           above_line = p_no_loftee > predicted_p,
           drop = p < p_no_loftee) %$% table(above_line, drop)
  
  loftee_compare_data %>%
    mutate(p = if_else(is.na(p) & !is.na(lof_z), 0, p) + 1e-10,
           p_no_loftee = if_else(is.na(p_no_loftee) & !is.na(lof_z), 0, p_no_loftee) + 1e-10
    ) %>%
    filter(!is.na(p) & !is.na(p_no_loftee)) %>%
    mutate(predicted_p = exp(lm(log(p_no_loftee) ~ oe_lof, .)$fitted.values),
           above_line = p_no_loftee > predicted_p,
           delta_p = p - p_no_loftee) %>%
    ggplot + aes(x = delta_p, group = above_line, fill = above_line) + 
    geom_density(alpha=0.5) + scale_x_log10() + theme_classic()
}

curation_results = function () {
  all_hom_data = read_delim(gzfile('data/all_homozygous_lofs.txt.bgz'), delim = '\t', 
                            col_types = list(locus=col_character()))
  
  curation_data = read_csv('data/Final_gnomad_LOF_curation.csv') %>%
    mutate(verdict=fct_relevel(verdict, 'LoF', 'likely_LoF',
                               'uncertain_LoF', 'likely_not_LoF', 'not_LoF'))
  gene_data = load_constraint_data()
  or_genes = gene_data %>%
    filter(grepl('^OR', gene)) %>%
    transmute(gene = gene, gene_list = 'Olfactory Genes')
  
  gene_lists = load_all_gene_list_data() %>%
    bind_rows(or_genes)
  
  curation_data %>%
    filter(verdict %in% c('likely_LoF', 'LoF')) %>%
    # filter(verdict %in% c('likely_not_LoF', 'not_LoF')) %>%
    select(-id, -variant_id_1, -curator_id) %>%
    summarize_if(is.numeric, sum) %>%
    gather('reason', 'count') %>%
    mutate(reason=fct_reorder(reason, count, .desc = TRUE)) %>%
    ggplot + aes(x = reason, y = count) + 
    geom_bar(stat='identity', fill='darkblue') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(l = 30))
  
  # Raw overall verdict
  curation_data %>%
    count(verdict)
  
  # Summary of similarities/differences
  curation_data_by_variant = curation_data %>%
    count(variant_id, verdict) %>%
    complete(variant_id, verdict, fill=list(n = 0)) %>%
    spread(verdict, n)
  
  curation_data_by_variant %>%
    count(LoF, likely_LoF, uncertain_LoF, likely_not_LoF, not_LoF) %>% 
    arrange(desc(n))
  
  curation_data_collapsed = curation_data %>%
    mutate(verdict=fct_recode(verdict,
                              'not_LoF' = 'likely_not_LoF',
                              'LoF' = 'likely_LoF'
    ),
    verdict=fct_relevel(verdict, 'LoF', 'uncertain_LoF', 'not_LoF'))
  
  # Overall verdict
  curation_data_collapsed %>%
    count(verdict)
  
  # Summarized similarities/differences
  curation_data_by_variant_collapsed = curation_data_collapsed %>%
    count(variant_id, verdict) %>%
    complete(variant_id, verdict, fill=list(n = 0)) %>%
    spread(verdict, n)
  
  curation_data_by_variant_collapsed %>%
    count(LoF, uncertain_LoF, not_LoF) %>% 
    arrange(desc(n))
  
  # Differing variants
  curation_data_by_variant_collapsed %>%
    filter(LoF > 0 & not_LoF > 0)
  
  confident_lofs = curation_data_by_variant_collapsed %>%
    filter(LoF > 0 & not_LoF == 0)
  
  ar_genes = load_all_gene_list_data() %>%
    filter(gene_list == 'all_ar') %$% gene
  
  all_hom_data %>%
    inner_join(confident_lofs) %>%
    count(gene_symbol = csq.gene_symbol) %$% gene_symbol -> hom_ko_genes
  write.table(hom_ko_genes, row.names=F, quote=F, col.names = F, file='data/hom_ko_genes.txt')
  
  all_hom_data %>%
    inner_join(confident_lofs) %>%
    count(gene_id = csq.gene_id) %$% gene_id -> hom_ko_genes
  
  gene_data %>%
    left_join(gene_lists) %>%
    mutate(ko = ifelse(gene_id %in% hom_ko_genes, 'ko', 'no_ko')) %>%
    count(gene_list, ko) %>% 
    spread(ko, n) %>%
    mutate(prop_ko = ko / (ko + no_ko)) %>%
    arrange(desc(prop_ko)) %>% data.frame
  
  gene_data %>%
    mutate(ar_gene = gene %in% ar_genes, 
           ko = gene_id %in% hom_ko_genes) %>%
    count(ar_gene, ko)
  
  # fisher.test(matrix(c(89, 1557, 1090, 16968), nrow=2))
  
  gene_data %>%
    mutate(ar_gene = gene %in% ar_genes, 
           ko = gene_id %in% hom_ko_genes) %>%
    filter(ar_gene & ko) %$% gene
  
  set.seed(42)
  sample_with_matched_loeuf_distr = function(gene_data) {
    data = gene_data %>%
      filter(!is.na(oe_lof_upper_bin)) %>%
      mutate(ar_gene = gene %in% ar_genes)
    
    data %>%
      filter(ar_gene) %$%
      density(oe_lof_upper) -> dens
    data %>%
      filter(!ar_gene) %$%
      density(oe_lof_upper) -> dens_bg
    
    data = data %>%
      mutate(prob = approxfun(dens)(oe_lof_upper) / approxfun(dens_bg)(oe_lof_upper))
    
    # data %>%
    #   filter(!ar_gene) %>%
    #   sample_n(sum(data$ar_gene), replace = F, weight = prob) %>%
    #   union(data %>% filter(ar_gene)) %T>%
    #   {print(group_by(., ar_gene) %>% summarize(mean_oe_lof_upper=mean(oe_lof_upper), n=n()))} %>%
    #   ggplot + aes(x = oe_lof_upper, group = ar_gene, fill = ar_gene) +
    #   geom_density(alpha=0.2) + theme_classic() + scale_x_log10()
    
    data %>%
      sample_n(sum(data$ar_gene), replace = F, weight = prob) %>%
      mutate(ko = gene_id %in% hom_ko_genes) %>%
      summarize(ko_and_ar = sum(ar_gene & ko),
                ko_not_ar = sum(!ar_gene & ko),
                not_ko_ar = sum(ar_gene & !ko),
                not_ko_not_ar = sum(!ar_gene & !ko)) %>%
      do(ft=fisher.test(matrix(c(.$ko_and_ar, .$ko_not_ar, .$not_ko_ar, .$not_ko_not_ar), nrow=2))) %>%
      tidy(ft)
    
  }
  test = bind_rows(pbreplicate(100, sample_with_matched_loeuf_distr(gene_data), simplify = F))
}

variants_per_individual = function(data_type='exomes', save_plot=F) {
  samples = data.frame(datatype='exomes', eas=9197, sas=15308, nfe=56885, fin=10824, asj=5040, afr=8128, amr=17296, oth=3070) %>%
    rbind(data.frame(datatype='genomes', eas=780, nfe=7718, fin=1738, asj=145, afr=4359, amr=424, sas=0, oth=513+31)) %>%
    gather('pop', 'n', -datatype) %>%
    filter(datatype==data_type) %>%
    do(bind_rows(., summarize(., datatype=data_type, pop='global', n=sum(n))))
  
  fname = get_or_download_file(base_fname = paste0('variants_per_sample_', data_type, '.txt.bgz'), 
                               subfolder = 'summary_results/', version = 'v1.1')
  vpi_data = read_delim(gzfile(fname), delim='\t', col_types = cols(lof=col_character()))
  
  proc_vpi_data = vpi_data %>%
    filter(bin != 'Not found') %>%
    mutate(lof=if_else(worst_csq %in% c('frameshift_variant', 'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'),
                       lof, NA_character_)) %>%
    group_by(worst_csq, lof, no_lof_flags, protein_coding, pop, bin, lcr, curated) %>%
    summarize(total=sum(total), total_hom=sum(total_hom)) %>% ungroup
  
  # LoFs per individual
  proc_vpi_data %>%
    filter(protein_coding & bin != '>95%') %>%
    group_by(worst_csq, lof, no_lof_flags, bin, pop, lcr, curated, no_lof_flags) %>%
    summarize(total=sum(total), total_hom=sum(total_hom)) %>% ungroup %>%
    left_join(samples) %>%
    mutate(var_per_ind = total / n,
           hom_per_ind = total_hom / n,
           bin = fct_relevel(as.factor(bin), 'Singleton', 'Doubleton', 'AC 3 - 5',
                             'AC 6 - 0.01%', '0.01% - 0.1%', '< 0.1%', '0.1% - 1%',
                             '1% - 10%', '10% - 95%', '>95%'),
           curated=is.na(curated) | curated) %>%
    filter(pop == 'global' & lof == 'HC') %>%
    arrange(worst_csq, bin) %>%
    group_by(lcr, curated, no_lof_flags) %>%
    summarize(var_per_ind = sum(var_per_ind), hom_per_ind = sum(hom_per_ind)) %>%
    ungroup %>%
    arrange(desc(no_lof_flags), lcr, desc(curated)) %>%
    mutate(var_per_ind = cumsum(var_per_ind), hom_per_ind=cumsum(hom_per_ind))
  
  if (save_plot) {
    pdf('variants_per_individual.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

efigure5 = function() {
  e5b = end_trunc_assessment()
  e5c = lof_frequency_spectrum()
  pdf('extended_data_figure5.pdf', height=3, width=8)
  ggarrange(e5b, e5c, nrow = 1, ncol = 2, align = 'hv', labels = c('b', 'c'), vjust = 1)
  dev.off()
  png('extended_data_figure5.png', height=3, width=8, units = 'in', res=300)
  ggarrange(e5b, e5c, nrow = 1, ncol = 2, align = 'hv', labels = c('b', 'c'), vjust = 1)
  dev.off()
}



