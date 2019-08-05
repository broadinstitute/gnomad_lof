gene_data = load_constraint_data()

plot_rare_disease = function(save_plot=F, phenotype = 'ddid', csqs_to_plot=c('Synonymous', 'pLoF'),
                             expected_cutoff = 10) {
  # Adapted from Beryl Cummings
  assertthat::assert_that(phenotype %in% c('ddid', 'ASD'))
  phenotype_names = c('ddid' = 'ID/DD',
                      'ASD' = 'ASD')
  
  sample_sizes = tribble(~group, ~n_indiv,
                         'ASD', 6430,
                         'ddid', 5305,
                         'Control', 2179)
  
  de_novo_data = get_de_novo_data()
  de_novo_results = de_novo_data %>%
    filter(group %in% c(phenotype, 'Control')) %>%
    mutate(csq = fct_relevel(case_when(csq %in% lof_like & lof == 'HC' ~ 'pLoF',
                                       csq == 'missense_variant' ~ 'Missense',
                                       csq == 'synonymous_variant' ~ 'Synonymous',
                                       TRUE ~ format_vep_category(csq)),
                             'pLoF', 'Missense', 'Synonymous')) %>%
    left_join(gene_data, by=c('ensg' = 'gene_id')) %>%
    filter(exp_lof >= expected_cutoff) %>%
    count(csq, group, symbol, oe_lof_upper_bin) %>%
    group_by(csq, oe_lof_upper_bin, group) %>%
    summarize(n_denovo=sum(n), n_genes=n()) %>%
    left_join(sample_sizes) %>%
    group_by(oe_lof_upper_bin, csq) %>%
    do(ps_tst(.)) %>% ungroup

  de_novo_results = de_novo_results %>%
    filter(csq %in% csqs_to_plot) %>%
    mutate(csq=fct_drop(csq))
  
  ymaxlim = ifelse(phenotype == 'ddid', 16, 7)
  p = de_novo_results %>%
    ggplot + aes(x = oe_lof_upper_bin, y = OR, color = csq, ymax = upper, ymin = lower) + theme_classic() +
    geom_point(position = position_dodge(width = 0.8), size = 3) + 
    geom_errorbar(position = position_dodge(width = 0.8), width = 0) +
    geom_segment(aes(x = oe_lof_upper_bin - 0.2, xend = oe_lof_upper_bin - 0.2, y = OR), 
                 yend = ymaxlim * 1.05, arrow = arrow(length = unit(0.2, "cm")), 
                 data = de_novo_results %>% filter(upper > ymaxlim)) +
    geom_hline(yintercept = 1, linetype = "dashed") + 
    oe_x_axis + ylab(paste('Rate ratio for de novo', 
                           paste('variants in', phenotype_names[[phenotype]], 'cases'),
                           'compared to controls',
                           sep='\n')) +
    coord_cartesian(ylim = c(0, ymaxlim)) +
    scale_color_manual(values = colors, guide=F) +
    geom_text(x = 9, y = ymaxlim * 0.93, vjust=1, hjust=1, label='synonymous', color = color_syn, size = 3) +
    geom_text(x = 9, y = ymaxlim, vjust=1, hjust=1, label='pLoF', color = color_lof, size = 3)
  
  de_novo_results = de_novo_data %>%
    filter(group %in% c(phenotype, 'Control')) %>%
    mutate(csq = fct_relevel(case_when(csq %in% lof_like & lof == 'HC' ~ 'pLoF',
                                       csq == 'missense_variant' ~ 'Missense',
                                       csq == 'synonymous_variant' ~ 'Synonymous',
                                       TRUE ~ format_vep_category(csq)),
                             'pLoF', 'Missense', 'Synonymous'
    )) %>%
    left_join(gene_data, by=c('ensg' = 'gene_id')) %>%
    mutate(percentile=ntile(oe_lof_upper, 100)) %>%
    filter(exp_lof >= expected_cutoff) %>%
    count(csq, group, percentile) %>%
    arrange(percentile) %>%
    complete(csq, group, percentile, fill=list(n=0)) %>%
    group_by(csq, group, percentile) %>%
    mutate(n_denovo=cumsum(n)) %>%
    # summarize(n_denovo=sum(n)) %>%
    left_join(sample_sizes) %>%
    group_by(percentile, csq) %>%
    do(ps_tst(.)) %>% ungroup
  
  # Percentiles
  p_percentile = de_novo_results %>%
    filter(csq %in% csqs_to_plot) %>%
    mutate(OR=ifelse(OR > 100, 100, OR), upper=ifelse(upper > 100, 100, upper)) %>%
    ggplot + aes(x = percentile, y = OR, color = csq, ymax = upper, ymin = lower) + theme_classic() +
    geom_point(position = position_dodge(width = 0.8), size = 3) + 
    geom_errorbar(position = position_dodge(width = 0.8), width = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") + 
    ylab(paste('Rate ratio for de novo', 
                           paste('variants in', phenotype_names[[phenotype]], 'cases'),
                           'compared to controls',
                           sep='\n')) +
    scale_color_manual(values = colors, guide=F) +
    annotate('text', x = 99, y = 93, vjust=1, hjust=1, label='synonymous', color = color_syn, size = 3) +
    annotate('text', x = 99, y = 100, vjust=1, hjust=1, label='pLoF', color = color_lof, size = 3)
  
  # Raw counts
  p_raw = de_novo_data %>%
    filter(group %in% c(phenotype, 'Control')) %>%
    mutate(csq = fct_relevel(case_when(csq %in% lof_like & lof == 'HC' ~ 'pLoF',
                                       csq == 'missense_variant' ~ 'Missense',
                                       csq == 'synonymous_variant' ~ 'Synonymous',
                                       TRUE ~ format_vep_category(csq)),
                             'pLoF', 'Missense', 'Synonymous'
    )) %>%
    left_join(gene_data, by=c('ensg' = 'gene_id')) %>%
    mutate(percentile=ntile(oe_lof_upper, 100)) %>%
    filter(exp_lof >= expected_cutoff) %>%
    count(csq, group, percentile) %>%
    arrange(percentile) %>%
    complete(csq, group, percentile, fill=list(n=0)) %>%
    count(csq, group, percentile, wt = n) %>%
    left_join(sample_sizes) %>%
    filter(csq %in% csqs_to_plot) %>%
    ggplot + aes(x = percentile, color = csq, shape = group, y = nn, fill = csq) + 
    geom_point() + 
    scale_color_manual(values = colors, guide=F) +
    scale_fill_manual(values = colors, guide=F) +
    annotate('text', x = 99, y = 25, vjust=1, hjust=1, label='synonymous', color = color_syn, size = 3) +
    annotate('text', x = 99, y = 30, vjust=1, hjust=1, label='pLoF', color = color_lof, size = 3) +
    theme_classic() + geom_smooth()
  
  # Scaled to n individuals
  scaled_data = de_novo_data %>%
    filter(group %in% c(phenotype, 'Control')) %>%
    mutate(csq = fct_relevel(case_when(csq %in% lof_like & lof == 'HC' ~ 'pLoF',
                                       csq == 'missense_variant' ~ 'Missense',
                                       csq == 'synonymous_variant' ~ 'Synonymous',
                                       TRUE ~ format_vep_category(csq)),
                             'pLoF', 'Missense', 'Synonymous'
    )) %>%
    left_join(gene_data, by=c('ensg' = 'gene_id')) %>%
    mutate(percentile=ntile(oe_lof_upper, 100)) %>%
    filter(exp_lof >= expected_cutoff) %>%
    count(csq, group, percentile) %>%
    arrange(percentile) %>%
    complete(csq, group, percentile, fill=list(n=0)) %>%
    count(csq, group, percentile, wt = n) %>%
    left_join(sample_sizes) %>%
    filter(csq %in% csqs_to_plot) %>%
    mutate(rate_ratio = n / n_indiv)
  
  scaled_data %>%
    ggplot + aes(x = percentile, color = csq, shape = group, y = rate_ratio, fill = csq) + 
    geom_point() + 
    scale_color_manual(values = colors, guide=F) +
    scale_fill_manual(values = colors, guide=F) +
    annotate('text', x = 99, y = 0.0093, vjust=1, hjust=1, label='synonymous', color = color_syn, size = 3) +
    annotate('text', x = 99, y = 0.01, vjust=1, hjust=1, label='pLoF', color = color_lof, size = 3) +
    theme_classic() + geom_smooth()
  
  textp = ifelse(phenotype == 'ddid', 45, 7)
  p_loess_ratio = scaled_data %>%
    nest(-csq, -group) %>%
    mutate(l=map(data, ~ loess(rate_ratio ~ percentile, data = .)),
           augmented=map(l, augment)) %>%
    unnest(augmented) %>%
    group_by(csq, percentile) %>%
    summarize(loess_ratio = sum(.fitted * as.integer(group == phenotype)) /
                sum(.fitted * as.integer(group == 'Control'))) %>%
    ggplot + aes(x = percentile, color = csq, y = loess_ratio, fill = csq) +
    geom_line(size=2) +
    scale_color_manual(values = colors, guide=F) +
    scale_fill_manual(values = colors, guide=F) +
    annotate('text', x = 99, y = 0.93*textp, vjust=1, hjust=1, label='synonymous', color = color_syn, size = 3) +
    annotate('text', x = 99, y = textp, vjust=1, hjust=1, label='pLoF', color = color_lof, size = 3) +
    theme_classic()
  
  if (save_plot) {
    pdf(paste0('rare_disease_', phenotype, '.pdf'), height=3, width=4)
    print(p)
    dev.off()
    png(paste0('rare_disease_', phenotype, '.png'), height=3*300, width=5*300, res=300)
    print(p)
    dev.off()
  }
  return(p)
}

calc95_ci <- function(df) {
  # Adapted from Jack Kosmicki
  # Requires n_denovo and n_indiv in dataset
  df %>%
    rowwise %>% mutate(ci_mult=qt(0.975, n_indiv)) %>%
    ungroup %>%
    mutate(rate = n_denovo / n_indiv,
           SE = sqrt(rate) / sqrt(n_indiv),
           ci = SE * ci_mult
    )
}

ps_tst <- function(f){
  # Adapted from Jack Kosmicki
  controls = subset(f, group == 'Control')
  case = subset(f, group != 'Control')
  ps = poisson.test(c(case$n_denovo, controls$n_denovo), c(case$n_indiv, controls$n_indiv))
  or = ps$estimate
  lower = ps$conf.int[1]
  upper = ps$conf.int[2]
  pval = ps$p.value
  case_control_n = paste(case$n_denovo, controls$n_denovo, sep ="/")
  data.frame(OR = or, lower = lower, upper = upper, pval = pval, n_denovo_breakdown = case_control_n)
}

generate_background_sets = function() {
  gene_data = load_constraint_data()
  
  gene_data %>%
    filter(!is.na(oe_lof_upper_bin_6)) %>%
    ggplot + aes(x = gene_length, group = oe_lof_upper_bin_6, fill = oe_lof_upper_bin_6) + 
    geom_density(alpha=0.5) + scale_x_log10() + theme_classic()
  
  gene_data %>%
    filter(!is.na(oe_lof_upper_bin_6)) %>%
    mutate(oe_lof_upper_bin_6 = as.factor(oe_lof_upper_bin_6)) %>%
    ggplot + aes(x = gene_length, y = oe_lof_upper_bin_6, fill = oe_lof_upper_bin_6) + 
    geom_density_ridges() + scale_x_log10(labels=comma) + theme_classic()
  
  set.seed(42)
  sample_with_matched_length_distr = function(gene_data) {
    data = gene_data %>% 
      filter(!is.na(oe_lof_upper_bin)) %>%
      mutate(constrained = oe_lof_upper_bin == 0,
             log_length = log10(gene_length))
    
    data %>%
      filter(constrained) %$%
      density(gene_length) -> dens
    data %>%
      filter(!constrained) %$%
      density(gene_length) -> dens_bg
    
    data = data %>%
      mutate(prob = approxfun(dens)(gene_length) / approxfun(dens_bg)(gene_length))
    
    # data %>%
    #   filter(!constrained) %>%
    #   sample_n(sum(data$constrained), replace = F, weight = prob) %>%
    #   union(data %>% filter(constrained)) %T>%
    #   {print(group_by(., constrained) %>% summarize(mean_length=mean(gene_length), n=n()))} %>%
    #   ggplot + aes(x = gene_length, group = constrained, fill = constrained) +
    #   geom_density(alpha=0.2) + theme_classic() + scale_x_log10()
    
    data %>%
      filter(!constrained) %>%
      sample_n(sum(data$constrained), replace = F, weight = prob) %$% gene
  }
  test = replicate(10, sample_with_matched_length_distr(gene_data))
  gene_data %>%
    filter(oe_lof_upper_bin == 0) %>%
    select(gene) %>%
    write.table(file='constrained_genes.tsv', quote=F, row.names=F, col.names = F, sep='\t')
  write.table(test, file='background_genes.tsv', quote=F, row.names=F, col.names = F, sep='\t')
}
  
partitioning_heritability_enrichment = function(save_plot=F, normalize=F) {
  sum_stats = load_sum_stats()
  
  # Broken in tidyr >= 0.8  
  # enrichment_data = sum_stats %>%
  #   filter(cond == 'all100' &
  #            grepl('oe_lof_upper_quantile_', name) &
  #            (is.na(cor_rm0.2) | cor_rm0.2 == 0)  # Select uncorrelated variables
  #   ) %>%
  #   group_by(name) %>%
  #   summarize(meta_enrichment = metagen(enrichment, enrichment_SE)$TE.random,
  #             meta_sd = metagen(enrichment, enrichment_SE)$seTE.random) %>%
  #   mutate(decile = as.integer(str_sub(name, -1)),
  #          enrichment = ifelse(normalize, meta_enrichment / mean(meta_enrichment), meta_enrichment)
  #   )
  
  filt_sum_stats = sum_stats %>%
    filter(cond == 'all100' &
             grepl('oe_lof_upper_quantile_', name) &
             (is.na(cor_rm0.2) | cor_rm0.2 == 0)  # Select uncorrelated variables
    ) 
  enrichment_data = ddply(filt_sum_stats, 'name', function(x) {
    enrichment = x$enrichment
    enrichment_SE = x$enrichment_SE
    meta_enrichment = metagen(enrichment, enrichment_SE)$TE.random
    meta_sd = metagen(enrichment, enrichment_SE)$seTE.random
    return(data.frame(meta_enrichment, meta_sd))
  }) %>%
    mutate(decile = as.integer(str_sub(name, -1)),
           enrichment = ifelse(normalize, meta_enrichment / mean(meta_enrichment), meta_enrichment)
    )
  
  p = enrichment_data %>%
    ggplot + aes(x = decile, y = meta_enrichment,
                 ymin = meta_enrichment - 1.96 * meta_sd, ymax = meta_enrichment + 1.96 * meta_sd) + 
    geom_pointrange() +
    theme_classic() + ylab("Partitioning heritability\nenrichment") + oe_x_axis +
    geom_hline(yintercept = 1, col = "darkgray", linetype = 'dashed')
  
  if (save_plot) {
    pdf('5b_partitioning_heritability_enrich.pdf', height=3, width=4)
    print(p)
    dev.off()
    png('5b_partitioning_heritability_enrich.png', height=3*300, width=4*300, res=300)
    print(p)
    dev.off()
  }
  return(p)
}
  
enriched_traits = function(save_plot=F) {
  sum_stats = load_sum_stats()
  
  continuous_data = sum_stats %>%
    filter(cond == 'all100' & name == 'oe_lof_upper_continuos') %>%
    mutate(logp=-log10(p), Domain=as.factor(Domain))
  
  bonf_threshold = -log10(0.05 / nrow(continuous_data))
  
  highlights = continuous_data %>%
    arrange(desc(logp)) %>%
    # {filter(., logp > bonf_threshold) %>% select(logp, description) %>% print} %>%
    head(5)
  
  p_thres = seq(2, 14, 2)
  set.seed(42)
  p = continuous_data %>%
    mutate(highlight=pheno %in% highlights$pheno,
           Domain=fct_lump(Domain, prop = 0.02),
           color_flag=case_when(
             highlight ~ color_lof,
             as.integer(Domain) %% 2 == 1 ~ "#56B1F7",
             TRUE ~ "#132B43")) %>%
    ggplot + aes(x = Domain, y = logp, col = color_flag) +
    geom_jitter(alpha=0.7, size=2, width=0.5) + 
    geom_text_repel(aes(label=description), data = highlights, color = 'black', size = 3) +
    theme_classic() + xlab(NULL) + 
    ylab('Trait enrichment\np-value') +
    scale_y_continuous(breaks=p_thres, 
                       labels=as.expression(sapply(p_thres, function(x) bquote(10^-.(x))))) +
    scale_color_identity(guide=FALSE) +
    geom_hline(yintercept = bonf_threshold, linetype='dashed', col = 'darkgray') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  continuous_data %>%
    select(-name, -cond) %>% 
    select(pheno, description, Domain, p, enrichment, enrichment_SE, everything()) %>%
    arrange(p) %T>%
    {write.table(., file='data/all_traits.tsv', quote = F, row.names=F, sep='\t')} %>%
    filter(., p < 0.0001) %>%
    select(description, p, enrichment, enrichment_SE, h2_liability) %>%
    mutate(p = format(p, scientific = T, digits=3)) %>%
    mutate_if(is.numeric, round, digits=3) %>%
    write.table(file='data/filtered_traits.tsv', quote = F, row.names=F, sep='\t')
  
   if (save_plot) {
    pdf('5c_enriched_traits.pdf', height=3, width=5)
    print(p)
    dev.off()
    png('5c_enriched_traits.png', height=3*300, width=5*300, res=300)
    print(p)
    dev.off()
  }
  return(p)
}

figure5 = function() {
  p5a = plot_rare_disease(phenotype = 'ddid')
  p5b = partitioning_heritability_enrichment()
  p5c = enriched_traits()
  pdf('figure5.pdf', height=7, width=4.5)
  ggarrange(p5a, p5b, p5c, ncol = 1, nrow = 3, align='v', labels = 'auto', vjust = 1)
  dev.off()
  png('figure5.png', height=7, width=4.5, units = 'in', res=300)
  ggarrange(p5a, p5b, p5c, ncol = 1, nrow = 3, align='v', labels = 'auto', vjust = 1)
  dev.off()
}