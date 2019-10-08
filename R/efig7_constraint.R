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
    pdf(paste0('e7_oe_distribution_', metric, '.pdf'), height=2, width=3)
    print(p)
    dev.off()
  }
  return(p)
}

oe_v_loeuf = function(save_plot=F, x_metric='oe_lof_upper', y_metric='oe_lof') {
  x_label = get_metric_common_name(x_metric)
  y_label = get_metric_common_name(y_metric)
  
  p = gene_data %>%
    filter(oe_lof < 1.5 & exp_lof > 0) %>%
    mutate(observed=if_else(obs_lof > 10, 10, obs_lof)) %>%
    ggplot + aes_string(x = x_metric, y = y_metric) + aes(color = observed) +
    geom_point_rast() + scale_color_gradient(breaks=seq(0, 10, 2)) +
    # scale_color_gradientn(colors = gradient_colors, values = rescale(gradient_values)) +
    # scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values)) +
    theme_classic() + xlab(x_label) + ylab(y_label)
  
  if (save_plot) {
    pdf(paste0('e7_', x_metric, '_v_', y_metric, '.pdf'), height=2, width=3)
    print(p)
    dev.off()
  }
  return(p)
}

oe_v_loeuf_lines = function(save_plot=F) {
  p = gene_data %>%
    filter(oe_lof < 1.5 & exp_lof > 0) %>%
    # sample_frac(0.1) %>%
    mutate(idx=dense_rank(-oe_lof)) %>%
    ggplot + aes(x = idx, y = oe_lof, ymin = oe_lof_lower, ymax = oe_lof_upper) +
    geom_pointrange(alpha=0.01) + coord_flip() + 
    # scale_color_gradientn(colors = gradient_colors, values = rescale(gradient_values)) +
    # scale_fill_gradientn(colors = gradient_colors, values = rescale(gradient_values)) +
    theme_classic() + ylab('Observed/Expected') + xlab('O/E Rank')
  
  if (save_plot) {
    pdf(paste0('e7_oe_cis.pdf'), height=2, width=3)
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
  print(ggarrange(e7row1, egg::set_panel_size(e7c, width=unit(2, 'in')), e7row34, align='hv',
            nrow = 3, ncol = 1,
            labels = c('', 'c', ''),
            heights = c(1, 1, 2)
  ))
  dev.off()
  png('extended_data_figure7.png', height=8, width=7.5, units = 'in', res=300)
  # egg::ggarrange(plots=list(a, b, c, d, e, f, g),
  #                nrow = 4, ncol = 2, align = 'v', 
  #                labels = c('a', 'b', 'c', 'd', 'e', 'f', 'g'), 
  #                label.args=list(gp=grid::gpar(font=2, cex=1.2)), vjust = 1)
  print(ggarrange(e7row1, egg::set_panel_size(e7c, width=unit(2, 'in')), e7row34, align='hv',
            nrow = 3, ncol = 1,
            labels = c('', 'c', ''),
            heights = c(1, 1, 2)
  ))
  dev.off()
 }

sfigure6 = function() {
  s6a = oe_v_loeuf()
  s6b = oe_v_loeuf(y_metric='pLI')
  s6c = oe_v_loeuf(x_metric='oe_lof', y_metric='pLI')
  s6d = oe_v_loeuf_lines()
  pdf('supplementary_figure6.pdf', height=6, width=8)
  print(ggarrange(s6a, s6b, s6c, s6d, ncol = 2, nrow = 2, labels = 'auto', align = 'v'))
  dev.off()
  png('supplementary_figure6.png', height=6, width=8, units = 'in', res=300)
  print(ggarrange(s6a, s6b, s6c, s6d, ncol = 2, nrow = 2, labels = 'auto', align = 'v'))
  dev.off()
}

get_expected_projection_data = function(force_recompute=F) {
  fname = 'supplementary_dataset_12_n_required.tsv.gz'
  if (!force_recompute) {
    fname = get_or_download_file(fname, subfolder = 'misc_files/')
  } else {
    ds_gene = load_downsampled_gene_data()
    gene_data = ds_gene %>%
      filter(canonical & pop == 'global' & downsampling >= 100) %>%
      group_by(gene) %>%
      filter(min(exp_lof) > 0 & min(exp_mis) > 0 & min(exp_syn) > 0) %>%
      mutate(log_exp_lof=log10(exp_lof), 
             log_exp_mis=log10(exp_mis),
             log_exp_syn=log10(exp_syn),
             log_n=log10(downsampling)) %>%
      do(lof_fit = lm(log_exp_lof ~ log_n, data = .),
         mis_fit = lm(log_exp_mis ~ log_n, data = .),
         syn_fit = lm(log_exp_syn ~ log_n, data = .))
    
    n_lof = 10
    gene_lof_fit = gene_data %>% tidy(lof_fit) %>% 
      summarize(slope=sum(estimate*(term == 'log_n')), intercept=sum(estimate*(term == '(Intercept)')))
    gene_mis_fit = gene_data %>% tidy(mis_fit) %>%
      summarize(slope=sum(estimate*(term == 'log_n')), intercept=sum(estimate*(term == '(Intercept)')))
    
    post_process_predictions <- function(data) {
      data %>%
        mutate(`5` = 10 ^ ((log10(5) - intercept) / slope),
               `10` = 10 ^ ((log10(10) - intercept) / slope),
               `20` = 10 ^ ((log10(20) - intercept) / slope),
               `50` = 10 ^ ((log10(50) - intercept) / slope),
               `100` = 10 ^ ((log10(100) - intercept) / slope)) %>%
        select(-slope, -intercept) %>%
        gather('n_variants', 'n_required', -gene) %>%
        mutate(n_variants = as.numeric(n_variants)) %>%
        group_by(n_variants) %>%
        mutate(rank=percent_rank(n_required)) %>% ungroup
    }
    samples_required_lof = post_process_predictions(gene_lof_fit)
    samples_required_mis = post_process_predictions(gene_mis_fit)
    output_file = gzfile(paste0('data/', fname), 'w')
    samples_required_lof %>% mutate(variant_type = 'pLoF') %>%
      union_all(samples_required_mis %>% mutate(variant_type = 'missense')) %>%
      write.table(output_file, quote=F, row.names = F, sep = '\t')
    close(output_file)
  }
  read_delim(gzfile(fname), delim = '\t')
}

expected_projections = function(projection_df, label='pLoF') {
  # samples_required_lof %>%
  #   filter(n_variants == 10) %>%
  #   ggplot + aes(y = rank, x = n_required) + geom_line(color=color_lof, size=2) +
  #   theme_classic() + scale_x_log10(label=comma) + 
  #   scale_y_continuous(label=percent) + 
  #   xlab('Sample size required') + ylab('Percent of human genes') +
  #   geom_vline(xintercept = 60706, linetype='dotted') +
  #   geom_vline(xintercept = 141456, linetype='dashed')
  projection_df %>%
    mutate(n_variants=fct_reorder(as.factor(n_variants), n_variants)) %>%
    filter(n_variants == 10 & rank <= 0.95) %>%
    arrange(desc(rank)) %>% head
  
  xlimits=c(100, 1e8)
  p = projection_df %>%
    mutate(n_variants=fct_reorder(as.factor(n_variants), n_variants)) %>%
    ggplot + aes(y = rank, x = n_required, color = n_variants) + 
    geom_line(size=2) + theme_classic() + scale_x_log10(label=comma, limits=xlimits) + 
    scale_y_continuous(label=percent) +
    scale_color_discrete(name='>= N variants\nexpected', h = c(40, 120)) +
    # scale_color_manual(name='>= N variants\nexpected', values=c('5' = 'blue', '10' = 'green')) +
    xlab('Sample size required') + ylab('Percent of human genes') +
    geom_vline(xintercept = 60706, linetype='dotted') +
    geom_vline(xintercept = 141456, linetype='dashed') + 
    annotate('text', x = xlimits[1], y = 1, hjust = 0, vjust = 1, label = label)
  return(p)
}

sfigure7 = function() {
  projection_data = get_expected_projection_data()
  s7a = expected_projections(projection_data %>% filter(variant_type == 'pLoF'))
  s7b = expected_projections(projection_data %>% filter(variant_type == 'missense'), 'Missense')
  pdf('supplementary_figure7.pdf', height=6, width=6)
  print(ggarrange(s7a, s7b, ncol = 1, nrow = 2, labels = 'auto', align = 'v'))
  dev.off()
  png('supplementary_figure7.png', height=6, width=6, units = 'in', res=300)
  print(ggarrange(s7a, s7b, ncol = 1, nrow = 2, labels = 'auto', align = 'v'))
  dev.off()
}

sfigure8 = function() {
  or_genes = gene_data %>%
    filter(grepl('^OR', gene)) %>%
    transmute(gene = gene, gene_list = 'Olfactory Genes')
  
  gene_lists = load_all_gene_list_data() %>%
    bind_rows(or_genes)
  cloud_data = gene_data %>%
    left_join(gene_lists) %>%
    filter(!is.na(gene_list) & !is.na(oe_lof_upper_bin)) %>%
    mutate(gene_list = if_else(grepl('haploinsufficiency', gene_list), 'Haploinsufficient', gene_list),
           gene_list = fct_recode(gene_list, 'Autosomal Recessive' = "all_ar", 'Autosomal Dominant' = "all_ad"),
           gene_list = fct_relevel(gene_list, 'Haploinsufficient'),
           oe_lof_upper_bin = label_function2(oe_lof_upper_bin)) %>%
    filter(gene_list %in% c('Haploinsufficient', 'Autosomal Dominant', 'Autosomal Recessive', 'Olfactory Genes')) %>%
    select(gene, gene_list, oe_lof_upper_bin) 
  
  # png('wordcloud.png', height=10, width=6.5, res=600)
  # print(p)
  # dev.off()
  # 
  # pdf('wordcloud_filtered.pdf', height=10, width=6.5)
  # p = cloud_data %>%
  #   ggplot + aes(label = gene) +
  #   geom_text_wordcloud_area(rm_outside = TRUE) +
  #   scale_size_area(max_size = 2) +
  #   facet_grid(oe_lof_upper_bin ~ gene_list)
  # print(p)
  # dev.off()
  
  f = gzfile('data/supplementary_dataset_13_gene_lists.tsv.gz', 'w')
  cloud_data %>% 
    write.table(f, row.names = F, quote=F, sep = '\t')
  close(f)
  
  p = cloud_data %>% 
    add_count(oe_lof_upper_bin, gene_list) %>%
    mutate(n = n + 0.1) %>%
    ggplot + aes(label = gene, size = 100 / log10(2 * n), color = gene_list) +
    geom_text_wordcloud_area(rm_outside = TRUE, grid_size = 0, max_grid_size = 0, grid_margin = 0, seed = 42) +
    scale_size_area(max_size = 29) + theme_minimal() +
    scale_color_manual(values=gene_list_colors, guide=F) +
    facet_grid(oe_lof_upper_bin ~ gene_list) +
    theme(strip.text = element_text(size = 48))
  pdf('supplementary_figure8.pdf', height=100, width=65)
  print(p)
  dev.off()
  res = 72
  png('supplementary_figure8.png', height=100*res, width=65*res, res=res)
  print(p)
  dev.off()
  
  # pdf('wordcloud_full.pdf', height=30, width=19.5)
  # p = cloud_data %>%
  #   ggplot + aes(label = gene) +
  #   geom_text_wordcloud_area() +
  #   scale_size_area(max_size = 2) + 
  #   facet_grid(oe_lof_upper_bin ~ gene_list)
  # print(p)
  # dev.off()
}

load_rvis_comparison_data = function() {
  or_genes = gene_data %>%
    filter(grepl('^OR', gene)) %>%
    transmute(gene = gene, gene_list = 'Olfactory Genes')
  
  ko_gene_lists = get_ko_gene_lists()
  gene_lists = load_all_gene_list_data() %>%
    bind_rows(or_genes) %>%
    bind_rows(ko_gene_lists)
  
  fname = paste0(data_dir, 'RVIS_Unpublished_ExACv2_March2017.txt')
  if (!file.exists(fname)) {
    download.file('http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt', fname)
  }
  rvis_data = read_delim(fname, delim='\t',
                         col_types = cols(
                           `OE-ratio_[ExAC v2]` = col_double(),
                           `%OE-ratio_[ExAC v2]` = col_double()))
  
  # s_het = read_xlsx('../misc_files/ng.3831-S2.xlsx')
  
  joined_data = gene_data %>%
    inner_join(rvis_data, by=c('gene' = 'CCDSr20')) #%>%
    # left_join(s_het, by=c('gene' = 'gene_symbol')) %>%
    # mutate(NFE_s_het = as.numeric(NFE_s_het))
  
  de_novo_data = get_de_novo_data() %>%
    filter(group %in% c('ddid', 'Control')) %>%
    mutate(csq = fct_relevel(case_when(csq %in% lof_like & lof == 'HC' ~ 'pLoF',
                                       csq == 'missense_variant' ~ 'Missense',
                                       csq == 'synonymous_variant' ~ 'Synonymous',
                                       TRUE ~ format_vep_category(csq)),
                             'pLoF', 'Missense', 'Synonymous')) %>%
    left_join(gene_data, by=c('ensg' = 'gene_id')) %>%
    filter(exp_lof >= 10) %>%
    count(csq, group, symbol, oe_lof_upper_bin) %>%
    spread('group', 'n', fill = 0) %>%
    mutate(dz_gene=Control == 0 & ddid > 1,
           control_gene=Control > 1 & ddid == 0) %>% 
    filter(!(dz_gene & control_gene))
  
  hi_data = joined_data %>%
    mutate(is_hi = gene %in% (gene_lists %>% filter(grepl('haploinsufficiency', gene_list)) %$% gene),
           is_mouse_het_lethal = gene %in% (gene_lists %>% filter(gene_list == 'mouse_het_lethal_genes') %$% gene),
           is_cell_essential = gene %in% (gene_lists %>% filter(gene_list == 'CEGv2') %$% gene),
           is_dd_gene = gene %in% (de_novo_data %>% filter(dz_gene) %$% symbol),
           is_dd_control_gene = gene %in% (de_novo_data %>% filter(control_gene) %$% symbol),
           is_essential = is_mouse_het_lethal | is_cell_essential,
           is_hi_or_essential = is_hi | is_essential)
  
  hi_data %>%
    filter(!is.na(`RVIS[pop_maf_0.05%(any)]`) & !is.na(cds_length)& !is.na(oe_lof_upper)) %>%
    select(oe_lof_upper, `RVIS[pop_maf_0.05%(any)]`, cds_length) %>%
    distinct() %>%
    cor
  return(hi_data)
}

rvis_compare_data = load_rvis_comparison_data()
rvis_comparisons = function(outcome_var='is_hi', add_pLI=F, add_exac_pLI=F, add_exac_s_het=F, 
                            drop_rvis=F, compute_partial_aucs=F,
                            save_plot=F) {
  rvis_compare_data$outcome = rvis_compare_data[[outcome_var]]
  if (outcome_var == 'is_dd_gene') {
    rvis_compare_data = rvis_compare_data %>% filter(is_dd_gene | is_dd_control_gene)
  }
  glm(outcome ~ 
        oe_lof_upper +
        pLI +
        `RVIS[pop_maf_0.05%(any)]` +
        # `OE-ratio_[ExAC v2]` +
        cds_length,
      data=rvis_compare_data, family='binomial') %>%
    summary
  
  if (compute_partial_aucs) {
    rvis_compare_data %>%
      mutate(pLI = -pLI,
             exac_pLI = -exac_pLI,
             RVIS = `RVIS[pop_maf_0.05%(any)]`,
             LOEUF = oe_lof_upper) -> formatted_comparison_data
  
    multiple_rocs = ldply(c(99:71/100, 7:1/10), function(x) {
      formatted_comparison_data %>%
        roc(outcome, LOEUF,
            partial.auc=c(1, x)
            ) -> loeuf_roc_run
      rvis_compare_data %>%
        roc(outcome, pLI,
            partial.auc=c(1, x)
            ) -> pLI_roc_run
      res = roc.test(loeuf_roc, pLI_roc)
      data.frame(spec=x, loeuf_roc_run=res$estimate[[1]], pli_roc_run=res$estimate[[2]], p=res$p.value)
    })
    p = multiple_rocs %>%
      mutate(significant = p < 0.05) %>%
      ggplot + aes(x = spec, shape = significant, y = loeuf_roc / pli_roc) +
      geom_point() + geom_line() + xlab('Specificity') + ylab('Ratio of AUCs')
  } else {
    plot_vars = c('outcome', 'LOEUF')
    if (add_pLI) plot_vars = c(plot_vars, 'pLI')
    if (add_exac_pLI) plot_vars = c(plot_vars, 'exac_pLI')
    if (add_exac_s_het) plot_vars = c(plot_vars, 'NFE_s_het')
    if (!drop_rvis) plot_vars = c(plot_vars, 'RVIS')
    metric_colors = c('LOEUF' = color_lof, 'RVIS' = '#F0810F', 'pLI' = 'darkblue', 'exac_pLI' = 'lightblue', 'NFE_s_het' = 'pink')
    p = rvis_compare_data %>%
      mutate(pLI = -pLI,
             exac_pLI = -exac_pLI,
             # NFE_s_het = -NFE_s_het,
             RVIS = `RVIS[pop_maf_0.05%(any)]`,
             LOEUF = oe_lof_upper) %>%
      select_at(vars(plot_vars)) %>%
      gather('Metric', 'score', -outcome) %>%
      mutate(Metric = fct_relevel(Metric, plot_vars[0:-1])) %>%
      ggplot + aes(m = -score, d = outcome, color = Metric) + 
      geom_roc(n.cuts = 0) + xlab('False Positive Fraction') +
      ylab('True Positive Fraction') + theme_classic() +
      scale_color_manual(values=metric_colors, guide=F)

    aucs = calc_auc(p)$AUC
    names(aucs) = plot_vars[0:-1]
    print(aucs)
    n_digits = 4
    y_locs = rev(plot_vars[0:-1])
    p = p + annotate('text', 1, 0.1 * which(y_locs == 'LOEUF') - 0.1, hjust=1, vjust=0, label=paste0('LOEUF (AUC: ', format(aucs['LOEUF'], digits = n_digits), ')'), 
               color=metric_colors['LOEUF'])
    
    if (!drop_rvis) {
      p = p + annotate('text', 1, 0, hjust=1, vjust=0, label=paste0('RVIS (AUC: ' , format(aucs['RVIS'], digits = n_digits), ')'), 
               color=metric_colors['RVIS'])
    }
    if (add_pLI) {
      p = p + 
        annotate('text', 1, 0.1 * which(y_locs == 'pLI') - 0.1, hjust=1, vjust=0, label=paste0('pLI (AUC: ', format(aucs['pLI'], digits = n_digits), ')'), 
                 color=metric_colors['pLI'])
    }
    if (add_exac_s_het) {
      p = p + 
        annotate('text', 1, 0.1 * which(y_locs == 'NFE_s_het') - 0.1, hjust=1, vjust=0, label=paste0('s het (AUC: ', format(aucs['NFE_s_het'], digits = n_digits), ')'), 
                 color=metric_colors['NFE_s_het'])
    }
    if (add_exac_pLI) {
      p = p + 
        annotate('text', 1, 0.1 * which(y_locs == 'exac_pLI') - 0.1, hjust=1, vjust=0, label=paste0('ExAC pLI (AUC: ', format(aucs['exac_pLI'], digits = n_digits), ')'), 
                 color=metric_colors['exac_pLI'])
    }
    
  }

  if (save_plot) {
    pdf('rvis_compare.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

synonymous_outliers = function() {
  mappability_coverage_summary_fname = get_or_download_file('gencode_grch37_gene_by_platform_coverage_summary.tsv.gz', subfolder = 'summary_gene_coverage/', version = 'v1.1')
  mappability_coverage_summary = read_delim(mappability_coverage_summary_fname, delim='\t')
  mappability_coverage_summary %<>% 
    filter(platform == 'all' & data_type == 'exomes') %>%
    separate(gene_id, c('gene_id', 'version')) %>%
    select(-gene_interval, -data_type, -platform, -version)
  
  fname = paste0(data_dir, 'gencode.v19.2wayconspseudos.gtf.gz')
  if (!file.exists(fname)) {
    download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.2wayconspseudos.gtf.gz', fname)
  }
  pseudogene_data = read_delim(fname, delim='\t', skip=5, col_names=F) %>%
    separate(X9, c('other', 'parent'), sep='parent_id "') %>%
    select(parent) %>%
    mutate(parent=gsub('";', '', parent, fixed = T)) %>% 
    distinct
  
  gene_data %>%
    left_join(mappability_coverage_summary) %>%
    mutate(pseudogene = gene_id %in% pseudogene_data$parent,
           syn_outlier = syn_z < -3.71) %>%
    filter(!is.na(mean_mappability) & !is.na(syn_outlier)) %T>%
    # {count(., pseudogene, syn_outlier) %>% print} %>%
    {count(., mean_mappability < 0.9, syn_outlier) %>% print} %>%
    ggplot + aes(x = mean_mappability) + geom_histogram(alpha=0.5) + 
    facet_grid(syn_outlier ~ ., scales = 'free')
  
  gene_data %>%
    left_join(mappability_coverage_summary) %>%
    mutate(syn_outlier = syn_z < -3.71) %>%
    filter(!is.na(mean_mappability) & !is.na(syn_outlier)) %>%
    ggplot + aes(x = mean_mappability) + geom_histogram(alpha=0.5) + 
    facet_grid(syn_outlier ~ ., scales = 'free')
}

sfigure9 = function() {
  s9a = rvis_comparisons('is_hi', add_pLI=T)
  s9b = rvis_comparisons('is_essential', add_pLI=T)
  pdf('supplementary_figure9.pdf', height=4, width=8)
  print(ggarrange(s9a, s9b, ncol = 2, nrow = 1, labels = 'auto', align = 'v'))
  dev.off()
  res = 300
  png('supplementary_figure9.png', height=4*res, width=8*res, res=res)
  print(ggarrange(s9a, s9b, ncol = 2, nrow = 1, labels = 'auto', align = 'v'))
  dev.off()
}

constraint_by_downsampling = function() {
  ds_data = load_downsampled_gene_data() %>%
    mutate(oe_lof_upper_bin = ntile(oe_lof_upper, 10)) %>%
    filter(canonical)

  or_genes = gene_data %>%
    filter(grepl('^OR', gene)) %>%
    transmute(gene = gene, gene_list = 'Olfactory Genes')
  
  gene_lists = load_all_gene_list_data() %>%
    bind_rows(or_genes) %>%
    bind_rows(get_ko_gene_lists())
  
  constraint_by_ds = ds_data %>%
    mutate(is_hi = gene %in% (gene_lists %>% filter(grepl('haploinsufficiency', gene_list)) %$% gene),
           is_mouse_het_lethal = gene %in% (gene_lists %>% filter(gene_list == 'mouse_het_lethal_genes') %$% gene),
           is_cell_essential = gene %in% (gene_lists %>% filter(gene_list == 'CEGv2') %$% gene),
           is_essential = is_mouse_het_lethal | is_cell_essential,
           is_hi_or_essential = is_hi | is_essential)
  return(constraint_by_ds)
}

constraint_by_ds = constraint_by_downsampling()
plot_constraint_by_downsampling = function(outcome_var='is_hi', ylabel='AUC (Haploinsufficient genes)', 
                                           filter_to_rvis_comparisons=F, save_plot=F) {
  constraint_by_ds$outcome = constraint_by_ds[[outcome_var]]
  
  ds = constraint_by_ds %>%
    filter(downsampling >= 1000 & pop == 'global')
  
  if (filter_to_rvis_comparisons) {
    ds = ds %>% semi_join(rvis_compare_data %>% select(gene, transcript))
  }
  p = ds %>%
    mutate(downsampling = as.factor(downsampling)) %>%
    mutate(LOEUF = -oe_lof_upper) %>%
    ggplot + aes(m = LOEUF, d = outcome, color = downsampling) + 
    geom_roc(n.cuts = 0) + xlab('False Positive Fraction') +
    ylab('True Positive Fraction')
  aucs = calc_auc(p)$AUC
  print(aucs)
  
  p = ds %>%
    count(downsampling) %>%
    cbind(aucs) %>%
    ggplot + aes(x = downsampling, y = aucs) + 
    geom_line(size=1) + xlab('Sample size') + 
    theme(plot.margin = margin(t = 5.5, r = 12.5, b = 5.5, l = 5.5)) + 
    scale_x_continuous(breaks=ds_breaks) + ylab(ylabel) + theme_classic()
  
  if (save_plot) {
    pdf('constraint_by_downsampling.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

sfigure10 = function() {
  s10a = plot_constraint_by_downsampling('is_hi', 'AUC (Haploinsufficient genes)')
  s10b = plot_constraint_by_downsampling('is_essential', 'AUC (Essential genes)')
  
  pdf('supplementary_figure10.pdf', height=4, width=8)
  print(ggarrange(s10a, s10b, ncol = 2, nrow = 1, labels = 'auto', align = 'v'))
  dev.off()
  res = 300
  png('supplementary_figure10.png', height=4*res, width=8*res, res=res)
  print(ggarrange(s10a, s10b, ncol = 2, nrow = 1, labels = 'auto', align = 'v'))
  dev.off()
}


