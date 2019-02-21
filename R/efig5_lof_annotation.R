gene_data = load_constraint_data()

pre_loftee_data = load_constraint_data(loftee=F) %>%
  transmute(gene, transcript, 
            p_no_loftee=p, classic_caf_no_loftee=classic_caf,
            n_sites_no_loftee=n_sites, max_af_no_loftee=max_af)

gene_data %>%
  left_join(pre_loftee_data) -> loftee_compare_data

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



