source('constants.R')

figure2b = function(save_plot=F) {
  loftee = load_maps_data(cut = 'loftee')
  maps_data = load_maps_data()
  
  syn_maps = maps_data %>% filter(protein_coding & csq == 'synonymous') %$% maps
  mis_maps = maps_data %>% filter(protein_coding & csq == 'missense') %$% maps
  
  loftee_filters = loftee %>%
    filter(!is.na(lof) & variant_count > 20) %>%
    mutate(
      # lof_filter = case_when(lof == 'HC' ~ 'HC',
      #                        lof == 'OS' ~ 'OS',
      #                        grepl(',', lof_filter) ~ 'MULTIPLE',
      #                        TRUE ~ lof_filter),
      lof_filter = case_when(lof == 'HC' ~ 'High Confidence',
                             lof == 'OS' ~ 'Other Splice',
                             grepl(',', lof_filter) ~ 'Multiple filters',
                             lof_filter %in% c('NON_DONOR_DISRUPTING', 'NON_ACCEPTOR_DISRUPTING') ~ 'Non-splice disrupting',
                             lof_filter %in% c('RESCUE_DONOR', 'RESCUE_ACCEPTOR') ~ 'Rescue splice',
                             lof_filter %in% c('3UTR_SPLICE', '5UTR_SPLICE') ~ 'UTR splice',
                             lof_filter == 'END_TRUNC' ~ 'Terminal truncation',
                             TRUE ~ lof_filter),
      lof_filter_simplified = fct_lump(lof_filter, prop=0.0002, w = variant_count)
    ) %>%
    regroup_maps(c('lof', 'lof_filter_simplified')) %>%
    mutate(lof_filter_simplified = fct_reorder(lof_filter_simplified, maps)) %>%
    arrange(maps)
  
  # loftee_colors = rep(color_syn, length(loftee_filters$lof_filter_simplified))
  loftee_colors = rep(color_lc_lof, length(loftee_filters$lof_filter_simplified))
  names(loftee_colors) = loftee_filters$lof_filter_simplified
  loftee_colors['High Confidence'] = color_lof
  loftee_colors['Other Splice'] = color_os
  
  loftee_alphas = rep(0.4, length(loftee_filters$lof_filter_simplified))
  names(loftee_alphas) = loftee_filters$lof_filter_simplified
  loftee_alphas['High Confidence'] = 1
  loftee_alphas['Other Splice'] = 1
  
  ylimits = c(syn_maps, max(loftee_filters$maps_upper))
  
  down_factor = 0.007
  loftee_filters_low_cutoff = loftee_filters %>%
    mutate(maps_lower = ifelse(maps_lower < syn_maps, syn_maps - down_factor, maps_lower))
  
  p2maps_orig = loftee_filters_low_cutoff %>%
    ggplot + theme_classic() +
    aes(x = lof_filter_simplified, y = maps, color = lof_filter_simplified,
        # alpha = lof_filter_simplified,
        ymin = maps_lower, ymax = maps_upper) + 
    geom_pointrange() +
    scale_color_manual(values=loftee_colors, guide='none') +
    # scale_alpha_manual(values=loftee_alphas, guide='none') +
    geom_hline(yintercept = mis_maps, color = color_mis, linetype = 'dashed') + 
    geom_hline(yintercept = syn_maps, color = color_syn, linetype = 'dashed') + 
    annotate('text', x = nrow(loftee_filters_low_cutoff), y = mis_maps + 0.007, size = 3, hjust = 1, color = color_mis, label = 'missense') + 
    annotate('text', x = nrow(loftee_filters_low_cutoff), y = syn_maps + 0.007, size = 3, hjust = 1, color = color_syn, label = 'synonymous') + 
    scale_x_discrete(labels=paste0(loftee_filters$lof_filter_simplified, '\n')) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), 
          plot.margin = margin(0, 5.5, 0, 35.5)) +
    geom_segment(aes(yend = maps_lower, xend = lof_filter_simplified),
                 arrow = arrow(length = unit(0.2, "cm")),
                 data = loftee_filters_low_cutoff %>% filter(maps_lower == syn_maps - down_factor)) +
    xlab(NULL) + ylab('MAPS') + coord_cartesian(ylim=ylimits)
  
  text_loc = -0.042
  if (save_plot) {
    pdf('2a_loftee.pdf', height=3.5, width=6)
    text_loc = -0.03
  }
  p2maps = p2maps_orig
  for (i in 1:nrow(loftee_filters)) {
    p2maps = p2maps + annotation_custom(
      text_grob(paste('n =', comma(loftee_filters$variant_count[i])),
                color = 'darkgray', rot = 30, hjust = 1, vjust = 1, size = 8),
      xmin = i + 1e-1, xmax = i + 1e-1, ymax = 0, ymin = text_loc)
  }
  gt2maps <- ggplot_gtable(ggplot_build(p2maps))
  gt2maps$layout$clip[gt2maps$layout$name == "panel"] <- "off"
  if (save_plot) {
    grid.draw(gt2maps)
    dev.off()
  }
  return(gt2maps)
}
 
figure2a = function(save_plot=F) {
  freq = load_freq_summary_data()
  
  freq_plot = freq %>%
    mutate(freq_bin = tools::toTitleCase(gsub('_', ' ', freq_bin)),
           freq_bin = gsub('Criteria Provided, ', '', freq_bin)) %>%
    filter(!is.na(fail_loftee) & worst_csq %in% lof_like) %>%
    group_by(freq_bin) %>%
    summarize(failed_filters=sum(n * fail_loftee),
              passed_filters=sum(n * pass_loftee),
              n=sum(n), prop_filtered = failed_filters / n) %>%
    mutate(freq_bin=ordered(fct_reorder(freq_bin, prop_filtered)),
           color=ifelse(freq_bin < 'Singleton', color_lof, dataset_colors[['ExAC']]))
  
  # All LOFTEE LoFs (372,970)
  freq_plot %>%
    filter(grepl('0', freq_bin) | grepl('eton', freq_bin)) %>% 
    summarize_if(is.numeric, sum)
  
  # Roughly consistent with canonical transcript:
  # test = load_constraint_data(an_adj=T)
  # sum(test$n_sites, na.rm=T)  # 335987
  # sum(test$n_sites > 0, na.rm=T)  # 16694
  
  bar_width = 0.8
  p = freq_plot %>%
    ggplot + aes(x = freq_bin, y = prop_filtered, fill = color) +
    geom_bar(stat='identity', width=bar_width) + theme_classic() +
    xlab(NULL) + ylab('Percent filtered by LOFTEE') +
    scale_y_continuous(labels=percent_format(accuracy = 1), expand = c(0,0), 
                       limits = c(0, max(freq_plot$prop_filtered) * 1.05)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), 
          plot.margin = margin(0, 5.5, 0, 20.5)) +
    geom_vline(xintercept = 5.5, linetype = 'dashed') + 
    annotate('text', 3, max(freq_plot$prop_filtered), label = 'ClinVar') + 
    annotate('text', 9, max(freq_plot$prop_filtered), label = 'gnomAD') +
    scale_fill_identity()
  # (p = freq_plot %>%
  #     ggplot + aes(x = freq_bin, y = prop_filtered) +
  #     geom_bar(stat='identity', width=bar_width) + theme_classic() +
  #     xlab(NULL) + ylab('Percent filtered by LOFTEE') +
  #     scale_y_continuous(labels=percent_format(accuracy = 1), expand = c(0,0), 
  #                        limits = c(0, max(freq_plot$prop_filtered) * 1.05)) +
  #     theme(axis.text.x = element_text(angle = 30, hjust = 1), 
  #           plot.margin = margin(0, 5.5, 0, 20.5),
  #           axis.line=element_blank()) +
  #     geom_vline(xintercept = 5.5, linetype = 'dashed') + 
  #     annotate('rect', xmin = 0, ymin = -Inf, xmax = 5.5, ymax = Inf, alpha=0.2, fill=color_lof) +
  #     annotate('rect', xmin = Inf, ymin = -Inf, xmax = 5.5, ymax = Inf, alpha=0.2, 
  #              fill=dataset_colors[['ExAC']]) +
  #     geom_bar(stat='identity', width=bar_width) +
  #     annotate('text', 3, max(freq_plot$prop_filtered), label = 'Clinvar') + 
  #     annotate('text', 9, max(freq_plot$prop_filtered), label = 'gnomAD'))
  
  freq %>%
    filter(!is.na(fail_loftee) & worst_csq %in% lof_like) %>%
    group_by(freq_bin) %>%
    summarize(failed_filters=sum(n * fail_loftee), n=sum(n), prop_filtered = failed_filters / n) %>%
    mutate(freq_bin=fct_reorder(freq_bin, prop_filtered)) %>%
    ggplot + aes(x = freq_bin, y = n) +
    geom_bar(stat='identity') + theme_classic() +
    xlab(NULL) + ylab('Total variants') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_y_log10(labels=comma)
  
  if (save_plot) {
    pdf('2b_loftee_frequency.pdf', height=3, width=5)
    print(p)
    dev.off()
  }
  return(p)
}

supp_loftee_flag = function() {
  loftee = load_maps_data(cut = 'loftee')
  maps_data = load_maps_data()
  
  syn_maps = maps_data %>% filter(protein_coding & csq == 'synonymous') %$% maps
  mis_maps = maps_data %>% filter(protein_coding & csq == 'missense') %$% maps
  
  loftee_flags = loftee %>%
    filter(!is.na(lof) & variant_count > 20) %>%
    mutate(lof_flags = case_when(is.na(lof_flags) & lof == 'HC' ~ 'HC',
                                 is.na(lof_flags) & lof == 'OS' ~ 'OS',
                                 is.na(lof_flags) & lof == 'LC' ~ 'LC',
                                 TRUE ~ lof_flags),
           lof_flags_simplified = fct_lump(paste(lof, lof_flags), prop = 0.0002, w = variant_count)
    ) %>%
    regroup_maps(c('lof', 'lof_flags_simplified')) %>%
    mutate(lof_flags_simplified = fct_reorder(lof_flags_simplified, maps)) %>%
    arrange(maps)
  
  loftee_colors = rep(color_syn, length(loftee_flags$lof_flags_simplified))
  names(loftee_colors) = loftee_flags$lof_flags_simplified
  loftee_colors['HC'] = color_lof
  loftee_colors['OS'] = color_os
  loftee_colors['LC'] = color_syn
  
  (p2flags = loftee_flags %>%
      ggplot + theme_classic() +
      # aes(x = lof_filter_simplified, y = maps, ymin = maps_lower, ymax = maps_upper) + geom_pointrange()
      aes(x = lof_flags_simplified, y = maps, color = lof) + geom_point() +
      scale_color_manual(values=loftee_colors, guide='none') + 
      geom_hline(yintercept = mis_maps, color = color_mis, linetype = 'dashed') + 
      geom_hline(yintercept = syn_maps, color = color_syn, linetype = 'dashed') + 
      scale_x_discrete(labels=paste0(loftee_flags$lof_flags_simplified, '\n', loftee_flags$variant_count))
  )
  p2flags + theme(axis.text.x = element_text(angle = 30, hjust = 1)) 
  
  
}
ds_melt = load_downsampled_data()
figure2c = function() {
  sumbounds = ds_melt %>%
    filter(variant_type == 'LoF' & func == 'sum') %>% group_by(obs_exp) %>% 
    summarize(x=max(downsampling), y=max(count))
  
  p2c = ds_melt %>%
    filter(variant_type == 'LoF' & func == 'sum') %>%
    ggplot + aes(x = downsampling, y = count, linetype = obs_exp, color = variant_type) + geom_line(color=color_lof) +
    theme_classic() + scale_linetype_manual(values=linetypes, name=NULL) +
    theme(legend.justification = c(0, 1), legend.position = c(0.01, 1)) +
    scale_x_continuous(label=comma, breaks=c(seq(0, 8e4, 4e4), max(sumbounds$x))) + 
    scale_y_continuous(label=comma) + #, breaks=c(seq(0, 4e5, 1e5), max(sumbounds$y))) +
    # scale_y_continuous(label=comma, breaks=c(0, 1e5, 3e5, 4e5, sumbounds$y)) +
    xlab('Sample size') + ylab('Total number of pLoF SNVs')
  return(p2c)
}

figure2d = function() {
  over10bounds = ds_melt %>%
    filter(variant_type == 'LoF' & func == 'over10') %>% group_by(obs_exp) %>% 
    summarize(x=max(downsampling), y=max(count))
  
  p2d = ds_melt %>%
    filter(variant_type == 'LoF' & func == 'over10') %>%
    ggplot + aes(x = downsampling, y = count, linetype = obs_exp, color = variant_type) + geom_line(color=color_lof) +
    theme_classic() + scale_linetype_manual(values=linetypes, name=NULL) + 
    theme(legend.justification = c(0, 1), legend.position = c(0.01, 1)) +
    scale_x_continuous(label=comma, breaks=c(seq(0, 8e4, 4e4), max(over10bounds$x))) +
    scale_y_continuous(label=percent, breaks=c(0, 0.2, 0.4, 0.6, max(over10bounds$y))) +
    xlab('Sample size') + ylab('Percent of genes with \u2265 10 pLoF SNVs')
  return(p2d)
}

supp_ds = function() {
  ds_melt %>% 
    filter(!is.na(variant_type) & func == 'sum') %>%
    ggplot + aes(x = downsampling, y = count, linetype = obs_exp, color = variant_type) + geom_line() +
    theme_classic() + scale_linetype_manual(values=linetypes, name='') + 
    scale_color_manual(values=colors, name='Variant Type')
  
  ds_melt %>% 
    filter(!is.na(variant_type) & func == 'over10') %>%
    ggplot + aes(x = downsampling, y = count, linetype = obs_exp, color = variant_type) + geom_line() +
    theme_classic() + scale_linetype_manual(values=linetypes, name='') + 
    scale_color_manual(values=colors, name='Variant Type')
  
  ds_melt %>% 
    filter(obs_exp == 'n' & func == 'sum') %>%
    ggplot + aes(x = downsampling, y = count) + geom_line(color=color_lof) +
    theme_classic() # + scale_x_log10() + scale_y_log10()
}

non_can_splice = function() {
  loftee = load_maps_data(cut = 'loftee_simplified_alleles')
  maps_data = load_maps_data()
  
  syn_maps = maps_data %>% filter(protein_coding & csq == 'synonymous') %$% maps
  mis_maps = maps_data %>% filter(protein_coding & csq == 'missense') %$% maps
  
  loftee_filters = loftee %>%
    filter(!is.na(lof) & (worst_csq == 'splice_donor_variant' | worst_csq == 'splice_acceptor_variants') & (is.na(lof_filter_simplified) | lof_filter_simplified == 'NON_CAN_SPLICE')) %>%
    mutate(lof_filter_simplified = case_when(lof == 'HC' ~ 'HC',
                                             lof == 'OS' ~ 'OS',
                                             TRUE ~ lof_filter_simplified)
    ) %>%
    regroup_maps(c('lof', 'worst_csq', 'lof_filter_simplified', 'alleles')) %>%
    arrange(maps)
  
  loftee_colors = rep(color_syn, length(loftee_filters$lof_filter_simplified))
  names(loftee_colors) = loftee_filters$lof_filter_simplified
  loftee_colors['HC'] = color_lof
  loftee_colors['OS'] = color_mis
  
  (p2maps = loftee_filters %>%
      ggplot + theme_classic() +
      # aes(x = lof_filter_simplified, y = maps, ymin = maps_lower, ymax = maps_upper) + geom_pointrange()
      aes(x = lof_filter_simplified, y = maps, color = lof_filter_simplified) + geom_point() +
      scale_color_manual(values=loftee_colors, guide='none') + 
      geom_hline(yintercept = mis_maps, color = color_mis, linetype = 'dashed') + 
      geom_hline(yintercept = syn_maps, color = color_syn, linetype = 'dashed') + 
      scale_x_discrete(labels=paste0(loftee_filters$lof_filter_simplified, '\n', loftee_filters$variant_count))
  )
  p2maps + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
}

figure2 = function() {
  p2a = figure2a() + theme(plot.margin = margin(12.5, 5.5, 5.5, 5.5))
  p2b = figure2b()
  p2c = figure2c() + theme(plot.margin = margin(5.5, 8.5, 5.5, 5.5))
  p2d = figure2d() + theme(plot.margin = margin(5.5, 8.5, 5.5, 5.5))
  pdf('figure2.pdf', height=6, width=9)
  ggarrange(p2a, p2c, p2b, p2d, ncol = 2, nrow = 2, 
            labels = c('a', 'c', 'b', 'd'), align = 'v', vjust = 1, widths=c(1.25, 1))
  dev.off()
  png('figure2.png', height=6, width=9, units = 'in', res=300)
  ggarrange(p2a, p2c, p2b, p2d, ncol = 2, nrow = 2, 
            labels = c('a', 'c', 'b', 'd'), align = 'v', vjust = 1, widths=c(1.25, 1))
  dev.off()
  
  #### Some summary statistics
  function() {
    gene_data %$%
      cor.test(obs_syn, exp_syn)
    
    gene_data %>% 
      summarize(old_over_10 = sum(exac_exp_lof >= 10, na.rm=T), 
                new_over_10 = sum(exp_lof >= 10, na.rm=T),
                old_over_0 = sum(exac_exp_lof >= 0, na.rm=T), 
                new_over_0 = sum(exp_lof >= 0, na.rm=T),
                n=n(),
                old_median = median(exac_exp_lof, na.rm=T),
                new_median = median(exp_lof, na.rm=T),
                prop_old_v2 = old_over_10 / old_over_0,
                prop_new_v2 = new_over_10 / new_over_0
      )
    gene_data %>% 
      summarize(median = median(oe_lof, na.rm=T),
                mean = mean(oe_lof, na.rm=T),
                median_u = median(oe_lof_upper, na.rm=T),
                mean_u = mean(oe_lof_upper, na.rm=T)
      )
  }
}
