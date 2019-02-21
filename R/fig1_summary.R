

summary_figure = function(data_type = 'exomes', group_splice=T, intergenic=F, legend=T,
                   y_label = T, y_axis_position = "left", maps_limits=c(NA, NA), po_limits=c(NA, NA)) {
  obs_poss = load_observed_possible_data(data_type) %>%
    filter(downsampling == max(downsampling))
  
  ordering = c('intergenic', 'intron', "5'UTR", "3'UTR", 'synonymous', 'missense', 
               'essential splice', 'splice donor', 'splice acceptor', 'nonsense')
  maps_data = load_maps_data(data_type=data_type)
  maps_data_plot = maps_data %>%
    filter(variant_count > 100 & (
      protein_coding & worst_csq %in% c(lof_like, 'missense_variant',
                                        'synonymous_variant', 'intron_variant',
                                        '5_prime_UTR_variant', '3_prime_UTR_variant')) | 
        (intergenic & worst_csq == 'intergenic_variant')
      )
    # mutate(csq=fct_reorder(csq, maps))
  
  obs_poss_plot = obs_poss %>%
    inner_join(maps_data_plot, by='csq')
  if (group_splice) {
    maps_data_plot = maps_data_plot %>%
      mutate(csq = fct_recode(csq, 'essential splice' = 'splice donor',
                              'essential splice' = 'splice acceptor')) %>%
      regroup_maps('csq')
    obs_poss_plot = obs_poss_plot %>%
      mutate(csq = fct_recode(csq, 'essential splice' = 'splice donor',
                              'essential splice' = 'splice acceptor')) %>%
      group_by(csq, coverage, variant_type) %>%
      summarize(observed = sum(observed), possible = sum(possible),
                proportion_observed = observed / possible,
                singletons=sum(singletons), 
                po_sem = 1.96 * sqrt(proportion_observed * (1 - proportion_observed) / possible)) %>% ungroup
  }
  maps_data_plot = maps_data_plot %>%
    mutate(csq=fct_relevel(csq, ordering))
  obs_poss_plot = obs_poss_plot %>%
    mutate(csq=fct_relevel(csq, ordering))
  p1maps = maps_data_plot %>%
    ggplot + aes(x = csq, y = maps, ymin = maps_lower, ymax = maps_upper, color = csq) + 
    geom_pointrange() + geom_point(size=2.5) + theme_classic() +
    scale_color_manual(values=variant_category_colors, guide='none') + 
    xlab(NULL) + 
    scale_y_continuous(position=y_axis_position, limits = maps_limits) +
    theme(axis.text.x = element_blank(), 
          plot.margin = margin(0, 5.5, 0, 5.5),
          axis.ticks.x.bottom = element_blank())
  if (y_label) {
    p1maps = p1maps + ylab('MAPS')
  } else {
    p1maps = p1maps + ylab(NULL)
  }
  top_legend = ifelse(is.na(maps_limits[2]), max(maps_data_plot$maps_upper), maps_limits[2])
  p1maps = p1maps + annotate('text', hjust = 0.5, vjust = 1, label = data_type,
                             x = nrow(maps_data_plot) / 2 + 0.5, y = top_legend)
  if (legend) {
    p1maps = p1maps + annotate('text', 0.75, top_legend, hjust=0, vjust=1, label='other', color = color_syn, size = 3) +
    annotate('text', 0.75, 0.9*top_legend, hjust=0, vjust=1, label='missense', color = color_mis, size = 3) +
    annotate('text', 0.75, 0.8*top_legend, hjust=0, vjust=1, label='pLoF', color = color_lof, size = 3)
  }
    
  
  obs_poss_plot %>%
    group_by(coverage) %>%
    summarize(observed = sum(observed), possible = sum(possible)) %>%
    ggplot + aes(x = coverage, y = observed / possible) + geom_line() +
    scale_y_continuous(labels=percent) +
    theme_classic() + xlab('Coverage') + ylab('Percent observed')
  
  p1prop = obs_poss_plot %>%
    filter(coverage >= 30) %>%
    group_by(csq, variant_type) %>%
    summarize(observed = sum(observed), possible = sum(possible),
              proportion_observed = observed / possible,
              singletons=sum(singletons), 
              po_sem = 1.96 * sqrt(proportion_observed * (1 - proportion_observed) / possible)) %>%
    filter(possible > 100) %>%
    ggplot + aes(x = csq, y = proportion_observed, color = variant_type, 
                 ymin = proportion_observed - po_sem, ymax = proportion_observed + po_sem) + 
    geom_pointrange() + geom_point(size=2.5) + xlab(NULL) +
    theme_classic() + scale_color_manual(values=variant_type_colors, guide='none') + 
    scale_y_continuous(labels=percent_format(accuracy = 1), position=y_axis_position, limits=po_limits) + 
    theme(axis.text.x = element_blank(), 
          plot.margin = margin(0, 5.5, 0, 5.5),
          axis.ticks.x.bottom = element_blank())
  
  if (y_label) {
    p1prop = p1prop + ylab('Percent observed')
  } else {
    p1prop = p1prop + ylab(NULL)
  }
  
  if (legend) {
    p1prop = p1prop +
      annotate('text', 0.75, 1, hjust=0, vjust=1, label='CpG transition', color = color_cpg, size = 3) +
      annotate('text', 0.75, 0.92, hjust=0, vjust=1, label='non-CpG transition', color = color_ti, size = 3) +
      annotate('text', 0.75, 0.84, hjust=0, vjust=1, label='transversion', color = color_tv, size = 3)
  }
  
  obs_poss_plot %>%
    filter(coverage >= 30) %>%
    group_by(csq) %>%
    summarize(observed = sum(observed), possible = sum(possible)) %>%
    ggplot + aes(x = csq, y = observed / possible, color = csq) + geom_point() + 
    theme_classic() + scale_color_manual(values=variant_category_colors, guide='none') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab(NULL) +
    scale_y_continuous(labels=percent) + ylab('Percent observed')
  
  p1hist_orig = obs_poss_plot %>%
    inner_join(maps_data_plot, by='csq') %>%
    filter(coverage >= 30 & csq %in% maps_data_plot$csq) %>%
    group_by(csq) %>%
    summarize(observed = sum(observed), possible = sum(possible)) %>%
    ggplot + aes(x = csq, y = observed, fill = csq) + geom_bar(stat='identity') +
    theme_classic() + scale_fill_manual(values=variant_category_colors, guide='none') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab(NULL) +
    scale_y_continuous(labels=comma)
  
  trans = function(x) pmin(x, breakpoint) + scaling_factor * pmax(x - breakpoint, 0)
  
  if (data_type == 'genomes') {
    scaling_factor = 0.01
    breakpoint = 1.5e6
    yticks = c(0, 2.5e5, 5e5, 7.5e5, 1e6, 1.25e6, 1.5e6, 2e7, 4e7, 6e7)
  } else {
    scaling_factor = 0.05
    breakpoint = 2.2e5
    yticks = c(0, 5e4, 1e5, 1.5e5, 2e5, 1e6, 2e6, 3e6, 4e6, 5e6)
  }
  
  break_end = 1.1
  p1hist = obs_poss_plot %>%
    filter(coverage >= 30) %>%
    group_by(csq) %>%
    summarize(observed = sum(observed), observed_t = trans(observed)) %>%
    ggplot + aes(x = csq, y = observed_t, fill = csq) + geom_bar(stat='identity') + 
    theme_classic() + scale_fill_manual(values=variant_category_colors, guide='none') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) + xlab(NULL) +
    # geom_rect(aes(xmin=0.4, xmax=Inf, ymin=breakpoint, ymax=breakpoint * 1.1), fill='white', 
    # color='gray', linetype='dashed') +
    geom_rect(aes(xmin=0.4, xmax=Inf, ymin=breakpoint, ymax=breakpoint * break_end), fill='white') +
    geom_hline(yintercept = breakpoint, color='gray', linetype='dashed') +
    geom_hline(yintercept = breakpoint * break_end, color='gray', linetype='dashed') +
    scale_y_continuous(limits=c(0, NA), breaks=trans(yticks), labels=comma(yticks),
                       expand = c(0, 0), position=y_axis_position) +
    theme(axis.line.y = element_blank()) 
  
  if (y_label) {
    p1hist = p1hist + ylab('Total observed')
  } else {
    p1hist = p1hist + ylab(NULL)
  }
  
  if (legend) {
    p1hist = p1hist + 
      annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = breakpoint) +
      annotate(geom = "segment", x = -Inf, xend = -Inf, y = breakpoint * 1.1, yend = Inf) +
      annotate(geom = "segment", x = -Inf, xend = -Inf, y = breakpoint, yend = breakpoint * 1.1,
               color = "white")
  }
  # (p1maps / p1prop) / p1hist
  return(list(p1maps, p1prop, p1hist))
}

figure1 = function() {
  p1e = summary_figure('exomes', group_splice = T, intergenic = F,
                     maps_limits=c(NA, 0.18), po_limits=c(NA, 1))
  p1g = summary_figure('genomes', group_splice = T, intergenic = T, legend = F, y_label = F, 
                     y_axis_position = 'right', maps_limits=c(NA, 0.18), po_limits=c(NA, 1))
  
  # g2 = summary_figure('genomes', group_splice = T, intergenic = T)
  # pdf('1_exomes.pdf', height=6, width=4)
  # ggarrange(plotlist = e, nrow = 3, ncol = 1, align = 'v', labels = 'auto', vjust = 1)
  # dev.off()
  # pdf('1_genomes.pdf', height=6, width=4)
  # ggarrange(plotlist = g2, nrow = 3, ncol = 1, align = 'v', labels = 'auto', vjust = 1)
  # dev.off()
  
  pdf('figure1.pdf', height=7, width=6, onefile = F)
  egg::ggarrange(plots = list(p1e, p1g) %>% transpose() %>% do.call(base::c, .),
                 nrow = 3, ncol = 2, align = 'v',
                 labels=c('c', 'd', 'e', 'f', 'g', 'h'),
                 label.args=list(gp=grid::gpar(font=2, cex=1.2)), vjust = 1)
  dev.off()
  png('figure1.png', height=7, width=6, units = 'in', res=300)
  egg::ggarrange(plots = list(p1e, p1g) %>% transpose() %>% do.call(base::c, .),
                 nrow = 3, ncol = 2, align = 'v',
                 labels=c('c', 'd', 'e', 'f', 'g', 'h'),
                 label.args=list(gp=grid::gpar(font=2, cex=1.2)), vjust = 1)
  dev.off()
}
