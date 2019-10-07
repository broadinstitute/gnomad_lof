
# Populations and percent ancestries
all_data = bind_rows(
  data.frame(dataset='1000\nGenomes', eas=523, sas=494, eur=514, afr=691, amr=355),
  data.frame(dataset='ESP', eur=4298, afr=2217),
  data.frame(dataset='ExAC', eas=4327, sas=8256, eur=33370+3307, afr=5203, amr=5789, oth=454),
  data.frame(dataset='gnomAD', eas=9197, sas=15308, eur=56885+10824, asj=5040, afr=8128, amr=17296, oth=3070),
  data.frame(dataset='gnomAD', eas=780, eur=7718+1738, asj=145, afr=4359, amr=424, oth=513+31),
  data.frame(dataset='discovEHR', unk=50000),
  # data.frame(dataset='Open\nSearch', oth=10545),
  # data.frame(dataset='GME', mde=2497),
  data.frame(dataset='BRAVO', unk=62784)
) %>% 
  gather(pop, count, -dataset) %>% filter(!is.na(count)) %>%
  group_by(dataset, pop) %>% summarize(count = sum(count)) %>% ungroup %>%
  mutate(dataset = fct_reorder(dataset, count, .fun = sum),
         pop = fct_relevel(pop, 'oth', 'unk', after=Inf))

subpop_data = tribble(
  ~dataset, ~pop, ~subpop, ~count,
  'genomes', 'eur', 'est', 2297,
  'genomes', 'eur', 'bgr', 6,
  'genomes', 'eur', 'seu', 53,
  'genomes', 'eur', 'onf', 1062,
  'genomes', 'eur', 'nwe', 4299,
  'genomes', 'eur', 'swe', 1,
  'genomes', 'eur', 'fin', 10824,
  'exomes', 'eur', 'est', 121,
  'exomes', 'eur', 'bgr', 1335,
  'exomes', 'eur', 'seu', 5752,
  'exomes', 'eur', 'onf', 15499,
  'exomes', 'eur', 'nwe', 21111,
  'exomes', 'eur', 'swe', 13067,
  'exomes', 'eur', 'fin', 3307,
  'exomes', 'eas', 'oea', 7212,
  'exomes', 'eas', 'kor', 1909,
  'exomes', 'eas', 'jpn', 76
) %>% group_by(pop, subpop) %>% summarize(count = sum(count)) %>% 
  ungroup %>% mutate(subpop = fct_relevel(subpop, 'onf', 'oea', after=Inf))

plot_panel_size = function(save_plot=F, return=F) {
  pop_names['onf'] = 'Other European'
  pop_names['unk'] = 'Unknown/Not Shared'
  pop_colors_border = pop_colors
  pop_colors_border['oth'] = 'gray'
  pop_colors['oth'] = 'white'
  angle = 0
  if (save_plot) {
    png('gnomad_barplot_compare.png',height=1600,width=2350,res=300)
  }
  full_full_p = all_data %>%
    mutate(dataset=fct_reorder(dataset, count, sum)) %>%
    ggplot + aes(x=dataset, y=count, fill=pop, color=pop) +
    geom_bar(stat='identity') + theme_classic() +
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_color_manual(name='', values=pop_colors_border, labels=pop_names) +
    scale_y_continuous(breaks=c(seq(0, 130000, 10000), 141456)) + ylab('') + xlab('') +
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_text(margin = ggplot2::margin(t = -10), size = 15, angle=angle),
          legend.position = c(0.025, 1), legend.justification = c('left', 'top'),
          legend.key.size = unit(0.5, 'cm')
    )
  if (return) {
    return(full_full_p)
  } else {
    print(full_full_p)
  }
  if (save_plot) dev.off()
  if (save_plot) {
    png('gnomad_barplot.png',height=1600,width=2350,res=300)
  }
  full_p = all_data %>%
    filter(!(dataset %in% c('discovEHR', 'BRAVO'))) %>%
    ggplot + aes(x=dataset, y=count, fill=pop) +
    geom_bar(stat='identity') + theme_classic() +
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_y_continuous(breaks=c(seq(0, 130000, 10000), 141456)) + ylab('') + xlab('') +
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_text(margin = ggplot2::margin(t = -10), size = 15, angle=angle),
          legend.position = c(0.025, 1), legend.justification = c('left', 'top'),
          legend.key.size = unit(0.5, 'cm')
    )
  split_p = all_data %>% filter(dataset == 'gnomAD') %>%
    ggplot + aes(x=pop, y=count, fill=pop) + 
    geom_bar(stat='identity') + theme_classic() + 
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_y_continuous(breaks=c(seq(0, 130000, 10000), 141456)) + ylab('') + xlab('') + 
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_text(margin = ggplot2::margin(t = -10), size = 15, angle=angle),
          legend.position = c(0.025, 1), legend.justification = c('left', 'top'),
          legend.key.size = unit(0.5, 'cm')
    )
  total_eas = subpop_data %>% filter(pop == 'eas') %$% sum(count)
  total_eur = subpop_data %>% filter(pop == 'eur') %$% sum(count)
  eur_p = subpop_data %>% filter(pop == 'eur') %>%
    ggplot + aes(x=pop, y=count, fill=subpop) + 
    geom_bar(stat='identity') + theme_classic() + 
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_y_continuous(breaks=c(seq(0, 70000, 10000), total_eur)) + ylab('') + xlab('') + 
    scale_x_discrete(labels=c("eur" = "gnomAD EUR")) +
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_text(margin = ggplot2::margin(t = -10), size = 15, angle=angle),
          # legend.position = c(0.025, 1), legend.justification = c('left', 'top'),
          # legend.key.size = unit(0.5, 'cm')
    )
  eas_p = subpop_data %>% filter(pop == 'eas') %>%
    ggplot + aes(x=pop, y=count, fill=subpop) + 
    geom_bar(stat='identity') + theme_classic() + 
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_y_continuous(breaks=c(seq(0, 8000, 2000), total_eas)) + ylab('') + xlab('') + 
    scale_x_discrete(labels=c("eas" = "gnomAD EAS")) +
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_text(margin = ggplot2::margin(t = -5), size = 15, angle=angle),
          # legend.position = c(0.025, 1), legend.justification = c('left', 'top'),
          # legend.key.size = unit(0.5, 'cm')
    )
  p = full_p + (eas_p / eur_p + plot_layout(heights = c(1, 3))) + plot_layout(widths = c(4, 1))
  # split_p + (eas_p / eur_p + plot_layout(heights = c(1, 3))) + plot_layout(widths = c(4, 1))
  if (return) {
    return(p)
  } else {
    print(p)
  }
  if (save_plot) dev.off()
}
plot_panel_size(save_plot=T)

subset_data = bind_rows(
  data.frame(dataset='genomes', subset="Non-neuro", afr=1694, amr=277, asj=123, eas=780, eur=582+6813, oth=367),
  data.frame(dataset='genomes', subset="Non-cancer", eas=780, eur=7718+1738, asj=145, afr=4359, amr=424, oth=513+31),
  data.frame(dataset='genomes', subset="Controls", afr=1287, amr=123, asj=19, eas=458, eur=581+2762, oth=212),
  data.frame(dataset='genomes', subset="Non-TOPMed", afr=4278, amr=405, asj=69, eas=761, eur=1738+5547, oth=506),
  data.frame(dataset='exomes', subset="Non-neuro", afr=8109, amr=15262, asj=3106, eas=6708, eur=8367+44779, oth=2433, sas=15304),
  data.frame(dataset='exomes', subset="Non-cancer", afr=7451, amr=17130, asj=4786, eas=8846, eur=10816+51377, oth=2810, sas=15263),
  data.frame(dataset='exomes', subset="Controls", afr=3582, amr=8556, asj=1160, eas=4523, eur=6697+21384, oth=957, sas=7845),
  data.frame(dataset='exomes', subset="Non-TOPMed", afr=6013, amr=17229, asj=4999, eas=9195, eur=10823+55840, oth=3032, sas=15308)
) %>% 
  gather(pop, count, -dataset, -subset) %>% filter(!is.na(count)) %>%
  mutate(pop = fct_relevel(pop, 'oth', after=Inf),
         subset = fct_relevel(subset, 'Non-TOPMed', 'Non-cancer', 'Non-neuro', 'Controls'))

dataset_colors = c('exomes' = '#4682B4', 'genomes' = '#73AB3D')

png('gnomad_barplot_subsets.png',height=1500,width=2650,res=300)
subset_data %>%
  ggplot + aes(x=subset, y=count, fill=dataset, color=dataset) +
  geom_bar(stat='identity') + theme_classic() +
  scale_fill_manual(name='', values=dataset_colors, labels=c('exomes', 'genomes')) +
  scale_color_manual(name='', values=dataset_colors, labels=c('exomes', 'genomes')) +
  scale_y_continuous(breaks=c(seq(0, 140000, 20000))) + ylab('') + xlab('') +
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        axis.text.y = element_text(size = 30, angle=0),
        axis.text.x = element_text(margin = ggplot2::margin(t = -10), size = 20, angle=-30),
        legend.justification = c('right', 'top'), legend.position = c(1, 1),
        legend.key.size = unit(0.8, 'cm'),
        legend.text = element_text(size=20),
        plot.margin = margin(r = 20, unit = "pt")
  ) + coord_flip()
dev.off()


plot_panel_size_angle = function(save_plot=F, return=F, angle=0) {
  all_data %>%
    ggplot + aes(x=dataset, y=count, fill=pop) + 
    geom_bar(stat='identity') + theme_classic() + 
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_y_continuous(breaks=seq(0, 140000, 10000)) + ylab('') + xlab('') + 
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_text(margin = ggplot2::margin(t = -10), angle = angle),
          legend.position = c(0.025, 1), legend.justification = c('left', 'top'),
          legend.key.size = unit(0.5, 'cm')
    )
  if (return) {
    return(p)
  } else {
    print(p)
  }
}
plot_panel_size_anne = function() {
  for (i in seq(0, 360, 10)) {
    print(plot_panel_size_angle(angle=i, return=T))
  }
}
saveGIF(plot_panel_size_anne(), movie.name = 'rotating_sample_sizes.gif', interval=0.1)

function() {
  all_data %>%
    ggplot + aes(x=dataset, y=count, fill=pop, frame=as.numeric(dataset), cumulative=T) + 
    geom_bar(stat='identity') + theme_classic() + 
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_y_continuous(breaks=seq(0, 140000, 10000)) + ylab('') + xlab('') + 
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_text(margin = ggplot2::margin(t = -10)),
          legend.position = c(0.025, 1), legend.justification = c('left', 'top'),
          legend.key.size = unit(0.5, 'cm')
    ) -> p
  animation::ani.options(ani.width = 1800, ani.height = 1200, interval = 0.5, loop=T, ani.res = 300)#, ani.dev = 'pdf', ani.type = 'pdf')
  gganimate(p, interval=0.5, title_frame = F, 'bars.gif')
  
}