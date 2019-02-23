source('constants.R')

model_colors = c(
  RF="#1f77b4",
  VQSR="#ff7f0e"
)
model_names = c(
  RF="Random Forests",
  VQSR="Variant Quality Score Recalibration"
)

n_trios = c(
  exomes=4568,
  genomes=212
)

# True = snv, False = indel
cutoffs = c(
  'TRUE'=90,
  'FALSE'=82
)

plot_titles = c(
  'exomes True' = "Exomes SNVs",
  'exomes False' = "Exomes Indels",
  'genomes True' = "Genomes SNVs",
  'genomes False' = "Genomes Indels"
)

labels = c('a','b','c','d','e','f','g','h','i','j','k','l','m')

# Load data
autosomes = read_tsv(get_or_download_file('autosomes.tsv.gz', subfolder = 'variant_qc/'))
chr20 = read_tsv(get_or_download_file('chr20.tsv.gz', subfolder = 'variant_qc/'))
concordance = read_tsv(get_or_download_file('concordance.tsv.gz', subfolder = 'variant_qc/'))

# General format
format_supp_plot = function(p, title=NA){
  p = p +
    scale_color_manual(name='', values=model_colors, labels=model_names) +
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) + 
    ggtitle("") +
    theme(legend.position = "None")
  return(p)
}

# Concordance plots format
format_concordance = function(p, same_scale){
  p =  p +
    scale_x_continuous(labels=percent_format(accuracy = 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if(same_scale){
    p = p +
      scale_y_continuous(labels=percent_format(accuracy = 1),  limits = c(0.8, 1.0))
  } else {
    p = p +
      scale_y_continuous(labels=percent_format(accuracy = .1))
  }
  return(p)
}

# Get a row
get_plot_row = function(plots, xlabel, ylabel, labels){
  plot = ggarrange(plotlist=plots, nrow=1, ncol=4,  labels=labels)
  return( plot %>%  
            annotate_figure(
              left = text_grob(ylabel, rot = 90),
              bottom = text_grob(xlabel)
            )
  )
}

get_header_row = function(nrows_y_label){
  row = ggarrange(
    text_grob("Exomes SNVs", face='bold', just='center'),
    text_grob("Exomes Indels", face='bold', just='center'),
    text_grob("Genomes SNVs", face='bold', just='center'),
    text_grob("Genomes Indels", face='bold', just='center'),
    nrow = 1,
    ncol = 4
  )
  
  if(nrows_y_label>0){
    row = annotate_figure(row, left=strrep("\n",nrows_y_label-1))
  }
  return(row)
}

get_legend_row = function(){
  p = ggplot(data.frame(x=c(0,1), y=c(0,1), model=c('RF','VQSR'))
             ,aes(x,y,col=model)) +  
    geom_point() +
    scale_color_manual(name='', values=model_colors, labels=model_names) +
    theme(legend.position="bottom")
  plot = ggarrange(get_legend(p))
  plot = ggarrange(plotlist = list(ggplot() + geom_blank(), get_legend(p),  ggplot() + geom_blank()), ncol=3)
  
  return(plot %>% annotate_figure(
    left = text_grob("\n", rot = 90)
  )
  )
}

## PR Plot
pr_plot = function(same_scale = F){
  truth_samples = c('NA12878' = 'NA12878', 'syndip' = 'Synthetic diploid')
  rows = list(get_header_row(2))
  for(s in names(truth_samples)){
    row=list()
    for(x in c('exomes', 'genomes')){
      for(y in c(TRUE, FALSE)){
        plot_data = concordance %>%
          filter(data_type == x & 
                   snv == y & 
                   rank_name == "truth_sample_rank" & 
                   truth_sample == s)
        p =  ggplot(plot_data,  aes(x=recall, y=precision, col=model)) +
          geom_point() + 
          geom_point(data=plot_data %>% filter(model=='RF'), aes(x=recall, y=precision, col=model)) +
          geom_vline(xintercept = plot_data %>% filter(model == 'RF' & bin==cutoffs[as.character(y)]) %$% recall, 
                     linetype='dashed') +
          geom_hline(yintercept = plot_data %>% filter(model == 'RF' & bin==cutoffs[as.character(y)]) %$% precision, 
                     linetype='dashed')
        p = format_supp_plot(p)
        row = c(row,list(format_concordance(p, same_scale)))
      }
    }
    rows = c(rows,
             list(get_plot_row(row,
                               paste(truth_samples[s],'recall'), 
                               paste(truth_samples[s],'\nprecision',sep=''),
                               labels=labels[(4*(length(rows)-1)+1):(4*(length(rows)-1)+4)]
             )))
  }
  rows = c(rows, list(get_legend_row()))
  return(ggarrange(plotlist=rows, nrow=length(rows), ncol = 1, heights = c(0.1, 1, 1, 0.2)))
}

# Rare variants metrics plots
rare_variants_metrics_plot = function(){
  
  #dnms
  dnms = autosomes %>%
    filter(rank_id == 'rank') %>%
    group_by(data_type, model, snv, bin) %>%
    summarise(n_de_novo=sum(n_de_novo)) %>%
    group_by(data_type, model, snv) %>%
    arrange(bin) %>%
    mutate(cum_dnm = cumsum(n_de_novo)/n_trios[data_type])
  
  dnm_plots=list()
  for(x in c('exomes', 'genomes')){
    for(y in c(TRUE, FALSE)){
      plot_data = dnms %>%
        filter(data_type == x & 
                 snv == y)
      p = ggplot(plot_data, aes(bin, cum_dnm, col=model))  + 
        geom_point() + 
        geom_point(data=plot_data %>% filter(model=='RF'), aes(bin, cum_dnm, col=model)) +
        geom_vline(xintercept=cutoffs[as.character(y)], linetype='dashed')
      dnm_plots = c(dnm_plots,
                    list(format_supp_plot(p))
      )
    }
  }
  dnm_row = get_plot_row(dnm_plots, 
                         'Model percentile', 
                         bquote(atop(~ italic('De novo') ~ 'calls per child', '(cumulative)')), 
                         labels = labels[1:4]
  )
  
  #trans singletons
  trans_singletons = chr20 %>%
    filter(rank_id == 'rank') %>%
    group_by(data_type, model, snv, bin) %>%
    summarise(n_trans_singletons=sum(n_trans_singletons)) %>%
    group_by(data_type, model, snv) %>%
    arrange(bin) %>%
    mutate(cum_trans_singletons = cumsum(n_trans_singletons))
  
  trans_singletons_plots=list()
  for(x in c('exomes', 'genomes')){
    for(y in c(TRUE, FALSE)){
      plot_data = trans_singletons %>% filter(data_type == x & snv == y)
      p = ggplot(plot_data, aes(bin, cum_trans_singletons, col=model))  + 
        geom_point() + 
        geom_point(data=plot_data %>% filter(model=='RF'), aes(bin, cum_trans_singletons, col=model)) +
        geom_vline(xintercept=cutoffs[as.character(y)], linetype='dashed')
      trans_singletons_plots=  c(trans_singletons_plots,
                                 list(format_supp_plot(p))
      )
    }
  }
  trans_singletons_row = get_plot_row(trans_singletons_plots, 
                                      'Model percentile', 
                                      'Transmitted singletons\n(chrom 20, cumulative)', 
                                      labels = labels[5:8]
  )
  
  #trans singletons
  validated_dnm = autosomes %>%
    filter(rank_id == 'rank') %>%
    group_by(data_type, model, snv, bin) %>%
    summarise(n_validated_de_novos=sum(n_validated_de_novos)) %>%
    group_by(data_type, model, snv) %>%
    arrange(bin) %>%
    mutate(cum_validated_de_novos = cumsum(n_validated_de_novos))
  
  validated_dnm_plots=list()
  for(x in c('exomes')){
    for(y in c(TRUE, FALSE)){
      plot_data = validated_dnm %>% filter(data_type == x & snv == y)
      p = ggplot(plot_data, aes(bin, cum_validated_de_novos, col=model))  + 
        geom_point() + 
        geom_point(data=plot_data %>% filter(model=='RF'), aes(bin, cum_validated_de_novos, col=model)) +
        geom_vline(xintercept=cutoffs[as.character(y)], linetype='dashed')
      validated_dnm_plots=  c(validated_dnm_plots,
                              list(format_supp_plot(p))
      )
    }
  }
  validated_dnm_row = get_plot_row(validated_dnm_plots, 
                                   'Model percentile', 
                                   bquote(atop('Validated' ~ italic('de novo') ~ 'mutations','(cumulative)')), 
                                   labels = labels[9:10]
  )
  plot = ggarrange(
    get_header_row(2),
    dnm_row, 
    trans_singletons_row, 
    validated_dnm_row,
    get_legend_row(),
    nrow = 5, 
    ncol = 1,
    heights = c(0.1, 1, 1, 1, 0.2))
  
  return(plot)
}

efigure2 = function() {
  e2 = pr_plot()
  
  pdf("extended_data_figure2.pdf", width = 8, height = 5)
  print(e2)
  dev.off()
  
  png("extended_data_figure2.png", width = 8, height = 5, units = 'in', res=300)
  print(e2)
  dev.off()
  
}

efigure3 = function() {
  e3 = rare_variants_metrics_plot()
  
  pdf("extended_data_figure3.pdf", width = 8, height = 7)
  print(e3)
  dev.off()
  
  png("extended_data_figure3.png", width = 8, height = 7, units = 'in', res=300)
  print(e3)
  dev.off()
  
}
