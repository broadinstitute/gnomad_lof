gene_data = load_constraint_data()
transcript_data = load_constraint_data(level = 'transcript')
tx_summary = load_tx_summary()
tx_disease_gene = load_gene_tissue_mapping()

produce_stats <- function(net, measures = NA, quickrun = T, runsilent = F){
  #Stats for the network, can supply list of igraph or tidygraph functions to run
  if(all(is.na(measures))){
    measures <- c('centrality_degree()',
                  'centrality_betweenness()',
                  'centrality_closeness()',
                  'centrality_eigen()',
                  'centrality_hub()')
    if(!quickrun){
      #These take some time - so only run if specifically asked for
      measures <- c(measures,
                    'centrality_katz()',
                    'centrality_pagerank()',
                    'centrality_random_walk()',
                    'centrality_power()',
                    'centrality_alpha()')
    }
  }
  name <- gsub("\\(.*\\)",'', measures)
  
  for(i in 1:length(measures)){
    if(!runsilent){
      cat('Working on ', name[i] ,'...\n')
    }
    tryCatch(
      net %<>% activate(nodes) %>%
        mutate(!!name[i] := !!parse_expr(measures[[i]])),
      error=function(e) {
        message(e, '\n')
        return(paste0(measures[i], ' failed'))
      })
  }
  return(net)
}

load_degree_data = function() {
  tmp_dir = "./tmpStringDB"
  suppressWarnings(dir.create(tmp_dir))
  string_db <- STRINGdb$new(version="10",
                            species=9606,
                            score_threshold=0,
                            input_directory=tmp_dir)
  
  # Wrangling -----------------------------------------------------------------
  gene_data_map = gene_data %>% 
    as.data.frame %>% #function doesn't work with tibbles
    string_db$map('transcript', removeUnmappedRows = F) %>%
    as_tibble()
  
  ints <- string_db$get_interactions(gene_data_map$STRING_id) %>% 
    select(from, to, combined_score) %>% 
    as_tibble()
  
  # Subset to high confidence
  ints.hc <- ints %>%
    filter(combined_score >= 700)
  # Init net
  net <- as_tbl_graph(ints.hc)
  
  #Stats
  net %<>% produce_stats() 
  # Multiple components, reduce to largest and rerun stats. 
  net <- to_components(net)[[1]]
  net %<>% produce_stats() %>%
    mutate(weighted_degree = centrality_degree() / local_ave_degree(),
           local_transitivity = local_transitivity())
  
  nds <- net %N>% as_tibble 
  nds %>% left_join(gene_data_map, by=c('name' = 'STRING_id')) %>%
    select(name, gene, transcript, everything()) %>%
    return
}
degree_data = load_degree_data()

mean_degree = function(save_plot=F) {
  #' ## Results
  #' 
  #' Check degree distribution
  
  cent.dat <- degree_data %>%
    filter(exp_lof > 10) %>%
    group_by(oe_lof_upper_bin) %>%
    summarise(Mean = mean(centrality_degree, na.rm = T),
              SD = sd(centrality_degree, na.rm = T),
              Median = median(centrality_degree, na.rm = T),
              medMAD = mad(centrality_degree),
              Mean.Betweenness = mean(log1p(centrality_betweenness), na.rm = T),
              Median.Betweenness = median(log1p(centrality_betweenness),
                                          na.rm = T),
              medMAD.bet = mad(log1p(centrality_betweenness)),
              medlocalTran = median(local_transitivity, na.rm = T),
              meanlocalTran = mean(local_transitivity, na.rm = T),
              total= n()) 
  
  degree_data %>%
    filter(exp_lof > 10) %>%
    lm(centrality_degree ~ oe_lof_upper + cds_length, .) %>%
    summary %>%
    print
  
  p = cent.dat %>%
    ggplot + aes(x = oe_lof_upper_bin, y = Mean, 
                 ymin = Mean - 1.96 * SD / sqrt(total), 
                 ymax = Mean + 1.96 * SD / sqrt(total)) +
    # geom_point(color = color_lof) +
    geom_pointrange(color = color_lof) + 
    theme_classic() + oe_x_axis + ylab('Mean number of\nprotein-protein interactions')
  
  if (save_plot) {
    pdf('4a_networks.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

constrained_v_unconstrained_transcript = function(save_plot=F, oe_threshold = NULL) {
  if (is.null(oe_threshold)) {
    oe_threshold = transcript_data %>% filter(oe_lof_upper_bin == 0) %$% max(oe_lof_upper)
  }
  
  # Summarize genes in terms of how many constrained and unconstrained transcripts they have 
  constraint_summary = transcript_data %>% 
    group_by(gene_id) %>%
    summarize(n_constrained_transcript = sum(oe_lof_upper < oe_threshold),
              n_unconstrained_transcript = sum(oe_lof_upper >= oe_threshold), 
              n_total_transcript = n()) %>% ungroup
  
  # Filter to genes where there is one constrained and one uncontrained transcripts 
  genes_with_constrained_unconstrained_transcripts = constraint_summary %>%
    filter(n_constrained_transcript >= 1 & n_unconstrained_transcript >= 1)
  print(paste('Found', nrow(genes_with_constrained_unconstrained_transcripts),
              'with at least one constrained and unconstrained at LOEUF =',
              oe_threshold))  # 2182 genes
  
  transcript_data %>%
    filter(gene_id %in% genes_with_constrained_unconstrained_transcripts$gene_id) %>%
    count(gene_id, constrained = oe_lof_upper < oe_threshold) %>%
    group_by(constrained) %>% 
    summarize(mean_num_transcripts = mean(n))
  # For genes with a constrained and unconstrained transcript, there is a mean of ~3.2 and 3.5, respectively
  
  transcript_data %>%
    filter(gene_id %in% genes_with_constrained_unconstrained_transcripts$gene_id) %>%
    left_join(tx_summary) %>%
    group_by(gene_id, transcript, cds_length,
             constrained = ifelse(oe_lof_upper < oe_threshold, 'constrained', 'unconstrained')) %>%
    summarize(mean_expression=mean(mean_expression, na.rm=T)) %>%
    lm(mean_expression ~ as.numeric(constrained == 'constrained') + cds_length, .) %>%
    summary %>%
    print
  
  transcript_data %>%
    filter(gene_id %in% genes_with_constrained_unconstrained_transcripts$gene_id) %>%
    left_join(tx_summary) %>%
    group_by(gene_id,
             constrained = ifelse(oe_lof_upper < oe_threshold, 'constrained', 'unconstrained')) %>%
    summarize(mean_cds_length = mean(cds_length, na.rm=T), mean_expression=mean(mean_expression, na.rm=T)) %>%
    # group_by(constrained) %>% summarize(mean(mean_expression))
    group_by(gene_id) %>% 
    summarize(delta = mean(mean_expression*(constrained == 'constrained')) - mean(mean_expression*(constrained == 'unconstrained'))) %>%
    summarize(mean_delta = mean(delta))
    # ggplot + aes(x = delta) + geom_density()
    
  
  # Histogram of difference
  for_paired_plot = transcript_data %>%
    filter(gene_id %in% genes_with_constrained_unconstrained_transcripts$gene_id) %>%
    left_join(tx_summary) %>%
    group_by(gene_id, 
             constrained = ifelse(oe_lof_upper < oe_threshold, 'constrained', 'unconstrained')) %>%
    # sample_n(1) %>%
    summarize(mean_expression=mean(mean_expression, na.rm=T)) %>%
    mutate(mean_expression = ifelse(mean_expression > 10, 10, mean_expression)) %>% 
    ungroup %>%
    spread('constrained', 'mean_expression')
  
  # Difference histogram
  p = for_paired_plot %>%
    filter(constrained + unconstrained > 0.5) %>%  # either is expressed at at least 1
    ggplot + aes(x = constrained/(constrained + unconstrained)) + 
    # ggplot + aes(x = (constrained - unconstrained)/(constrained + unconstrained)) + 
    # xlab("Difference in expression of constrained and uncontrained transcript / Sum of expression of the two transcripts") +
    geom_histogram(color = "steelblue", fill = "steelblue", alpha = 0.5) + 
    theme_classic() + ylab('Number of genes') + scale_x_continuous(labels=percent) +
    xlab("Percent of expression\nfrom constrained transcript")
  
  if (save_plot) {
    pdf('4_constrained_v_unconstrained_transcript.pdf', height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

plot_expression_metrics_by_constraint = function(save_plot=F, only_canonical=T, metric='n_tissues_expressed',
                                                 canonical_all_comparison=F, x_axis='oe_lof_upper_bin') {
  all_tx = transcript_data %>%
    left_join(tx_summary) %>%
    rename(transcript_oe_lof_upper_bin = oe_lof_upper_bin) %>%
    left_join(gene_data %>% select(gene_id, oe_lof_upper_bin)) %>%
    mutate_(metric = metric,
            x_axis = x_axis)
  
  if (metric == 'mean_expression') {
    all_tx %<>% mutate(metric=ifelse(metric > 5, 5, metric))
  }
  if (x_axis == 'transcript_oe_lof_upper_bin') {
    xlabel = paste(oe_x_label, '(transcript)')
  } else {
    xlabel = oe_x_label
  }
  
  if (!canonical_all_comparison) {
    if (metric == 'n_tissues_expressed') {
      ylabel = paste0('Number of tissues where\n',
                      if_else(only_canonical, 'canonical ', ''), 'transcript is expressed')
    } else {
      ylabel = paste('Mean', if_else(only_canonical, 'canonical ', ''), 'transcript expression')
    }
    
    lm(n_tissues_expressed ~ oe_lof_upper + cds_length, all_tx) %>% summary
    
    plot_data = all_tx %>%
      filter(!only_canonical | canonical)
    
    lm(n_tissues_expressed ~ oe_lof_upper + cds_length, plot_data) %>% summary
    
    cor.test(plot_data$metric, plot_data$oe_lof_upper)
    cor.test(plot_data$metric, plot_data$oe_lof_upper, method='spearman')
    all_tx %>%
      lm(metric ~ oe_lof_upper, .) %>% 
      summary
  
    collapsed_tx = plot_data %>%
      group_by(x_axis) %>%
      summarize(mean_metric = mean(metric, na.rm=T),
                sd_metric = sd(metric, na.rm=T),
                n = n(),
                sem_metric = 1.96 * sd_metric / sqrt(n))
  
    p = ggplot() +
      geom_violin(aes(x = x_axis, y = metric, group = x_axis),
                  data = plot_data,
                  fill = color_lof, alpha = 0.2, color = F) +
      theme_classic() + oe_x_axis + xlab(xlabel) +
      geom_pointrange(aes(x = x_axis, y = mean_metric,
                          ymin = mean_metric - sem_metric,
                          ymax = mean_metric + sem_metric),
                      data = collapsed_tx, color = color_lof) +
      ylab(ylabel)
  } else {
    if (metric == 'n_tissues_expressed') {
      ylabel = 'Number of tissues where\ntranscript is expressed'
    } else {
      ylabel = 'Mean transcript expression'
    }

    canonical_plot_data = all_tx %>%
      filter(canonical) %>% 
      group_by(x_axis) %>% 
      summarize(mean_metric = mean(metric, na.rm=T),
                sd_metric = sd(metric, na.rm=T),
                n = n(),
                sem_metric = 1.96 * sd_metric / sqrt(n),
                transcripts='canonical')
    all_plot_data = all_tx %>%
      group_by(x_axis) %>% 
      summarize(mean_metric = mean(metric, na.rm=T),
                sd_metric = sd(metric, na.rm=T),
                n = n(),
                sem_metric = 1.96 * sd_metric / sqrt(n),
                transcripts='all')
    
    full_plot_data = canonical_plot_data %>%
      bind_rows(all_plot_data)
    
    tx_colors = c('all' = 'dodgerblue1',
                  'canonical' = 'navy')
    p = full_plot_data %>% 
      ggplot() + theme_classic() +
      aes(x = x_axis, y = mean_metric, 
          ymin = mean_metric - sem_metric, 
          ymax = mean_metric + sem_metric,
          color = transcripts) +
      geom_pointrange() +
      oe_x_axis + xlab(xlabel) +
      scale_color_manual(values=tx_colors, guide=F) +
      ylab(ylabel) + 
      annotate('text', x = 9.5, y = max(full_plot_data$mean_metric), hjust = 1, label='Canonical transcripts', 
                color=tx_colors[['canonical']]) +
      annotate('text', x = 9.5, y = max(full_plot_data$mean_metric) * 0.92, hjust = 1, label='All transcripts', 
                color=tx_colors[['all']])
  }
  
  if (save_plot) {
    pdf(paste0('4b_', metric, '.pdf'), height=3, width=4)
    print(p)
    dev.off()
  }
  return(p)
}

figure4 = function() {
  p4a = mean_degree()
  p4b = plot_expression_metrics_by_constraint()
  p4c = constrained_v_unconstrained_transcript()
  pdf('figure4.pdf', height=7, width=3.5)
  print(ggarrange(p4a, p4b, p4c, ncol = 1, nrow = 3, labels = 'auto', align = 'v', vjust = 1))
  dev.off()
  png('figure4.png', height=7., width=3.5, units = 'in', res=300)
  print(ggarrange(p4a, p4b, p4c, ncol = 1, nrow = 3, labels = 'auto', align = 'v', vjust = 1))
  dev.off()
}
