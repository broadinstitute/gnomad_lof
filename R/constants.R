# install.packages(c('tidyverse', 'Hmisc', 'devtools', 'broom', 'plotly'))
# install.packages(c('slackr', 'magrittr', 'gapminder', 'readr', 'purrr'))
# install.packages(c('skimr', 'gganimate', 'gghighlight', 'plotROC', 'naniar'))
# install.packages(c('BiocManager', 'cowplot', 'corrplot', 'corrr', 'ggridges'))
# install.packages(c('ggpubr', 'meta', 'tidygraph', 'pbapply', 'RMySQL'))
# BiocManager::install('STRINGdb', version = '3.8')
# devtools::install_github('VPetukhov/ggrastr')
# devtools::install_github('hafen/trelliscopejs')
# devtools::install_github('thomasp85/patchwork')

options(stringsAsFactors = F)
library(Hmisc)
library(RCurl)
library(plyr)
library(tidyverse)
library(broom)
library(plotly)
library(slackr)
library(magrittr)
library(scales)
library(trelliscopejs)
library(gapminder)
library(skimr)
library(gghighlight)
library(plotROC)
library(naniar)
library(patchwork) # devtools::install_github("thomasp85/patchwork")
library(corrplot)
library(corrr)
library(readxl)
library(ggridges)
library(grid)
library(ggpubr)
library(ggrastr)
library(grDevices)
library(meta)
library(STRINGdb)
library(tidygraph)
library(rlang)
library(pbapply)
library(ggrepel)
library(RMySQL)
library(cowplot)
library(ggwordcloud)

data_dir = './data/'
suppressWarnings(dir.create(data_dir))
data_versions = c('v1.1', 'v1.0')
get_data_url = function(version = 'v1.0') {
  return(paste0('https://storage.googleapis.com/gnomad-public/papers/2019-flagship-lof/', version, '/'))
}

slackr_setup(channel='#constraint')
post_slack = FALSE

ggplot2 = function(...) {
  return(ggplot(...) + theme_bw() + theme(panel.border=element_blank(), legend.key = element_blank()))
}
color_syn = '#AAAAAA'
color_mis = '#FF6103'
color_os = '#74099E'
color_lof = '#9D1309'
color_lc_lof = '#EE799F'
var_type_aliases = c('syn' = 'Synonymous', 'mis' = 'Missense', 'lof' = 'pLoF', 'os' = 'Other splice')
colors = c('synonymous_variant' = color_syn,
           'missense_variant' = color_mis,
           'stop_gained' = color_lof,
           'Synonymous' = color_syn,
           'Missense' = color_mis,
           'synonymous' = color_syn,
           'missense' = color_mis,
           'nonsense' = color_lof,
           'Other splice' = color_os,
           'LoF' = color_lof,
           'pLoF' = color_lof)
linetypes = c('Expected' = 'dashed',
              'Observed' = 'solid')
mut_cols = c("transversion" = "#EA4444", "non-CpG transition" = "#458B00", "CpG transition" = '#2E9FFE',
             'ACG' = '#422efe', 'CCG' = '#2e73fe', 'GCG' = '#2ec3fe', 'TCG' = '#2efef7')
dataset_colors = c('ExAC' = '#4682B4', 'gnomAD' = '#73AB3D', 
                   'exomes' = '#4682B4', 'genomes' = '#73AB3D')

color_amr = k_amr = '#ED1E24'
color_eur = k_eur = '#6AA5CD'
color_afr = k_afr = '#941494'
color_sas = k_sas = '#FF9912'
color_eas = k_eas = '#108C44'
color_oth = k_oth = '#ABB9B9'
color_mde = k_mde = '#33CC33'
color_asj = k_asj = 'coral'

color_nfe = k_nfe = color_eur
color_fin = k_fin = '#002F6C'

color_syn = k_syn = '#AAAAAA'
color_mis = k_mis = '#FF6103'
color_lof = k_lof = '#9D1309'

color_cpg = '#2E9FFE'
color_ti = '#458B00'
color_tv = '#EA4444'

lof_like = c('frameshift_variant','essential_splice','stop_gained','splice_donor_variant','splice_acceptor_variant')
mis_like = c('missense_variant','inframe_indel','stop_lost',
             'mature_miRNA_variant','start_lost')
syn_like = c('synonymous_variant','3_prime_UTR_variant','5_prime_UTR_variant', 'splice_region_variant',
             'extended_splice','stop_retained_variant','non_coding_transcript_exon_variant',
             'upstream_gene_variant', 'downstream_gene_variant',
             'intron_variant','intergenic_variant','regulatory_region_variant')

format_vep_category = function(category_list) {
  return(category_list %>%
           gsub("_"," ", .) %>%
           gsub('stop gained', 'nonsense', .) %>%
           gsub("inframe indel", "in-frame indel", .) %>%
           gsub("initiator codon", "start lost", .) %>%
           gsub(" variant", "", .) %>%
           gsub("transcript exon", "transcript", .) %>%
           gsub(" prime ","'", .) %>%
           gsub("probably damaging", "prob damaging", .) %>%
           gsub("possibly damaging", "poss damaging", .))
}

variant_category_colors = c(rep(color_lof, length(lof_like)),
                            rep(color_mis, length(mis_like)),
                            rep(color_syn, length(syn_like)))
names(variant_category_colors) = c(format_vep_category(lof_like),
                                   format_vep_category(mis_like),
                                   format_vep_category(syn_like))

variant_type_colors = c(color_tv, color_ti, color_cpg, color_cpg)
names(variant_type_colors) = c('transversion', 'non-CpG transition', 'CpG transition', 'CpG')

pops <- c('afr', 'amr', 'eas', 'fin', 'nfe', 'sas')
pop_colors = c('afr' = color_afr,
               'amr' = color_amr,
               'eas' = color_eas,
               'fin' = color_fin,
               'eur' = color_nfe,
               'nfe' = color_nfe,
               'oth' = color_oth,
               'sas' = color_sas,
               'mde' = color_mde,
               'asj' = color_asj,
               'uniform' = 'pink',
               'consanguineous' = 'pink',
               'sas_non_consang' = 'orange',
               'exac' = 'gray',
               'est' = 'black',
               'bgr' = '#66C2A5',
               'est' = '#4891D9',
               'nwe' = '#C60C30',
               'seu' = '#3ca021',
               'swe' = 'purple',
               'onf' = color_nfe,
               'kor' = '#4891D9',
               'jpn' = '#BC002D',
               'oea' = color_eas,
               'unk' = color_oth)
pop_names = c('oth' = 'Other',
              'afr' = 'African/African-American',
              'amr' = 'Latino',
              'eas' = 'East Asian',
              'fin' = 'Finnish',
              'eur' = 'European',
              'nfe' = 'European',
              'sas' = 'South Asian',
              'mde' = 'Middle Eastern',
              'asj' = 'Ashkenazi Jewish',
              'uniform' = 'Uniform',
              'sas_non_consang' = 'South Asian (F < 0.05)',
              'consanguineous' = 'South Asian (F > 0.05)',
              'exac' = 'ExAC',
              'bgr' = 'Bulgarian',
              'est' = 'Estonian',
              'nwe' = 'Northern European',
              'seu' = 'Southern European',
              'swe' = 'Swedish',
              'onf' = 'Other Non-Finnish European',
              'kor' = 'Korean',
              'jpn' = 'Japanese',
              'oea' = 'Other East Asian',
              'unk' = 'Unknown')

agilent = c("Exome Express", "Standard Exome Sequencing", "Standard Exome Sequencing v2",
            "WHOLE EXOME HYB SEL & SEQ", "HYB SEL & SEQ")
ice = c("Standard Exome Sequencing v3", "Standard Exome Sequencing v4", 
        "Nextera Exome", "G4L WES + Array v1", "Express Human WES (Standard Coverage) v1",
        "Exome Express v2", "Exome Express v3",
        "Express Human WES (Deep Coverage) v1")
ice150 = c("G4L WES + Array v2", "Standard Exome Sequencing v5") # ICE with 150bp reads
wgs = c('PCR Free', 'PCR Plus')

platform_names = c('gnomAD' = 'gnomAD',
                   'agilent' = 'Agilent',
                   'ice' = 'Illumina',
                   'ice150' = 'Illumina (150 bp reads)',
                   'nimblegen' = 'Nimblegen',
                   'multiple' = 'Multiple',
                   'unknown' = 'Unknown')

platform_alphas = c('agilent' = 0.8, 'ice' = 0.8, 'ice150' = 0.8, 'nimblegen' = 0.8, 'gnomAD' = 0.8, 'multiple' = 0.5, 'unknown' = 0.5)
platform_colors = c('agilent' = '#F8766D', 'ice' = '#00BA38', 'ice150' = '#619CFF', 'nimblegen' = '#B79F00', 'gnomAD' = '#F564E3', 'multiple' = '#F564E3', 'unknown' = 'gray')

# gradient_colors = c('#FF2600', '#FF9300', '#FFC000', '#E1C977', '#C2D1ED', '#A8E577')
# gradient_values = c(0, 0.4, 0.6, 0.8, 1, 1.6)

color_benign = '#87A4DC'
color_dominant = '#F0810F'
color_recessive = '#FFBB00'
gradient_colors = c('#FF2600', '#FF9300', '#FFC000', '#E1C977', '#C2D1ED')
gradient_values = c(0, 0.4, 0.6, 0.8, 1)
gradient_value_deciles = c(0, 2, 5, 7, 9)

gradient_colors2 = c('#FF2600', '#FF9300', '#FFC000', '#E1C977', '#C2D1ED', '#A8E577')
gradient_values2 = c(0, 0.4, 0.6, 0.8, 1, 1.2)

gene_list_colors = c(
  'Haploinsufficient' = color_lof,
  'Essential Genes' = '#CB0000',
  'Autosomal Dominant' = color_dominant,
  'Autosomal Recessive' = color_recessive,
  'Olfactory Genes' = color_benign,  # '#598234',
  'Background' = 'lightgray',
  'Universe' = 'lightgray'
)

annotate_variant_types = function(data) {
  return(data %>%
           mutate(transition = (ref == 'A' & alt == 'G') | (ref == 'G' & alt == 'A') | (ref == 'C' & alt == 'T') | (ref == 'T' & alt == 'C'),
                  cpg = (ref == 'C' & alt == 'T' & substr(context, 2, 3) == 'CG') | (ref == 'G' & alt == 'A' & substr(context, 1, 2) == 'CG'),
                  variant_type = case_when(
                    !transition ~ 'transversion',
                    cpg ~ 'CpG transition',
                    !cpg & transition ~ 'non-CpG transition',
                    TRUE ~ 'wut'
                  )))
}

generate_ci = function(obs, exp, alpha=0.05) {
  distr = data.frame(x = seq(0, 2, by=0.001)) %>%
    mutate(y = dpois(obs, x*exp),
           y_cumulative = cumsum(y)) %>%
    mutate(y_cumulative_norm = y_cumulative/max(y_cumulative))
  
  low =  distr %>% filter(y_cumulative_norm < alpha) %$% max(x)
  high = distr %>% filter(y_cumulative_norm > 1 - alpha) %$% min(x)
  return(data.frame(lower=low, upper=high))
}


constraint_metric_name = 'LOEUF'
oe_x_label = paste(constraint_metric_name, 'decile')
label_function = function(x) {
  y = 10 * x + 5
  ifelse(y %% 20 == 0, paste0(y, '%'), "")
}
oe_x_axis = list(xlab(oe_x_label),
                 scale_x_continuous(labels=label_function, breaks=seq(-0.5, 9.5, 1), limits=c(-0.5, 9.5)))

get_metric_common_name = function(x) {
  return(case_when(x == 'oe_lof' ~ 'Observed/Expected',
                   x == 'oe_lof_upper' ~ constraint_metric_name,
                   TRUE ~ x))
}

label_function2 = function(x) {
  paste0(x*10, '-', x*10 + 10, '%')
}

ds_breaks = c(0, 4e4, 8e4, 125748)
ds_breaks_log = c(1e2, 1e3, 1e4, 125748)
downsampling_x_axis = function(log=T) {
  if (log) {
    list(xlab('Sample size'), scale_x_log10(labels=comma, breaks=ds_breaks_log))
  } else {
    list(xlab('Sample size'), scale_x_continuous(labels=comma, breaks=pds_breaks))
  }
}

get_or_download_file = function(base_fname, subfolder='', local_name='', version='') {
  fname = paste0(data_dir, ifelse(local_name != '', local_name, base_fname))
  if (!file.exists(fname)) {
    if (version == '') {
      for (version in data_versions) {
        url = paste0(get_data_url(version), subfolder, base_fname)
        if (url.exists(url)) {
          break
        }
      }
    }
    url = paste0(get_data_url(version), subfolder, base_fname)
    download.file(url, fname)
  }
  return(fname)
}

get_de_novo_data = function() {
  base_fname = 'asd_ddid_de_novos.tx_annotated.021819.tsv.bgz'
  fname = paste0(data_dir, base_fname)  
  if (!file.exists(fname)) {
    url = paste0('https://storage.googleapis.com/gnomad-public/papers/2019-tx-annotation/results/de_novo_variants/', base_fname)
    download.file(url, fname)
  }
  new_de_novo_data = read_delim(gzfile(fname), delim = '\t', 
                                col_types = list(CHROM=col_character(), locus=col_character()))
  
  fname = get_or_download_file('ASC_v15_callset1.5_de_novo_calls_with_annotations.raw.2017-12-13.txt',
                               subfolder = 'misc_files/')
  quality_data = read_delim(fname, delim = '\t')
  
  new_de_novo_data %>%
    left_join(quality_data %>% transmute(CHROM, POSITION, REF, ALT, filters = Filters)) %>%
    filter(is.na(filters) | filters == 'PASS') %>%
    return
}

load_downsampled_gene_data = function() {
  fname = get_or_download_file('gnomad.v2.1.1.lof_metrics.downsamplings.txt.bgz')
  return(read_delim(gzfile(fname), delim = '\t', col_types = list(canonical=col_logical())))
}

load_downsampled_data = function() {
  ds_data = load_downsampled_gene_data()
  
  collapsed_ds = ds_data %>% filter(canonical) %>%
    group_by(downsampling, pop) %>%
    summarize_at(c('exp_syn', 'obs_syn', 'exp_mis', 'obs_mis', 'exp_lof', 'obs_lof', 'n_sites', 'caf'),
                 .funs = list(sum = function(x) sum(x, na.rm=T),
                              over0 = function(x) sum(x >= 0, na.rm=T),
                              total = function(x) length(x),
                              over5raw = function(x) sum(x >= 5, na.rm=T)/length(x),
                              over10raw = function(x) sum(x >= 10, na.rm=T)/length(x),
                              over20raw = function(x) sum(x >= 20, na.rm=T)/length(x),
                              over5 =  function(x) sum(x >= 5, na.rm=T)/sum(x > 0, na.rm=T),
                              over10 = function(x) sum(x >= 10, na.rm=T)/sum(x > 0, na.rm=T),
                              over20 = function(x) sum(x >= 20, na.rm=T)/sum(x > 0, na.rm=T))
    ) %>% ungroup
  
  ds_melt = collapsed_ds %>% filter(pop == 'global') %>% select(-pop) %>%
    gather(key='metric', value='count', -downsampling) %>%
    separate(metric, into = c('obs_exp', 'variant_type', 'func')) %>%
    mutate(variant_type = if_else(variant_type == 'sites', NA_character_, variant_type)) %>%
    mutate(variant_type = fct_recode(variant_type, 'Synonymous' = 'syn', 'Missense' = 'mis', 'LoF' = 'lof'),
           obs_exp = fct_recode(obs_exp, 'Expected' = 'exp', 'Observed' = 'obs'))
  return(ds_melt)
}

load_constraint_data = function(level='gene', loftee=T) {
  fname = paste0('gnomad.v2.1.1.lof_metrics.by_', level, '.txt.bgz')
  subfolder = ''
  local_name = ''
  if (!loftee) {
    subfolder = 'other_cuts/no_loftee/'
    local_name = paste0('gnomad.v2.1.1.lof_metrics.no_loftee.by_', level, '.txt.bgz')
  }
  fname = get_or_download_file(fname, subfolder, local_name)
  if (level == 'transcript') {
    col_list = list(canonical=col_logical())
  } else {
    col_list = list()
  }
  read_delim(gzfile(fname), delim = '\t',
             col_types = do.call(cols, col_list)) %>%
    mutate(s = mu_lof / p) %>%
    return
}

get_ko_gene_lists = function(list_dir = 'ko_gene_lists/') {
  gene_lists = map_df(list.files(list_dir, 'list_.+\\.tsv'), 
                      function(x) read_tsv(paste(list_dir, x, sep = '/'), 
                                           col_names = F, col_types = cols()
                                           ) %>% 
                        transmute(gene = toupper(X1), gene_list=str_sub(x, 6, -5)))
  return(gene_lists)
}

load_gene_tissue_mapping = function() {
  transcript_data = load_constraint_data(level = 'transcript')
  transcript_data = transcript_data %>%
    filter(!is.na(oe_lof_upper)) %>%
    group_by(gene) %>%
    mutate(most_constrained_transcript = min(oe_lof_upper) == oe_lof_upper) %>% ungroup
  
  gene_tissue_fname = get_or_download_file('NIHMS950518-supplement-table_s2.xlsx', subfolder = 'misc_files/')
  gene_tissue = read_excel(gene_tissue_fname, skip=4) %>%
    filter(!is.na(gene))
  
  fname = get_or_download_file('GTEx.v7.median_expression_per_tx_per_tissue.021018.tsv.bgz', subfolder = 'misc_files/')
  tx_expression = read_delim(gzfile(fname), delim = '\t') %>%
    select(-v, -agg_expression) %>%
    rename(transcript = transcript_id)
  
  tissue_mapping = data.frame(gtex_tissue = sort(names(tx_expression)[3:55]),
                              tissue = c(rep("Adipose_Tissue", 2),
                                         "Adrenal_Gland", rep("Blood_Vessel",3),"Bladder",
                                         rep("Brain",13), "Breast", rep("Cells",2),
                                         rep("Cervix",2), rep("Colon",2),
                                         rep("Esophagus",3), "Fallopian",
                                         rep("Heart",2), "Kidney", "Liver",
                                         "Lung", "MinorSalivaryGland", "Muscle",
                                         "Nerve","Ovary", "Pancreas","Pituitary","Prostate",
                                         rep("Skin", 2), "Small_Intestine",
                                         "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", 'Blood'))
  
  most_expressed_transcript_by_gene_tissue = tx_expression %>%
    gather('gtex_tissue', 'expression', -transcript, -gene_id) %>%
    left_join(tissue_mapping) %>%
    select(-gtex_tissue) %>%
    group_by(transcript, gene_id, tissue) %>% summarize(expression = max(expression)) %>%  # comment out to remove collapsing by tissue
    group_by(tissue, gene_id) %>%
    mutate(most_expressed = expression == max(expression)) %>% 
    ungroup %>%
    left_join(transcript_data %>% select(transcript, gene_id, gene))
  
  gene_tissue_mapping = gene_tissue %>% select(-mim.ids, -mim.titles, -num.affected.tissues) %>%
    gather('tissue', 'num_diseases', -gene) %>%
    mutate(tissue = gsub(" ", "_", tissue)) %>%
    filter(num_diseases > 0)
  
  most_expressed_transcript_by_gene_tissue_disease_only = most_expressed_transcript_by_gene_tissue %>% 
    inner_join(gene_tissue_mapping)
  
  tx_disease_gene = transcript_data %>%
    left_join(most_expressed_transcript_by_gene_tissue_disease_only)
  
  return(tx_disease_gene)
}

load_maps_data = function(cut = 'plain', data_type = 'exomes') {
  col_list = list()
  if (cut == 'plain') {
    col_list[['protein_coding']] = col_logical()
  }
  fname = get_or_download_file(paste0('maps_', cut, '_', data_type, '.txt.bgz'), subfolder = 'summary_results/')
  read_delim(gzfile(fname), delim = '\t',
             col_types = do.call(cols, col_list)) %>%
    mutate(csq=format_vep_category(worst_csq),
           maps_upper=maps + 1.96 * maps_sem,
           maps_lower=maps - 1.96 * maps_sem)
}

regroup_maps = function(data, maps_grouping) {
  maps = data %>%
    group_by_at(vars(maps_grouping)) %>% 
    summarize(singleton_count=sum(singleton_count),
              expected_singletons=sum(expected_singletons),
              variant_count=sum(variant_count),
              ps=singleton_count / variant_count,
              maps=(singleton_count - expected_singletons)/variant_count,
              maps_sem=sqrt(ps * (1 - ps) / variant_count),
              maps_upper=maps + 1.96 * maps_sem,
              maps_lower=maps - 1.96 * maps_sem) %>% ungroup
  return(maps)
}

cumulative_maps = function(data, prefix='maps') {
  maps = data %>%
    mutate(temptemp_singleton_count=cumsum(singleton_count),
           temptemp_expected_singletons=cumsum(expected_singletons),
           temptemp_variant_count=cumsum(variant_count),
           temptemp_ps=temptemp_singleton_count / temptemp_variant_count,
           temptemp=(temptemp_singleton_count - temptemp_expected_singletons)/temptemp_variant_count,
           temptemp_sem=sqrt(temptemp_ps * (1 - temptemp_ps) / temptemp_variant_count),
           temptemp_upper=temptemp + 1.96 * temptemp_sem,
           temptemp_lower=temptemp - 1.96 * temptemp_sem) %>% 
    rename_at(vars(temptemp_singleton_count:temptemp_lower), 
              function(x) paste0(prefix, str_sub(as.character(x), 9)))
  return(maps)
}

load_observed_possible_data = function(data_type='exomes') {
  fname = paste0('observed_possible_expanded_', data_type, '.txt.bgz')
  fname = get_or_download_file(fname, subfolder = 'summary_results/')
  read_delim(gzfile(fname), 
             delim = '\t', col_types=list(possible=readr::col_number(),
                                          observed=readr::col_number(),
                                          singletons=readr::col_number()
             )) %>%
    mutate(csq=format_vep_category(worst_csq)) %>%
    return
}

load_observed_possible_sites_data = function(data_type='exomes') {
  fname = paste0('observed_possible_sites_', data_type, '.txt')
  fname = get_or_download_file(fname, subfolder = 'summary_results/', version='v1.1')
  read_delim(fname, delim = '\t', col_types=list(coverage=readr::col_number())) %>%
    return
}

load_indel_data = function(data_type='exomes') {
  fname = paste0('indels_summary_', data_type, '.txt.bgz')
  fname = get_or_download_file(fname, subfolder = 'summary_results/')
  read_delim(gzfile(fname), delim = '\t') %>%
    mutate(csq=format_vep_category(worst_csq)) %>%
    return
}

load_freq_summary_data = function() {
  fname = get_or_download_file('freq_loftee_exomes.txt.bgz', subfolder = 'summary_results/')
  read_delim(gzfile(fname), delim = '\t',
             col_types = list(fail_loftee = col_logical(),
                              loftee_os = col_logical(),
                              pass_loftee = col_logical(),
                              pass_loftee_with_flags = col_logical())) %>%
    return
}

load_sv_data = function() {
  fname = get_or_download_file('gnomAD-SV_v2_rev1_releasable_forFlagship.rare_biallelic_LoF_deletions_per_gene.obs_exp.txt', 
                               subfolder = 'misc_files/', version = 'v1.1')
  read_delim(fname, delim = '\t') %>%
    rename(gene = '#gene') %>%
    return
}

load_omim_by_year_data = function() {
  omim_data = read_delim('data/forKonrad_cleaned_gene_discovery_years_2018-10-09.tsv', delim = '\t') %>%
    filter(yearDiscovered > 1990)
  
  # omim_data %>%
  #   count(yearDiscovered) %>%
  #   ggplot + aes(x = yearDiscovered, y = n) + geom_bar(stat='identity')
  
  omim_data %>%
    group_by(Ensembl.Gene.ID) %>%
    summarize(year = min(yearDiscovered), 
              discoverybyNGS = any(discoverybyNGS)
    ) %>%
    return
}

load_omim_data = function() {
  omim_data = read_delim('../../gene_lists/other_data/omim.use.tsv', delim = '\t',
                         col_types = list(chromosome = col_character())) %>%
    unnest(gene = strsplit(genes, "|", fixed=T)) %>%
    return
}

load_all_gene_list_data = function(list_dir = 'gene_lists/lists/') {
  all_files = list.files(list_dir, '.+tsv')
  if (length(all_files) == 0) {
    if (length(all_files) == 0) system('git clone https://github.com/macarthur-lab/gene_lists.git')
    all_files = list.files(list_dir, '.+tsv')
  }
  gene_lists = map_df(all_files, 
                      function(x) read_tsv(paste(list_dir, x, sep = '/'), 
                                           col_names = F, col_types = cols()
                      ) %>% 
                        transmute(gene = toupper(X1), gene_list=str_sub(x, end=-5)))
  return(gene_lists)
}

load_sum_stats = function() {
  fname = get_or_download_file('reshA3.tsv.gz', subfolder = 'misc_files/')
  read_tsv(fname, col_types = list(var_id=col_character(), prevalence=col_number())) %>%
    mutate(name = as.factor(name)) %>%
    mutate(description = gsub('Schizoprenia', 'Schizophrenia', gsub('\\*', '', description)))
}

load_tx_summary = function(expression_cutoff=0.3) {
  fname = get_or_download_file('GTEx.v7.median_expression_per_tx_per_tissue.021018.tsv.bgz', subfolder = 'misc_files/')
  tx_expression = read_delim(gzfile(fname), delim = '\t') %>%
    select(-c(v, agg_expression, Bladder, Brain_Spinalcord_cervicalc_1_, Brain_Substantianigra,
              Cervix_Ectocervix,Cervix_Endocervix, FallopianTube, Kidney_Cortex,
              MinorSalivaryGland, Uterus, Ovary,Testis, Vagina,
              Cells_EBV_transformedlymphocytes, Cells_Transformedfibroblasts, Prostate)) %>%
    rename(transcript = transcript_id)
  
  tx_melt = tx_expression %>%
    gather('tissue', 'expression', -transcript, -gene_id)
  
  tx_melt %>%
    group_by(transcript, gene_id) %>%
    summarize(mean_expression=mean(expression, na.rm=T),
              brain_mean_expression=mean(expression * ifelse(grepl('Brain', tissue), 1, NA), na.rm=T),
              n_tissues_expressed=sum(expression > expression_cutoff, na.rm=T),
              max_expression=max(expression, na.rm=T),
              max_expressed_tissue=tissue[which.max(expression)]
    ) %>% ungroup %>% 
    return
}

load_go_data = function() {
  read_delim('../misc_files/goa_human.gaf.gz', delim = '\t', comment = '!', 
             col_names = c("DB", "DB_Object_ID", "DB_Object_Symbol", "DB_Object_Name",
                           "DB_Object_Synonym", "DB_Object_Type", "Taxon",
                           "Parent_Object_ID", "DB_Xref", "Properties"),
             col_types=list(X16 = col_character())) %>%
    return
}
