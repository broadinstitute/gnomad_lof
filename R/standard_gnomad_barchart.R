source('constants.R')

# Make stacked barplot for exomes vs genomes

input_datasets = bind_rows(
  data.frame(dataset='1000\nGenomes', genomes=T, eas=523, sas=494, eur=514, afr=691, amr=355),
  data.frame(dataset='ESP', genomes=F, eur=4298, afr=2217),
  data.frame(dataset='ExAC', genomes=F, eas=4327, sas=8256, eur=33370+3307, afr=5203, amr=5789, unk=454),
  data.frame(dataset='gnomAD\nv2', genomes=F, eas=9197, sas=15308, eur=56885+10824, asj=5040, afr=8128, amr=17296, unk=3070),
  data.frame(dataset='gnomAD\nv2', genomes=T, eas=780, eur=7718+1738, asj=145, afr=4359, amr=424, unk=513+31),
#  data.frame(dataset='discovEHR', genomes=F, unk=50000),
  data.frame(dataset='GME', genomes=F, mde=2497),
  data.frame(dataset='BRAVO', genomes=T, unk=62784),
  data.frame(dataset='gnomAD\nv3', genomes=T, afr=21042, amr=6835, asj=1662, eas=1567, eur=32299+5244+450, unk=1077, sas=1526), # eur = nfe + fin + ami
  data.frame(dataset='dummy', genomes=T, unk=0)
)

all_data = input_datasets %>%
  gather(pop, count, -dataset, -genomes) %>% filter(!is.na(count)) %>%
  group_by(dataset, genomes, pop) %>% summarize(count=sum(count)) %>% ungroup %>%
  mutate(dataset = fct_reorder(dataset, count, .fun=sum),
         pop = fct_relevel(pop, 'oth', 'unk', after=Inf))

bars = input_datasets %>%
  gather(pop, count, -dataset, -genomes) %>% filter(!is.na(count)) %>%
  group_by(dataset, genomes) %>% summarize(count=sum(count)) %>% ungroup %>% group_by(genomes) %>%
  arrange(count, .by_group=T)

# Compute breaks for 'Exomes' and 'Genomes' labels
labels = gsub("dummy", "", bars$dataset)
dummy_idx = which(bars$dataset == "dummy")
nbars = dim(bars)[1]
label_breaks = c(dummy_idx/2, (nbars+dummy_idx+1)/2)

pop_names['onf'] = 'Other European'
pop_names['unk'] = 'Unknown/Not Shared'
pop_colors_border = pop_colors
pop_colors_border['oth'] = 'gray'
pop_colors['oth'] = 'white'


pdf('gnomad_barplot_compare.pdf', width=14)
	all_data %>%
    mutate(dataset=fct_reorder(interaction(dataset, genomes), count, sum)) %>%
    ggplot + aes(x=interaction(dataset, genomes), y=count, fill=pop) +
    geom_bar(stat='identity', size=1.1) + theme_classic(base_size=24) +
    scale_fill_manual(name='', values=pop_colors, labels=pop_names) +
    scale_y_continuous(breaks=c(seq(0, 130000, 10000), 141456)) + ylab('') + xlab('') +
    annotate("text", label_breaks, y = -25000, label = c("Exomes", "Genomes"), size=7, fontface = "bold") +
    annotate("text", seq_len(length(labels)), y=-10000, label=labels, size=6) +
    theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_blank(), axis.ticks.x = element_blank(),
          axis.ticks.length = unit(0.3, "cm"),
          axis.text.x = element_blank(),
          legend.position = "bottom",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-25,0,0,0)
    )
dev.off()