#!/usr/bin/env python

# the purpose of this script is to test for clonal hematopoiesis by checking for differences in the age and allele balance distibutions for LoF and synonymous variants

import hail as hl
import hail.expr.aggregators as agg
from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
import pandas as pd
import numpy as np
from scipy import stats


# upload gnomad sites vcf for exomes 
ht = hl.read_table(public_exomes_ht_path())

# subset to chromosome 20 to create test dataset
ht20 = ht.filter(ht.locus.contig == "20")
ht20e = ht20.explode(ht20.vep.transcript_consequences)

# filter to pass-only variants
ht20e = ht20e.filter(ht20e.filters.length() == 0, keep=True)

# pull out one example of bin edges and take the floor
bin_edges_ab = ht20e.aggregate(hl.agg.take(ht20e.ab_hist_alt.bin_edges, 1))
bin_floors_ab = bin_edges_ab[0][:-1]

bin_edges_age = ht20e.aggregate(hl.agg.take(ht20e.age_hist_het[0].bin_edges, 1))
bin_floors_age = bin_edges_age[0] # don't use -1 because want to keep upper bound as floor of n_larger
bin_floors_age.insert(0,25) # add bin floor of 25 for <30 age bin found in n_smaller


# aggregate histograms and group by transcript ID and variant severity

# set the variant severity based on consequence and lof annotations (only set to lof if high confidence)
ht20e = ht20e.annotate(variant_severity = hl.cond(ht20e.vep.transcript_consequences.consequence_terms.contains("synonymous_variant"),
                                                  "synonymous",
                                                  hl.cond(ht20e.vep.transcript_consequences.consequence_terms.contains("missense_variant"),
                                                          "missense",
                                                          "none")))
               
ht20e = ht20e.annotate(variant_severity = hl.cond(hl.is_missing(ht20e.vep.transcript_consequences.lof),
                                                  ht20e.variant_severity,
                                                  hl.cond(ht20e.vep.transcript_consequences.lof == "HC",
                                                          "lof",
                                                          ht20e.variant_severity)))
# sum het and hom age histograms per variant
ht20e = ht20e.annotate(age_middle = np.add(ht20e.age_hist_het[0].bin_freq, ht20e.age_hist_hom[0].bin_freq),
                       age_n_smaller = hl.sum([ht20e.age_hist_het[0].n_smaller, ht20e.age_hist_hom[0].n_smaller]),
                       age_n_larger = hl.sum([ht20e.age_hist_het[0].n_larger, ht20e.age_hist_hom[0].n_larger]))

# add n_smaller and n_larger to ends of age histogram
ht20e = ht20e.annotate(age_hist = hl.array([ht20e.age_n_smaller]).extend(ht20e.age_middle).append(ht20e.age_n_larger))
                                         
# aggregate histograms and group by transcript ID and variant severity
table_result = (ht20e.group_by(ht20e.vep.transcript_consequences.transcript_id, ht20e.vep.transcript_consequences.gene_symbol, ht20e.variant_severity).
                aggregate(summed_ab_hists = hl.agg.array_agg(lambda x: hl.agg.sum(x), ht20e.ab_hist_alt.bin_freq),
                          summed_age_hists = hl.agg.array_agg(lambda x: hl.agg.sum(x), ht20e.age_hist)))

# compare to a median of the summed histograms as an alternative (take median of summed histograms for synonymous variants in transcripts)
median_ab_hist = table_result.aggregate(hl.agg.filter(table_result.variant_severity == "synonymous", hl.agg.array_agg(lambda x: hl.median(hl.agg.collect(x)), table_result.summed_ab_hists)))
median_age_hist = table_result.aggregate(hl.agg.filter(table_result.variant_severity == "synonymous", hl.agg.array_agg(lambda x: hl.median(hl.agg.collect(x)), table_result.summed_age_hists)))


# convert to pandas dataframe
table_result = table_result.key_by()
tr_ht = table_result.select('transcript_id','variant_severity','gene_symbol','summed_ab_hists','summed_age_hists')
df = tr_ht.to_pandas()


# expand histograms and assign to bin floor values
def expand_hist(repeat_numbers,bin_numbers):
    """
    Generates an expanded histogram with the values of bin_numbers repeated x times, where x is specified in repeat_numbers
    :param list repeat_numbers: list of number of times to repeat the respective elements in bin_numbers
    :param list bin_numbers: list of values to repeat
    :return: list of the expanded histogram
    :rtype: list
    """
    summed_hist = np.repeat(bin_numbers,repeat_numbers)
    return summed_hist

df['expanded_hists_ab'] = (df['summed_ab_hists']).apply(expand_hist, args=[bin_floors_ab])
df['expanded_hists_age'] = (df['summed_age_hists']).apply(expand_hist, args=[bin_floors_age])

# preserve gene symbols matched to transcript id
transcript_to_gene = df[['transcript_id','gene_symbol']].drop_duplicates()

# spread pandas dataframe
df = df.pivot(index='transcript_id', columns='variant_severity')[['expanded_hists_ab', 'expanded_hists_age', 'summed_ab_hists', 'summed_age_hists']]
df['transcript_id'] = df.index

# add back in gene symbol
df = pd.merge(df, transcript_to_gene[['transcript_id','gene_symbol']], on='transcript_id', how='left')

# rename tuple columns
def rename_cols(column):
    """
    Reformats tuple column names to str format
    :param col: column name
    :return: column name in str format
    :rtype: str
    """
    if isinstance(column, tuple):
        column = '_'.join(str(x) for x in column)
    return column
df.columns = map(rename_cols, df.columns)


# run Kolmogorov-Smirnov test
def test_ks(a1,a2):
    """
    Runs the two sample Kolmogorov-Smirnov test on the two supplied arrays, requires ndarrays with a length of at least 1
    :param list a1: array 1
    :param list a2: array 2
    :return: [KS statisic, KS p-value]
    :rtype: list
    """
    if isinstance(a1,np.ndarray) and isinstance(a2,np.ndarray):
        if len(a1) > 0 and len(a2) > 0: # should set a length threshold here
            ks_result = stats.ks_2samp(a1,a2)
            return pd.Series((ks_result[0], ks_result[1]))
        else:
            return pd.Series((np.NaN, np.NaN))
    else:
        return pd.Series((np.NaN, np.NaN))

# test lof vs synonymous distributions for ab and age
df[['ks_statistic_ab','ks_p_value_ab']] = df[['expanded_hists_ab_lof','expanded_hists_ab_synonymous']].apply(lambda x: test_ks(*x), axis=1)
df[['ks_statistic_age','ks_p_value_age']] = df[['expanded_hists_age_lof','expanded_hists_age_synonymous']].apply(lambda x: test_ks(*x), axis=1)

# test lof vs synonymous distributions for ab and age using the median histogram data
expanded_hist_ab_median = expand_hist(median_ab_hist,bin_floors_ab)
expanded_hist_age_median = expand_hist(median_age_hist,bin_floors_age)

df[['ks_statistic_ab_median','ks_p_value_ab_median']] = (df['expanded_hists_ab_lof']).apply(test_ks, args=[expanded_hist_ab_median])
df[['ks_statistic_age_median','ks_p_value_age_median']] = (df['expanded_hists_age_lof']).apply(test_ks, args=[expanded_hist_age_median])

# sort by p-value for ks test on ab
df = df.sort_values(by=['ks_p_value_ab'])

# get lengths of arrays that were used for ks test
df['len_expanded_hists_ab_lof'] = df['expanded_hists_ab_lof'].apply(lambda x: len(x) if isinstance(x,np.ndarray) else np.NaN)
df['len_expanded_hists_ab_synonymous'] = df['expanded_hists_ab_synonymous'].apply(lambda x: len(x) if isinstance(x,np.ndarray) else np.NaN)
df['len_expanded_hists_age_lof'] = df['expanded_hists_age_lof'].apply(lambda x: len(x) if isinstance(x,np.ndarray) else np.NaN)
df['len_expanded_hists_age_synonymous'] = df['expanded_hists_age_synonymous'].apply(lambda x:  len(x) if isinstance(x,np.ndarray) else np.NaN)
df['len_expanded_hist_ab_median'] = len(expanded_hist_ab_median)
df['len_expanded_hist_age_median'] = len(expanded_hist_age_median)


# retain only certain columns
df_out = df[['transcript_id',
             'gene_symbol',
             'summed_ab_hists_lof',
             'summed_age_hists_lof',
             'summed_ab_hists_synonymous',
             'summed_age_hists_synonymous',
             'ks_statistic_ab',
             'ks_p_value_ab',
             'ks_statistic_age',
             'ks_p_value_age',
             'ks_statistic_ab_median',
             'ks_p_value_ab_median',
             'ks_statistic_age_median',
             'ks_p_value_age_median',
             'len_expanded_hists_ab_lof',
             'len_expanded_hists_ab_synonymous',
             'len_expanded_hists_age_lof',
             'len_expanded_hists_age_synonymous',
             'len_expanded_hist_ab_median',
             'len_expanded_hist_age_median']]



