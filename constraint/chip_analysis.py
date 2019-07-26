#!/usr/bin/env python

# the purpose of this script is to test for clonal hematopoiesis by checking for differences in the age and allele balance distibutions for LoF and synonymous variants

from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
from scipy import stats
import numpy as np
import argparse


logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("chip")
logger.setLevel(logging.INFO)


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


def test_ks(a1,a2):
    """
    Runs the two sample Kolmogorov-Smirnov test on the two supplied arrays, requires ndarrays with at least 2 distinct values
    :param list a1: array 1
    :param list a2: array 2
    :return: [KS statisic, KS p-value]
    :rtype: list
    """
    if isinstance(a1,np.ndarray) and isinstance(a2,np.ndarray):
        if len(set(a1)) > 1 and len(set(a2)) > 1:  # require at least 2 distinct values
            ks_result = stats.ks_2samp(a1,a2)
            return pd.Series((ks_result[0], ks_result[1]))
        else:
            return pd.Series((np.NaN, np.NaN))
    else:
        return pd.Series((np.NaN, np.NaN))
    
    
def test_mood(a1,a2):
    """
    Runs the two sample Kolmogorov-Smirnov test on the two supplied arrays, requires ndarrays with at least 2 distinct values
    :param list a1: array 1
    :param list a2: array 2
    :return: [KS statisic, KS p-value]
    :rtype: list
    """
    if isinstance(a1,np.ndarray) and isinstance(a2,np.ndarray):
        if len(set(a1)) > 1 and len(set(a2)) > 1:  # require at least 2 distinct values
            stat,p,grand_median,contigency_table = stats.median_test(a1,a2)
            return p
        else:
            return pd.Series((np.NaN, np.NaN))
    else:
        return pd.Series((np.NaN, np.NaN))


def main(args):
    
    # upload gnomad sites ht for exomes 
    logger.info("Reading in hail table...")
    ht = hl.read_table(public_exomes_ht_path())

    # explode hail table
    ht = ht.explode(ht.vep.transcript_consequences)

    # filter to pass-only variants
    ht = ht.filter(ht.filters.length() == 0, keep=True)
    
    # filter to canonical-only variants
    ht = ht.filter(ht.vep.transcript_consequences.canonical == 1, keep=True)
    
    
    # pull out one example of bin edges and take the floor
    logger.info("Obtaining bin floors...")
    bin_edges_ab = ht.aggregate(hl.agg.take(ht.ab_hist_alt.bin_edges, 1))
    bin_floors_ab = bin_edges_ab[0][:-1]

    bin_edges_age = ht.aggregate(hl.agg.take(ht.age_hist_het[0].bin_edges, 1))
    bin_floors_age = bin_edges_age[0]  # don't use -1 because want to keep upper bound as floor of n_larger
    bin_floors_age.insert(0,25)  # add bin floor of 25 for <30 age bin found in n_smaller
    
    
    # aggregate histograms and group by transcript ID and variant severity
    logger.info("Aggregating histograms...")
    ht = ht.annotate(chrom = ht.locus.contig)
    
    # set the variant severity based on consequence and lof annotations (only set to lof if high confidence)
    ht = ht.annotate(variant_severity = (hl.case()
                                         .when(hl.is_defined(ht.vep.transcript_consequences.lof) & (ht.vep.transcript_consequences.lof == "HC" ), "lof")
                                         .when(ht.vep.transcript_consequences.consequence_terms.contains("synonymous_variant"), "synonymous")
                                         .when(ht.vep.transcript_consequences.consequence_terms.contains("missense_variant"), "missense")
                                         .default("none")))

    # sum het and hom age histograms per variant
    ht = ht.annotate(age_middle = ht.age_hist_het[0].bin_freq + ht.age_hist_hom[0].bin_freq,
                           age_n_smaller = ht.age_hist_het[0].n_smaller + ht.age_hist_hom[0].n_smaller,
                           age_n_larger = ht.age_hist_het[0].n_larger + ht.age_hist_hom[0].n_larger)

    # add n_smaller and n_larger to ends of age histogram
    ht = ht.annotate(age_hist = hl.array([ht.age_n_smaller]).extend(ht.age_middle).append(ht.age_n_larger))
                      
    # aggregate histograms and group by transcript ID and variant severity
    table_result = (ht.group_by(ht.chrom, ht.vep.transcript_consequences.transcript_id, ht.vep.transcript_consequences.gene_symbol, ht.vep.transcript_consequences.gene_id, ht.variant_severity).
                    aggregate(summed_ab_hists = hl.agg.array_sum(ht.ab_hist_alt.bin_freq),
                              summed_age_hists = hl.agg.array_sum(ht.age_hist)))

    # compare to a median of the summed histograms as an alternative test (take median of summed histograms for synonymous variants in transcripts)
    median_ab_hist = table_result.aggregate(hl.agg.filter(table_result.variant_severity == "synonymous", hl.agg.array_agg(lambda x: hl.median(hl.agg.collect(x)), table_result.summed_ab_hists)))
    median_age_hist = table_result.aggregate(hl.agg.filter(table_result.variant_severity == "synonymous", hl.agg.array_agg(lambda x: hl.median(hl.agg.collect(x)), table_result.summed_age_hists)))
    
    
    # convert to pandas dataframe
    logger.info("Converting to pandas dataframe...")
    table_result = table_result.key_by()
    tr_ht = table_result.select('chrom', 'transcript_id', 'variant_severity', 'gene_symbol', 'gene_id', 'summed_ab_hists', 'summed_age_hists')
    df = tr_ht.to_pandas()

    # test lof vs synonymous distributions for ab and age using the median histogram data
    expanded_hist_ab_median = expand_hist(median_ab_hist,bin_floors_ab)
    expanded_hist_age_median = expand_hist(median_age_hist,bin_floors_age)

    subsets_to_combine = []
    for chromosome in set(df['chrom']):
        logger.info("Processing chromosome {chromosome}...".format(**locals()))
        logger.info("Expanding histograms...")
        subset_df = df[df['chrom'] == chromosome]
        subset_df['expanded_hists_ab'] = (subset_df['summed_ab_hists']).apply(expand_hist, args=[bin_floors_ab])
        subset_df['expanded_hists_age'] = (subset_df['summed_age_hists']).apply(expand_hist, args=[bin_floors_age])
        
        
        logger.info("Merging gene names and pivoting table...")
        # preserve gene symbols/ids matched to transcript id
        transcript_to_gene = subset_df[['transcript_id', 'gene_symbol', 'gene_id']].drop_duplicates()

        # spread pandas dataframe
        subset_df = subset_df.pivot(index='transcript_id', columns='variant_severity')[['expanded_hists_ab', 'expanded_hists_age', 'summed_ab_hists', 'summed_age_hists']]
        subset_df['transcript_id'] = subset_df.index

        # add back in gene symbol
        subset_df = pd.merge(subset_df, transcript_to_gene[['transcript_id','gene_symbol', 'gene_id']], on='transcript_id', how='left')
        
        # rename columns
        subset_df.columns = map(rename_cols, subset_df.columns)
    
    
        logger.info("Running ks test...")
        # test lof vs synonymous distributions for ab and age
        subset_df[['ks_statistic_ab','ks_p_value_ab']] = subset_df[['expanded_hists_ab_lof','expanded_hists_ab_synonymous']].apply(lambda x: test_ks(*x), axis=1)
        subset_df[['ks_statistic_age','ks_p_value_age']] = subset_df[['expanded_hists_age_lof','expanded_hists_age_synonymous']].apply(lambda x: test_ks(*x), axis=1)
        subset_df[['ks_statistic_ab_median','ks_p_value_ab_median']] = (subset_df['expanded_hists_ab_lof']).apply(test_ks, args=[expanded_hist_ab_median])
        subset_df[['ks_statistic_age_median','ks_p_value_age_median']] = (subset_df['expanded_hists_age_lof']).apply(test_ks, args=[expanded_hist_age_median])

        # sort by p-value for ks test on ab
        subset_df = subset_df.sort_values(by=['ks_p_value_ab'])
    
    
        logger.info("Running Mood's median test...")
        # test lof vs synonymous distributions for ab and age
        subset_df['mood_p_ab'] = subset_df[['expanded_hists_ab_lof','expanded_hists_ab_synonymous']].apply(lambda x: test_mood(*x), axis=1)
        subset_df['mood_p_age'] = subset_df[['expanded_hists_age_lof','expanded_hists_age_synonymous']].apply(lambda x: test_mood(*x), axis=1)

        
        logger.info("Calculating length of arrays...")
        # get lengths of arrays that were used for ks test
        subset_df['len_expanded_hists_ab_lof'] = subset_df['expanded_hists_ab_lof'].apply(lambda x: len(x) if isinstance(x,np.ndarray) else np.NaN)
        subset_df['len_expanded_hists_ab_synonymous'] = subset_df['expanded_hists_ab_synonymous'].apply(lambda x: len(x) if isinstance(x,np.ndarray) else np.NaN)
        subset_df['len_expanded_hists_age_lof'] = subset_df['expanded_hists_age_lof'].apply(lambda x: len(x) if isinstance(x,np.ndarray) else np.NaN)
        subset_df['len_expanded_hists_age_synonymous'] = subset_df['expanded_hists_age_synonymous'].apply(lambda x:  len(x) if isinstance(x,np.ndarray) else np.NaN)

    
        logger.info("Calculating medians...")
        subset_df['med_ab_lof'] = subset_df['expanded_hists_ab_lof'].apply(lambda x: np.median(x) if isinstance(x,np.ndarray) else np.NaN)
        subset_df['med_ab_synonymous'] = subset_df['expanded_hists_ab_synonymous'].apply(lambda x: np.median(x) if isinstance(x,np.ndarray) else np.NaN)
        subset_df['med_age_lof'] = subset_df['expanded_hists_age_lof'].apply(lambda x: np.median(x) if isinstance(x,np.ndarray) else np.NaN)
        subset_df['med_age_synonymous'] = subset_df['expanded_hists_age_synonymous'].apply(lambda x: np.median(x) if isinstance(x,np.ndarray) else np.NaN)

        # retain only certain columns
        df_out = subset_df[['transcript_id',
                 'gene_symbol',
                 'gene_id',
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
                 'mood_p_ab',
                 'mood_p_age',
                 'med_ab_lof',
                 'med_ab_synonymous',
                 'med_age_lof',
                 'med_age_synonymous']]

        subsets_to_combine.append(df_out)
    
    # concatenate results from each chromosome
    df = pd.concat(subsets_to_combine)
    
    df['len_expanded_hist_ab_median'] = len(expanded_hist_ab_median)
    df['len_expanded_hist_age_median'] = len(expanded_hist_age_median)  

    # write out results file
    with hl.hadoop_open(args.output, 'w') as out:
        df.to_csv(out, sep='\t')


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='The purpose of this script is to test for clonal hematopoiesis by checking for differences in the age and allele balance distibutions for LoF and synonymous variants')
    parser.add_argument('-o', '--output', required=True, help='path/name for output file')
    
    args = parser.parse_args()
    main(args)

    

