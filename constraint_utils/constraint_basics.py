__author__ = 'konradk'

from .generic import *
import statsmodels.formula.api as smf
import pickle

root = 'gs://gnomad-resources/constraint/hail-0.2'

# Unprocessed files
fasta_path = "{}/reference/Homo_sapiens_assembly19.fasta".format(root)
gerp_annotations_path = 'gs://annotationdb/hail-0.2/hail-tables/GRCh37/gerp.scores.GRCh37.ht'  # Gerp is in here as S
regional_variation_raw_path = 'gs://gnomad-resources/constraint/source/whole_genome_regional_variation_in_mutation_rate.50kb_bins.txt'

# Mu exploratory analyses
all_possible_summary_pickle = 'gs://konradk/tmp/all_possible_counts_by_context.pckl'
all_possible_summary_unfiltered_pickle = 'gs://gnomad-resources/constraint/hail-0.2/exploratory/all_possible_counts_by_context_unfiltered.pckl'

# Input datasets
context_mt_path = 'gs://gnomad-resources/context/hail-0.2/Homo_sapiens_assembly19.fasta.snps_only.mt'
split_context_mt_path = 'gs://gnomad-resources/context/hail-0.2/Homo_sapiens_assembly19.fasta.snps_only.split.mt'
# processed_genomes_ht_path = f'{root}/genomes_processed.ht'
processed_exomes_ht_path = f'{root}/exomes_processed.ht'
processed_genomes_ht_path = f'{root}/genomes_processed.mt'
# processed_exomes_ht_path = f'{root}/exomes_processed.mt'

# mutation_rate_ht_path = '{}/standard/mutation_rate.ht'.format(root)
mutation_rate_ht_path = f'{root}/exploratory/mutation_rate_methylation_nocommon.ht'
# po_coverage_ht_path = '{}/standard/prop_observed_by_coverage.ht'.format(root)
po_coverage_ht_path = f'{root}/exploratory/prop_observed_by_coverage.ht'
po_ht_path = '{}/standard/prop_observed.ht'.format(root)
po_syn_ht_path = '{}/standard/prop_observed_syn.ht'.format(root)

po_coverage_x_ht_path = '{}/standard/prop_observed_by_coverage_x.ht'.format(root)
po_x_ht_path = '{}/standard/prop_observed_x.ht'.format(root)
po_syn_x_ht_path = '{}/standard/prop_observed_syn_x.ht'.format(root)

HIGH_COVERAGE_CUTOFF = 0.9
VARIANT_TYPES_FOR_MODEL = ('ACG', 'TCG', 'CCG', 'GCG', 'non-CpG')


# Data loads
def get_old_mu_data() -> hl.Table:
    old_mu_data = hl.import_table('gs://gnomad-resources/constraint/source/fordist_1KG_mutation_rate_table.txt',
                                  delimiter=' ', impute=True)
    return old_mu_data.transmute(context=old_mu_data['from'], ref=old_mu_data['from'][1],
                                 alt=old_mu_data.to[1]).key_by('context', 'ref', 'alt')


def load_regional_data() -> hl.Table:
    ht = hl.import_table(regional_variation_raw_path, impute=True)
    return ht.transmute(filters=ht['filter'],
                        interval=hl.locus_interval(ht.chr.replace('chr', ''), ht.start + 1, ht.end, includes_end=True)
                        ).key_by('interval')


def load_all_possible_summary(filtered: bool = True) -> Dict[hl.Struct, int]:
    fname = all_possible_summary_pickle if filtered else all_possible_summary_unfiltered_pickle
    with hl.hadoop_open(fname, 'rb') as f:
        return pickle.load(f)


# Pre-process
def export_fasta(hc) -> None:
    # Done with an 0.1 jar
    hc.read('gs://gnomad-resources/Homo_sapiens_assembly19.fasta.snps_only.vep.vds').export_variants('gs://gnomad-resources/context/source/Homo_sapiens_assembly19.fasta.snps_only.vep.txt.bgz', 'v, va.context, va.vep', types=True, parallel=True)


def import_fasta() -> None:
    # Works with --jar gs://hail-common/builds/devel/jars/hail-devel-47de006f2c62-Spark-2.2.0.jar
    ht = hl.import_table('gs://gnomad-resources/context/source/Homo_sapiens_assembly19.fasta.snps_only.vep.txt.bgz/*', no_header=True,
                         types={'f0': hl.tstr,
                                'f1': hl.tstr,
                                'f2': hl.tstruct(
                                    assembly_name=hl.tstr,
                                    allele_string=hl.tstr,
                                    ancestral=hl.tstr,
                                    colocated_variants=hl.tarray(
                                        hl.tstruct(aa_allele=hl.tstr, aa_maf=hl.tfloat, afr_allele=hl.tstr, afr_maf=hl.tfloat, allele_string=hl.tstr, amr_allele=hl.tstr, amr_maf=hl.tfloat, clin_sig=hl.tarray(hl.tstr),
                                                   end=hl.tint, eas_allele=hl.tstr, eas_maf=hl.tfloat, ea_allele=hl.tstr, ea_maf=hl.tfloat, eur_allele=hl.tstr, eur_maf=hl.tfloat, exac_adj_allele=hl.tstr, exac_adj_maf=hl.tfloat,
                                                   exac_allele=hl.tstr, exac_afr_allele=hl.tstr, exac_afr_maf=hl.tfloat, exac_amr_allele=hl.tstr, exac_amr_maf=hl.tfloat, exac_eas_allele=hl.tstr, exac_eas_maf=hl.tfloat,
                                                   exac_fin_allele=hl.tstr, exac_fin_maf=hl.tfloat, exac_maf=hl.tfloat, exac_nfe_allele=hl.tstr, exac_nfe_maf=hl.tfloat, exac_oth_allele=hl.tstr, exac_oth_maf=hl.tfloat,
                                                   exac_sas_allele=hl.tstr, exac_sas_maf=hl.tfloat, id=hl.tstr, minor_allele=hl.tstr, minor_allele_freq=hl.tfloat, phenotype_or_disease=hl.tint, pubmed=hl.tarray(hl.tint),
                                                   sas_allele=hl.tstr, sas_maf=hl.tfloat, somatic=hl.tint, start=hl.tint, strand=hl.tint),
                                    ),
                                    context=hl.tstr, end=hl.tint, id=hl.tstr, input=hl.tstr,
                                    intergenic_consequences=hl.tarray(
                                        hl.tstruct(allele_num=hl.tint, consequence_terms=hl.tarray(hl.tstr), impact=hl.tstr, minimised=hl.tint, variant_allele=hl.tstr)),
                                    most_severe_consequence=hl.tstr, motif_feature_consequences=hl.tarray(
                                        hl.tstruct(allele_num=hl.tint, consequence_terms=hl.tarray(hl.tstr), high_inf_pos=hl.tstr,impact=hl.tstr, minimised=hl.tint, motif_feature_id=hl.tstr, motif_name=hl.tstr, motif_pos=hl.tint, motif_score_change=hl.tfloat, strand=hl.tint, variant_allele=hl.tstr)),
                                    regulatory_feature_consequences=hl.tarray(hl.tstruct(
                                        allele_num=hl.tint, biotype=hl.tstr, consequence_terms=hl.tarray(hl.tstr), impact=hl.tstr, minimised=hl.tint, regulatory_feature_id=hl.tstr, variant_allele=hl.tstr)
                                    ), seq_region_name=hl.tstr, start=hl.tint, strand=hl.tint, transcript_consequences=hl.tarray(hl.tstruct(
                                        allele_num=hl.tint, amino_acids=hl.tstr, biotype=hl.tstr, canonical=hl.tint, ccds=hl.tstr, cdna_start=hl.tint, cdna_end=hl.tint, cds_end=hl.tint, cds_start=hl.tint,
                                        codons=hl.tstr, consequence_terms=hl.tarray(hl.tstr), distance=hl.tint, domains=hl.tarray(hl.tstruct(db=hl.tstr, name=hl.tstr)), exon=hl.tstr,
                                        gene_id=hl.tstr, gene_pheno=hl.tint, gene_symbol=hl.tstr, gene_symbol_source=hl.tstr, hgnc_id=hl.tint, hgvsc=hl.tstr, hgvsp=hl.tstr, hgvs_offset=hl.tint,
                                        impact=hl.tstr, intron=hl.tstr, lof=hl.tstr, lof_flags=hl.tstr, lof_filter=hl.tstr, lof_info=hl.tstr, minimised=hl.tint, polyphen_prediction=hl.tstr, polyphen_score=hl.tfloat,
                                        protein_end=hl.tint, protein_start=hl.tint, protein_id=hl.tstr, sift_prediction=hl.tstr, sift_score=hl.tfloat, strand=hl.tint, swissprot=hl.tstr, transcript_id=hl.tstr, trembl=hl.tstr, uniparc=hl.tstr, variant_allele=hl.tstr
                                    )), variant_class=hl.tstr
                                )}
    ).rename({'f0': 'v', 'f1': 'context', 'f2': 'vep'})
    ht.transmute(**hl.parse_variant(ht.v)).key_by('locus', 'alleles').write(context_mt_path)
    # ht = ht.transmute(**hl.parse_variant(ht.v)).key_by('locus', 'alleles')
    # hl.MatrixTable.from_rows_table(ht, partition_key='locus').write(context_mt_path)


def split_context_mt(raw_context_mt_path: str, coverage_ht_paths: Dict[str, str], methylation_ht_path: str,
                     split_context_mt_path: str, overwrite: bool = False) -> None:
    raw_context_mt = split_multi_dynamic(hl.read_matrix_table(raw_context_mt_path))

    coverage_hts = {loc: hl.read_table(coverage_ht_path) for loc, coverage_ht_path in coverage_ht_paths.items()}
    coverage_hts = {loc: coverage_ht.drop('#chrom', 'pos') if '#chrom' in list(coverage_ht.row) else coverage_ht
                    for loc, coverage_ht in coverage_hts.items()}
    methylation_ht = hl.read_table(methylation_ht_path)
    gerp_ht = hl.read_table(gerp_annotations_path)

    raw_context_mt = raw_context_mt.annotate_rows(
        methylation=methylation_ht[raw_context_mt.locus],
        coverage=hl.struct(**{loc: coverage_ht[raw_context_mt.locus] for loc, coverage_ht in coverage_hts.items()}),
        gerp=gerp_ht[raw_context_mt.locus].S)

    raw_context_mt.write(split_context_mt_path, overwrite)


def pre_process_data(ht_path: str, split_context_mt_path: str,
                     output_ht_path: str, overwrite: bool = False) -> None:
    # ht = hl.read_table(ht_path)
    # mt = hl.MatrixTable.from_rows_table(ht, partition_key='locus')
    # mt.write(ht_path.replace('.ht', '.mt'))
    mt = hl.read_matrix_table(ht_path.replace('.ht', '.mt'))
    context_mt = hl.read_matrix_table(split_context_mt_path).drop('a_index', 'was_split')
    context_mt = context_mt.annotate_rows(vep=context_mt.vep.drop('colocated_variants'))
    context_mt = hl.filter_intervals(context_mt, [hl.parse_locus_interval('1-22')])
    mt.annotate_rows(**context_mt[mt.row_key, :]).write(output_ht_path, overwrite)


def prepare_ht(ht, trimer: bool = False):
    if trimer:
        ht = trimer_from_heptamer(ht)
    str_len = 3 if trimer else 7

    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht))
        ht = ht.annotate(methylation_level=hl.cond(ht.cpg, hl.int(ht.methylation.MEAN * 10), -1),
                         # hl.or_missing(ht.cpg, hl.int(ht.methylation.MEAN * 10)),
                         csq=ht.vep.most_severe_consequence,
                         exome_coverage=hl.int(ht.coverage.exomes.mean * 10))
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht))
        ht = ht.annotate_rows(methylation_level=hl.cond(ht.cpg, hl.int(ht.methylation.MEAN * 10), -1),
                              # hl.or_missing(ht.cpg, hl.int(ht.methylation.MEAN * 10)),
                              csq=ht.vep.most_severe_consequence,
                              exome_coverage=hl.int(ht.coverage.exomes.mean * 10))

    return ht


# Explore mu
def calculate_mu_by_downsampling(genome_ht: hl.Table, raw_context_ht: hl.MatrixTable,
                                 recalculate_all_possible_summary: bool = True,
                                 recalculate_all_possible_summary_unfiltered: bool = False,
                                 omit_methylation: bool = False, count_singletons: bool = False,
                                 remove_common: bool = True,
                                 grouping_variables: Optional[Tuple[str]] = (), summary_file: str = None
                                 ) -> hl.Table:

    context_ht = filter_to_autosomes(remove_coverage_outliers(raw_context_ht))
    genome_ht = remove_coverage_outliers(genome_ht)

    context_ht = context_ht.select_rows('context', 'ref', 'alt', 'methylation_level', 'gerp', 'csq')
    genome_ht = genome_ht.select_rows('context', 'ref', 'alt', 'methylation_level', 'gerp', 'csq', 'freq')

    if grouping_variables:
        context_ht = context_ht.filter_rows(hl.all(lambda x: x, [hl.is_defined(context_ht[x]) for x in grouping_variables]))
        genome_ht = genome_ht.filter_rows(hl.all(lambda x: x, [hl.is_defined(genome_ht[x]) for x in grouping_variables]))
    else:
        context_ht = filter_for_mu(context_ht)
        genome_ht = filter_for_mu(genome_ht)

    if remove_common:
        ac_cutoff = 5
        downsampling_level = 1000
        # Removes anything with AC > ac_cutoff from denominator
        genome_freq = genome_ht[context_ht.row_key, :].freq
        context_ht = context_ht.filter_rows(
            hl.is_missing(genome_freq) |
            hl.any(lambda f: (f.AC[1] <= ac_cutoff) & (f.meta.get('downsampling') == str(downsampling_level)) &
                             (f.meta.get('pop') == 'global') & (f.meta.get('group') == 'adj') & (f.meta.size() == 3),
                   genome_freq)
        ).rows()
        # Keep only AC < ac_cutoff + 1 in numerator
        genome_ht = filter_by_frequency(genome_ht, allele_count=ac_cutoff + 1, downsampling=downsampling_level,
                                        direction='below', adj=True).rows()

    if not summary_file: summary_file = all_possible_summary_pickle
    all_possible_dtype = count_variants(context_ht, omit_methylation=omit_methylation,
                                        additional_grouping=grouping_variables, return_type_only=True)
    if recalculate_all_possible_summary:
        all_possible = count_variants(context_ht, omit_methylation=omit_methylation,
                                      additional_grouping=grouping_variables).variant_count
        with hl.hadoop_open(summary_file, 'wb') as f:
            pickle.dump(all_possible, f)
    with hl.hadoop_open(summary_file, 'rb') as f:
        all_possible = pickle.load(f)

    if recalculate_all_possible_summary_unfiltered:
        all_possible_unfiltered = count_variants(raw_context_ht.rows(), omit_methylation=omit_methylation).variant_count
        all_possible_unfiltered = {x: y for x, y in list(all_possible_unfiltered.items()) if x.context is not None}
        with hl.hadoop_open(all_possible_summary_unfiltered_pickle, 'wb') as f:
            pickle.dump(all_possible_unfiltered, f)
    all_possible_unfiltered = load_all_possible_summary(filtered=False)

    genome_ht = count_variants(genome_ht, count_downsamplings=['global', 'nfe', 'afr'],
                               count_singletons=count_singletons,
                               additional_grouping=grouping_variables, omit_methylation=omit_methylation)

    old_mu_data = get_old_mu_data()
    ht = genome_ht.annotate(possible_variants=hl.literal(all_possible, dtype=all_possible_dtype)[genome_ht.key],
                            # possible_variants_unfiltered=hl.literal(all_possible_unfiltered, dtype=all_possible_dtype)[genome_ht.key],
                            old_mu_snp=old_mu_data[hl.struct(
                                context=genome_ht.context, ref=genome_ht.ref, alt=genome_ht.alt)].mu_snp)
    ht = ht.persist()
    total_bases = ht.aggregate(hl.agg.sum(ht.possible_variants)) // 3
    # total_bases_unfiltered = ht.aggregate(hl.agg.sum(ht.possible_variants_unfiltered)) // 3
    # total_mu = ht.aggregate(hl.agg.sum(ht.old_mu_snp * ht.possible_variants_unfiltered) / total_bases_unfiltered)
    total_mu = 1.1926e-08  # (should be 1.2e-08)

    correction_factors = ht.aggregate(total_mu / (hl.agg.array_sum(ht.downsampling_counts_global) / total_bases))
    correction_factors_nfe = ht.aggregate(total_mu / (hl.agg.array_sum(ht.downsampling_counts_nfe) / total_bases))
    correction_factors_afr = ht.aggregate(total_mu / (hl.agg.array_sum(ht.downsampling_counts_afr) / total_bases))

    ht = annotate_variant_types(ht.annotate(downsamplings_frac_observed=ht.downsampling_counts_global / ht.possible_variants,
                                            downsamplings_mu_snp=hl.literal(correction_factors) * ht.downsampling_counts_global / ht.possible_variants,
                                            downsamplings_mu_nfe=hl.literal(correction_factors_nfe) * ht.downsampling_counts_nfe / ht.possible_variants,
                                            downsamplings_mu_afr=hl.literal(correction_factors_afr) * ht.downsampling_counts_afr / ht.possible_variants))
    downsamplings = ht.downsamplings.collect()[0]
    index_1kg = downsamplings.index(1000)
    return ht.annotate(mu_snp=ht.downsamplings_mu_snp[index_1kg],
                       mu_snp_nfe=ht.downsamplings_mu_nfe[index_1kg],
                       mu_snp_afr=ht.downsamplings_mu_afr[index_1kg])


def get_proportion_observed_by_coverage(exome_ht: hl.Table, context_ht: hl.MatrixTable, mutation_ht: hl.Table,
                                        recompute_possible: bool = False) -> hl.Table:

    context_ht = fast_filter_vep(context_ht).select_rows('context', 'ref', 'alt', 'methylation_level', 'exome_coverage')
    context_ht = context_ht.filter_rows(hl.is_defined(context_ht.exome_coverage))

    exome_ht = fast_filter_vep(exome_ht).select('context', 'ref', 'alt', 'methylation_level', 'exome_coverage', 'freq')

    grouping = ('exome_coverage',)

    possible_file = 'gs://konradk/tmp/possible.pckl'
    if recompute_possible:
        possible_variants = count_variants(context_ht.rows(), additional_grouping=grouping)
        with hl.hadoop_open(possible_file, 'wb') as f:
            pickle.dump(possible_variants, f)

    with hl.hadoop_open(possible_file, 'rb') as f:
        possible_variants = pickle.load(f)

    variant_counts_dtype = count_variants(context_ht.rows(), additional_grouping=grouping, return_type_only=True)
    pprint(variant_counts_dtype)
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=100, count_downsamplings=('global',))
    ht = ht.annotate(possible_variants=hl.literal(possible_variants.variant_count, dtype=variant_counts_dtype)[ht.key],
                     mu_snp=mutation_ht[hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt, methylation_level=ht.methylation_level)].mu_snp)
    return ht


# Model building
def calculate_coverage_model(coverage_ht: hl.Table) -> (float, float):
    """
    Calibrates coverage model (returns slope and intercept)
    """
    coverage_pd = coverage_ht.to_pandas()
    lm = smf.ols(formula='log_mean_scaled_proportion_observed ~ log_coverage', data=coverage_pd).fit()
    slope = lm.params['log_coverage']
    intercept = lm.params['Intercept']
    return slope, intercept


def build_plateau_models(ht: hl.Table) -> dict:
    """
    Calibrates high coverage model (returns slope and intercept)
    """
    output = {}
    for variant_type_model in VARIANT_TYPES_FOR_MODEL:
        high_coverage_pd = ht.filter(ht.variant_type_model == variant_type_model).to_pandas()
        lm = smf.ols(formula='high_coverage_proportion_observed ~ mutation_rate', data=high_coverage_pd).fit()
        output[variant_type_model] = dict(lm.params)
    return output


# Misc
def maps(ht: hl.Table, mutation_ht: hl.Table, vep_root: str = 'vep') -> hl.Table:
    ht = ht.annotate(csq=ht[vep_root]['most_severe_consequence'])
    ht = count_variants(ht, count_singletons=True, additional_grouping=('csq', ))
    ht = ht.annotate(mu=mutation_ht[ht.key].mu_snp,
                     ps=ht.singleton_count / ht.variant_count)
    syn_ps_pd = ht.filter(ht.csq == 'synonymous_variant').to_pandas()

    lm = smf.ols(formula='ps ~ mu', data=syn_ps_pd).fit()
    slope = lm.params['mu']
    intercept = lm.params['Intercept']
    ht = ht.annotate(expected_singletons=(ht.mu * slope + intercept) * ht.n_variants)

    agg_ht = (ht.group_by('csq')
              .aggregate(n_singletons=hl.agg.sum(ht.singleton_count),
                         expected_singletons=hl.agg.sum(ht.expected_singletons),
                         n_variants=hl.agg.sum(ht.variant_count)))
    agg_ht = agg_ht.annotate(ps=agg_ht.singleton_count / agg_ht.variant_count,
                             maps=(agg_ht.singleton_count - agg_ht.expected_singletons) / agg_ht.variant_count)
    agg_ht = agg_ht.annotate(sem_ps=(agg_ht.ps * (1 - agg_ht.ps) / agg_ht.variant_count) ** 0.5)
    return agg_ht


# Plotting
def old_new_compare(source, axis_type='log'):
    p1 = figure(title="Mutation rate comparison", x_axis_type=axis_type, y_axis_type=axis_type, tools=TOOLS)
    p1.xaxis.axis_label = 'Old mutation rate'
    p1.yaxis.axis_label = 'New mutation rate'

    p1.scatter(x='old_mu_snp', y='mu_snp', fill_color='colors', line_color='colors', legend='variant_type',
               source=source)

    p1.select_one(HoverTool).tooltips = [(x, f'@{x}') for x in
                                         ('context', 'ref', 'alt', 'methylation_level', 'old_mu_snp', 'mu_snp') if
                                         x in list(source.data)]
    # p1.ray(x=[1e-9], y=[1e-9], length=0, angle=[45], angle_units="deg", color="#FB8072", line_width=2)
    p1.legend.location = "top_left"
    return p1