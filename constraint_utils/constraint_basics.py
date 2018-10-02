__author__ = 'konradk'

from .generic import *
import pickle
import copy

root = 'gs://gnomad-resources/constraint/hail-0.2'

# Unprocessed files
fasta_path = "{}/reference/Homo_sapiens_assembly19.fasta".format(root)
gerp_annotations_path = 'gs://gnomad-resources/annotations/gerp.scores.GRCh37.ht'  # Gerp is in here as S
regional_variation_raw_path = 'gs://gnomad-resources/constraint/source/whole_genome_regional_variation_in_mutation_rate.50kb_bins.txt'

# Mu exploratory analyses
all_possible_summary_pickle = 'gs://konradk/tmp/all_possible_counts_by_context.pckl'
all_possible_summary_unfiltered_pickle = 'gs://gnomad-resources/constraint/hail-0.2/exploratory/all_possible_counts_by_context_unfiltered.pckl'

# Input datasets
unsplit_context_mt_path = 'gs://gnomad-resources/context/hail-0.2/Homo_sapiens_assembly19.fasta.snps_only_20180907.mt'
processed_genomes_ht_path = f'{root}/genomes_processed.mt'
processed_exomes_ht_path = f'{root}/exomes_processed.mt'

# location = 'standard'
location = 'weighted'
mutation_rate_ht_path = f'{root}/standard/mutation_rate_methylation_bins.ht'
po_coverage_ht_path = f'{root}/standard/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
po_ht_path = f'{root}/{location}/prop_observed.ht'
constraint_ht_path = f'{root}/{location}/constraint.ht'

HIGH_COVERAGE_CUTOFF = 40
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
    ht.transmute(**hl.parse_variant(ht.v)).key_by('locus', 'alleles').write(unsplit_context_mt_path)
    # ht = ht.transmute(**hl.parse_variant(ht.v)).key_by('locus', 'alleles')
    # hl.MatrixTable.from_rows_table(ht).write(unsplit_context_mt_path)


def split_context_mt(raw_context_mt_path: str, coverage_ht_paths: Dict[str, str], methylation_ht_path: str,
                     split_context_mt_path: str, overwrite: bool = False) -> None:
    # raw_context_mt = hl.split_multi_hts(hl.read_matrix_table(raw_context_mt_path))
    # raw_context_mt.write(f'{raw_context_mt_path}.split.mt', overwrite)
    raw_context_mt = hl.read_matrix_table(f'{raw_context_mt_path}.split.mt')

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


def vep_context_mt():
    mt_full = hl.read_matrix_table('gs://gnomad-resources/context/hail-0.2/Homo_sapiens_assembly19.fasta.snps_only.mt')
    chunks = ('1-2', '3-4', '5-6',  '7-9', '10-12', '13-15', '16-18', '19-22', 'X-Y')
    for i in chunks:
        print(i)
        mt = hl.filter_intervals(mt_full, [hl.parse_locus_interval(i)])
        hl.vep(mt, vep_config).write(f'gs://gnomad-resources/context/hail-0.2/parts/Homo_sapiens_assembly19.fasta.snps_only.{i}.mt', overwrite=True, stage_locally=True)
    print('Done! Combining...')
    hts = [hl.read_matrix_table(f'gs://gnomad-resources/context/hail-0.2/parts/Homo_sapiens_assembly19.fasta.snps_only.{i}.mt') for i in chunks]
    ht = hts[0].union_rows(*hts[1:])
    ht.write(unsplit_context_mt_path, True)


def pre_process_data(ht_path: str, rf_path: str, split_context_mt_path: str,
                     output_ht_path: str, overwrite: bool = False) -> None:
    # ht = hl.read_table(ht_path)
    # mt = hl.MatrixTable.from_rows_table(ht, partition_key='locus')
    # mt.write(ht_path.replace('.ht', '.mt'), overwrite)
    mt = hl.read_matrix_table(ht_path.replace('.ht', '.mt'))
    context_mt = hl.read_matrix_table(split_context_mt_path).drop('a_index', 'was_split')
    context_mt = context_mt.annotate_rows(vep=context_mt.vep.drop('colocated_variants'))
    rf_ht = hl.read_table(rf_path)
    mt.annotate_rows(**context_mt[mt.row_key, :],
                     pass_filters=(hl.len(rf_ht[mt.row_key].filters) == 0) & (mt.freq[0].AC[1] > 0)
                     ).write(output_ht_path, overwrite)


def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    if trimer:
        ht = trimer_from_heptamer(ht)
    str_len = 3 if trimer else 7

    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    annotation = {
        'methylation_level': hl.case().when(
            ht.cpg & (ht.methylation.MEAN > 0.6), 2
        ).when(
            ht.cpg & (ht.methylation.MEAN > 0.2), 1
        ).default(0)
        # 'methylation_level': hl.cond(ht.cpg, hl.int(ht.methylation.MEAN * 10), -1),
    }
    if annotate_coverage:
        # coverage_binning = 100
        annotation['exome_coverage'] = ht.coverage.exomes.median
        # hl.int(ht.coverage.exomes.over_20 * coverage_binning) / coverage_binning
    return ht.annotate(**annotation) if isinstance(ht, hl.Table) else ht.annotate_rows(**annotation)


def filter_to_autosomes_par(ht: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    return ht.filter(ht.locus.in_autosome_or_par())


def annotate_with_mu(ht: hl.Table, mutation_ht: hl.Table, output_loc: str = 'mu_snp',
                     keys: Tuple[str] = ('context', 'ref', 'alt', 'methylation_level')) -> hl.Table:
    mu = hl.literal(mutation_ht.aggregate(hl.dict(hl.agg.collect(
        (hl.struct(**{k: mutation_ht[k] for k in keys}), mutation_ht.mu_snp)))))
    mu = mu.get(hl.struct(**{k: ht[k] for k in keys}))
    return ht.annotate(**{output_loc: hl.case().when(hl.is_defined(mu), mu).or_error('Missing mu')})


def calculate_mu_by_downsampling(genome_ht: hl.Table, raw_context_ht: hl.MatrixTable,
                                 recalculate_all_possible_summary: bool = True,
                                 recalculate_all_possible_summary_unfiltered: bool = False,
                                 omit_methylation: bool = False, count_singletons: bool = False,
                                 remove_common_downsampled: bool = True, remove_common_ordinary: bool = False,
                                 grouping_variables: Optional[Tuple[str]] = (), summary_file: str = None
                                 ) -> hl.Table:

    context_ht = filter_to_autosomes(remove_coverage_outliers(raw_context_ht))
    genome_ht = filter_to_autosomes(remove_coverage_outliers(genome_ht))

    if grouping_variables:
        context_ht = context_ht.filter_rows(hl.all(lambda x: x, [hl.is_defined(context_ht[x]) for x in grouping_variables]))
        genome_ht = genome_ht.filter_rows(hl.all(lambda x: x, [hl.is_defined(genome_ht[x]) for x in grouping_variables]))
    else:
        context_ht = filter_for_mu(context_ht)
        genome_ht = filter_for_mu(genome_ht)

    context_ht = context_ht.select_rows('context', 'ref', 'alt', 'methylation_level', *grouping_variables)
    genome_ht = genome_ht.select_rows('context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping_variables)

    genome_join = genome_ht[context_ht.row_key, :]
    if remove_common_downsampled:
        ac_cutoff = 5
        downsampling_level = 1000
        # In the denominator, only keep variants not in the genome dataset, or with AC <= ac_cutoff and passing filters
        context_ht = context_ht.filter_rows(
            hl.is_missing(genome_join) |
            (hl.any(lambda f: (f.AC[1] <= ac_cutoff) & (f.meta.get('downsampling') == str(downsampling_level)) &
                              (f.meta.get('pop') == 'global') & (f.meta.get('group') == 'adj') & (f.meta.size() == 3),
                    genome_join.freq) & genome_join.pass_filters)
        ).rows()
        # Keep only AC < ac_cutoff + 1 in numerator
        genome_ht = filter_by_frequency(genome_ht, 'below', allele_count=ac_cutoff + 1,
                                        downsampling=downsampling_level, adj=True)
    elif remove_common_ordinary:
        af_cutoff = 0.001
        context_ht = context_ht.filter_rows(
            hl.is_missing(genome_join) | (
                    # hl.any(lambda f: (f.AF[1] <= af_cutoff) & (f.meta.get('group') == 'adj') & (f.meta.size() == 1), genome_join.freq) &
                    (genome_join.freq[0].AF[1] <= af_cutoff) &
                    genome_join.pass_filters
            )
        ).rows()
        # Keep only AF < af_cutoff in numerator
        genome_ht = filter_by_frequency(genome_ht, 'below', frequency=af_cutoff, adj=True)
    else:
        context_ht = context_ht.filter_rows(
            hl.is_missing(genome_join) | genome_join.pass_filters
        ).rows()
    genome_ht = genome_ht.filter_rows(genome_ht.pass_filters).rows()

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
    # Uncomment this line and next few commented lines to back-calculate total_mu from the old mutation rate dataset
    # all_possible_unfiltered = load_all_possible_summary(filtered=False)

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
    total_mu = 1.2e-08

    correction_factors = ht.aggregate(total_mu / (hl.agg.array_sum(ht.downsampling_counts_global) / total_bases))
    correction_factors_nfe = ht.aggregate(total_mu / (hl.agg.array_sum(ht.downsampling_counts_nfe) / total_bases))
    correction_factors_afr = ht.aggregate(total_mu / (hl.agg.array_sum(ht.downsampling_counts_afr) / total_bases))

    ht = annotate_variant_types(ht.annotate(downsamplings_frac_observed=ht.downsampling_counts_global / ht.possible_variants,
                                            downsamplings_mu_snp=hl.literal(correction_factors) * ht.downsampling_counts_global / ht.possible_variants,
                                            downsamplings_mu_nfe=hl.literal(correction_factors_nfe) * ht.downsampling_counts_nfe / ht.possible_variants,
                                            downsamplings_mu_afr=hl.literal(correction_factors_afr) * ht.downsampling_counts_afr / ht.possible_variants))
    downsamplings = ht.downsamplings.collect()[0]
    index_1kg = downsamplings.index(1000)
    return ht.annotate(proportion_observed_1kg=ht.downsampling_counts_global[index_1kg] / ht.possible_variants,
                       proportion_observed=ht.variant_count / ht.possible_variants,
                       mu_snp=ht.downsamplings_mu_snp[index_1kg],
                       mu_snp_nfe=ht.downsamplings_mu_nfe[index_1kg],
                       mu_snp_afr=ht.downsamplings_mu_afr[index_1kg])


def get_proportion_observed_by_coverage(exome_ht: hl.MatrixTable, context_ht: hl.MatrixTable, mutation_ht: hl.Table,
                                        recompute_possible: bool = False, remove_from_denominator: bool = True,
                                        remove_filtered: bool = True, remove_ac0: bool = False) -> hl.Table:

    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    context_ht = add_most_severe_csq_to_tc_within_ht(context_ht)

    context_ht = fast_filter_vep(context_ht).select_rows('context', 'ref', 'alt', 'methylation_level', 'exome_coverage')
    context_ht = context_ht.filter_rows(hl.is_defined(context_ht.exome_coverage))

    exome_ht = fast_filter_vep(exome_ht).select_rows('context', 'ref', 'alt', 'methylation_level', 'exome_coverage', 'freq', 'pass_filters')

    grouping = ('exome_coverage',)
    af_cutoff = 0.001

    exome_join = exome_ht[context_ht.row_key, :]
    keep_criteria = True
    if remove_from_denominator:
        keep_criteria &= hl.any(lambda f: (f.AF[1] <= af_cutoff) &  # (f.AF[1] > 0) &
                                          (f.meta.get('group') == 'adj') & (f.meta.size() == 1),
                                exome_join.freq)
    if remove_filtered:
        keep_criteria &= exome_join.pass_filters
    context_ht = context_ht.filter_rows(hl.is_missing(exome_join) | keep_criteria).rows()

    if remove_filtered:
        exome_ht = exome_ht.filter_rows(exome_ht.pass_filters)

    if remove_ac0:
        exome_ht = filter_by_frequency(exome_ht, 'above', allele_count=0)
    exome_ht = filter_by_frequency(exome_ht, 'below', frequency=af_cutoff).rows()

    possible_file = 'gs://konradk/tmp/possible_coverage.ht'
    if recompute_possible:
        possible_ht = count_variants(context_ht, additional_grouping=grouping, force_grouping=True)
        possible_ht = annotate_with_mu(possible_ht, mutation_ht)
        possible_ht.transmute(possible_variants=possible_ht.variant_count).write(possible_file, True)

    possible_ht = hl.read_table(possible_file)
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=100, count_downsamplings=('global',))
    ht = ht.join(possible_ht, 'outer')
    return ht


def build_models(coverage_ht: hl.Table, trimers: bool = False, weighted: bool = False) -> Tuple[Tuple[float, float], Dict[str, Tuple[float, float]]]:
    keys = ['context', 'ref', 'alt', 'methylation_level', 'mu_snp']

    all_high_coverage_ht = coverage_ht.filter(coverage_ht.exome_coverage >= HIGH_COVERAGE_CUTOFF)
    high_coverage_ht = all_high_coverage_ht.group_by(*keys).aggregate(
        observed_variants=hl.agg.sum(all_high_coverage_ht.variant_count),
        possible_variants=hl.agg.sum(all_high_coverage_ht.possible_variants))

    plateau_models = build_plateau_models(
        annotate_variant_types(
            high_coverage_ht.annotate(
                high_coverage_proportion_observed=high_coverage_ht.observed_variants / high_coverage_ht.possible_variants
            ), not trimers),
        weighted=weighted
    )

    high_coverage_scale_factor = all_high_coverage_ht.aggregate(
        hl.agg.sum(all_high_coverage_ht.variant_count) /
        hl.agg.sum(all_high_coverage_ht.possible_variants * all_high_coverage_ht.mu_snp))

    all_low_coverage_ht = coverage_ht.filter((coverage_ht.exome_coverage < HIGH_COVERAGE_CUTOFF) &
                                             (coverage_ht.exome_coverage > 0))

    low_coverage_ht = all_low_coverage_ht.group_by(log_coverage=hl.log10(all_low_coverage_ht.exome_coverage)).aggregate(
        low_coverage_obs_exp=hl.agg.sum(all_low_coverage_ht.variant_count) /
                             (high_coverage_scale_factor * hl.agg.sum(all_low_coverage_ht.possible_variants * all_low_coverage_ht.mu_snp)))
    coverage_model = build_coverage_model(low_coverage_ht)
    # TODO: consider weighting here as well

    return coverage_model, plateau_models


def add_most_severe_csq_to_tc_within_ht(t):
    annotation = t.vep.annotate(transcript_consequences=t.vep.transcript_consequences.map(
        add_most_severe_consequence_to_consequence))
    return t.annotate_rows(vep=annotation) if isinstance(t, hl.MatrixTable) else t.annotate(vep=annotation)


def take_one_annotation_from_tc_within_ht(t):
    annotation = t.vep.annotate(transcript_consequences=t.vep.transcript_consequences[0])
    return t.annotate_rows(vep=annotation) if isinstance(t, hl.MatrixTable) else t.annotate(vep=annotation)


def get_proportion_observed(exome_ht: hl.MatrixTable, context_ht: hl.MatrixTable, mutation_ht: hl.Table,
                            plateau_models: Dict[str, Tuple[float, float]], coverage_model: Tuple[float, float],
                            recompute_possible: bool = False,
                            remove_from_denominator: bool = True, custom_model: str = None) -> hl.Table:

    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    context_ht = add_most_severe_csq_to_tc_within_ht(context_ht)

    if custom_model == 'syn_canonical':
        context_ht = take_one_annotation_from_tc_within_ht(fast_filter_vep(context_ht))
        exome_ht = take_one_annotation_from_tc_within_ht(fast_filter_vep(exome_ht))
    elif custom_model == 'worst_csq':
        context_ht = process_consequences(context_ht)
        context_ht = context_ht.explode_rows(context_ht.vep.worst_csq_by_gene)
        exome_ht = process_consequences(exome_ht)
        exome_ht = exome_ht.explode_rows(exome_ht.vep.worst_csq_by_gene)
    else:
        context_ht = context_ht.explode_rows(context_ht.vep.transcript_consequences)
        exome_ht = exome_ht.explode_rows(exome_ht.vep.transcript_consequences)

    context_ht, _ = annotate_constraint_groupings(context_ht, custom_model=custom_model)
    exome_ht, grouping = annotate_constraint_groupings(exome_ht, custom_model=custom_model)

    context_ht = context_ht.filter_rows(hl.is_defined(context_ht.exome_coverage)).select_rows(
        'context', 'ref', 'alt', 'methylation_level', *grouping)
    exome_ht = exome_ht.select_rows(
        'context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping)

    af_cutoff = 0.001

    exome_join = exome_ht[context_ht.row_key, :]
    if remove_from_denominator:
        context_ht = context_ht.filter_rows(
            hl.is_missing(exome_join) |
            (hl.any(lambda f: (f.AF[1] <= af_cutoff) & (f.AF[1] > 0) & (f.meta.get('group') == 'adj') & (f.meta.size() == 1),
                    exome_join.freq) & exome_join.pass_filters)
        ).rows()
    else:
        context_ht = context_ht.rows()

    exome_ht = filter_by_frequency(exome_ht.filter_rows(exome_ht.coverage > 0), 'above', allele_count=0)
    exome_ht = filter_by_frequency(exome_ht.filter_rows(exome_ht.pass_filters), 'below', frequency=af_cutoff).rows()

    possible_file = f'{root}/{location}/possible_data/possible_transcript_{custom_model}.ht'
    if recompute_possible:
        ht = count_variants(context_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True)
        ht = annotate_with_mu(ht, mutation_ht)
        ht = ht.transmute(possible_variants=ht.variant_count)
        ht = annotate_variant_types(ht.annotate(mu_agg=ht.mu_snp * ht.possible_variants))
        model = hl.literal(plateau_models)[ht.cpg]
        ht = ht.annotate(adjusted_mutation_rate=ht.mu_agg * model[1] + model[0],
                         coverage_correction=hl.case()
                         .when(ht.coverage == 0, 0)
                         .when(ht.coverage >= HIGH_COVERAGE_CUTOFF, 1)
                         .default(coverage_model[1] * hl.log10(ht.coverage) + coverage_model[0]))
        ht = ht.annotate(expected_variants=ht.adjusted_mutation_rate * ht.coverage_correction)
        ht.write(possible_file, True)

    possible_variants_ht = hl.read_table(possible_file)
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True)  # , count_downsamplings=('global',))
    ht = ht.join(possible_variants_ht, 'outer')
    ht.write(f'{root}/{location}/possible_data/all_data_transcript_{custom_model}.ht', True)
    ht = hl.read_table(f'{root}/{location}/possible_data/all_data_transcript_{custom_model}.ht')

    grouping.remove('coverage')
    ht = ht.group_by(*grouping).partition_hint(1000).aggregate(variant_count=hl.agg.sum(ht.variant_count),
                                                               adjusted_mutation_rate=hl.agg.sum(ht.adjusted_mutation_rate),
                                                               possible_variants=hl.agg.sum(ht.possible_variants),
                                                               expected_variants=hl.agg.sum(ht.expected_variants))
    return ht.annotate(obs_exp=ht.variant_count / ht.expected_variants)


def finalize_dataset(po_ht: hl.Table, skip_transcript: bool = False, n_partitions: int = 100) -> hl.Table:
    keys = ['gene']
    if not skip_transcript: keys.extend(['transcript', 'canonical'])
    po_ht = po_ht.repartition(n_partitions).persist()

    # Getting classic LoF annotations (no LOFTEE)
    classic_lof_annotations = hl.literal({'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'})
    lof_ht_classic = po_ht.filter(classic_lof_annotations.contains(po_ht.annotation) &
                                  ((po_ht.modifier == 'HC') | (po_ht.modifier == 'LC')))
    lof_ht_classic = collapse_lof_ht(lof_ht_classic, keys)
    lof_ht_classic = lof_ht_classic.rename({x: f'{x}_classic' for x in list(lof_ht_classic.row_value)})

    # Getting classic LoF annotations (with LOFTEE)
    lof_ht_classic_hc = po_ht.filter(classic_lof_annotations.contains(po_ht.annotation) & (po_ht.modifier == 'HC'))
    lof_ht_classic_hc = collapse_lof_ht(lof_ht_classic_hc, keys)
    lof_ht_classic_hc = lof_ht_classic_hc.rename({x: f'{x}_classic_hc' for x in list(lof_ht_classic_hc.row_value)})

    # TODO: assess diff of this and classic HC counts to view depletions
    # Getting all LoF annotations (LOFTEE HC + new splice variants)
    lof_ht = po_ht.filter(po_ht.modifier == 'HC')
    po_ht = po_ht.drop('possible_variants', 'adjusted_mutation_rate', 'obs_exp')
    lof_ht = collapse_lof_ht(lof_ht, keys)
    lof_ht = lof_ht.annotate(**oe_confidence_interval(lof_ht, lof_ht.obs_lof, lof_ht.exp_lof)[lof_ht.key])

    mis_ht = po_ht.filter(po_ht.annotation == 'missense_variant')
    mis_ht = mis_ht.group_by(*keys).aggregate(obs_mis=hl.agg.sum(mis_ht.variant_count),
                                              exp_mis=hl.agg.sum(mis_ht.expected_variants),
                                              oe_mis=hl.agg.sum(mis_ht.variant_count) / hl.agg.sum(mis_ht.expected_variants)
                                              )

    # TODO: Aggregate to fix XG problem
    pphen_mis_ht = po_ht.filter(po_ht.modifier == 'probably_damaging').key_by(*keys).drop('modifier', 'annotation')
    pphen_mis_ht = pphen_mis_ht.transmute(obs_mis_pphen=pphen_mis_ht.variant_count,
                                          exp_mis_pphen=pphen_mis_ht.expected_variants,
                                          oe_mis_pphen=pphen_mis_ht.variant_count / pphen_mis_ht.expected_variants)

    # TODO: change this to aggregate to get the non-None ones?
    # TODO: Aggregate to fix XG problem
    syn_ht = po_ht.filter((po_ht.annotation == 'synonymous_variant') & (po_ht.modifier == 'None')).key_by(*keys).drop('modifier', 'annotation')
    syn_ht = syn_ht.transmute(obs_syn=syn_ht.variant_count, exp_syn=syn_ht.expected_variants,
                              oe_mis_syn=syn_ht.variant_count / syn_ht.expected_variants)

    ht = lof_ht_classic.annotate(**mis_ht[lof_ht_classic.key], **pphen_mis_ht[lof_ht_classic.key],
                                 **syn_ht[lof_ht_classic.key], **lof_ht[lof_ht_classic.key],
                                 **lof_ht_classic_hc[lof_ht_classic.key])
    ht = ht.annotate(**oe_confidence_interval(ht, ht.obs_lof_classic_hc, ht.exp_lof_classic_hc,
                                              prefix='oe_classic_hc')[ht.key])
    return ht.annotate(**oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof)[ht.key])


def collapse_lof_ht(lof_ht: hl.Table, keys: Tuple[str]) -> hl.Table:
    lof_ht = lof_ht.group_by(*keys).aggregate(obs_lof=hl.agg.sum(lof_ht.variant_count),
                                              exp_lof=hl.agg.sum(lof_ht.expected_variants),
                                              possible_lof=hl.agg.sum(lof_ht.possible_variants),
                                              adjusted_mu_lof=hl.agg.sum(
                                                  lof_ht.adjusted_mutation_rate)).persist()
    lof_ht = lof_ht.filter(lof_ht.exp_lof > 0)
    return lof_ht.annotate(
        **pLI(lof_ht, lof_ht.obs_lof, lof_ht.exp_lof)[lof_ht.key],
        oe_lof=lof_ht.obs_lof / lof_ht.exp_lof).key_by(*keys)


def annotate_constraint_groupings(ht: Union[hl.Table, hl.MatrixTable],
                                  custom_model: str = None) -> Tuple[Union[hl.Table, hl.MatrixTable], List[str]]:
    """
    HT must be exploded against whatever axis

    Need to add `'coverage': ht.exome_coverage` here (which will get corrected out later)
    """
    if custom_model == 'worst_csq':
        groupings = {
            'annotation': ht.vep.worst_csq_by_gene.most_severe_consequence,
            'modifier': hl.case()
                .when(hl.is_defined(ht.vep.worst_csq_by_gene.polyphen_prediction),
                      ht.vep.worst_csq_by_gene.polyphen_prediction)
                .when(hl.is_defined(ht.vep.worst_csq_by_gene.lof),
                      ht.vep.worst_csq_by_gene.lof)
                .default('None'),
            'gene': ht.vep.worst_csq_by_gene.gene_symbol,
            'coverage': ht.exome_coverage
        }
    else:
        groupings = {
            'annotation': ht.vep.transcript_consequences.most_severe_consequence,
            'modifier': hl.case()
                .when(hl.is_defined(ht.vep.transcript_consequences.polyphen_prediction),
                      ht.vep.transcript_consequences.polyphen_prediction)
                .when(hl.is_defined(ht.vep.transcript_consequences.lof),
                      ht.vep.transcript_consequences.lof)
                .default('None'),
            'transcript': ht.vep.transcript_consequences.transcript_id,
            'gene': ht.vep.transcript_consequences.gene_symbol,
            'canonical': hl.or_else(ht.vep.transcript_consequences.canonical == 1, False),
            'coverage': ht.exome_coverage
        }
        if custom_model == 'splice_region':
            groupings['distance_splice'] = ht.vep.transcript_consequences
    ht = ht.annotate(**groupings) if isinstance(ht, hl.Table) else ht.annotate_rows(**groupings)
    return ht, list(groupings.keys())


# Model building
def build_coverage_model(coverage_ht: hl.Table) -> (float, float):
    """
    Calibrates coverage model (returns intercept and slope)
    """
    return tuple(coverage_ht.aggregate(hl.agg.linreg(coverage_ht.low_coverage_obs_exp, [1, coverage_ht.log_coverage])).beta)


def build_plateau_models(ht: hl.Table, weighted: bool = False) -> Dict[str, Tuple[float, float]]:
    """
    Calibrates high coverage model (returns intercept and slope)
    """
    # TODO: try square weighting
    return ht.aggregate(hl.agg.group_by(ht.cpg,
                                        hl.agg.linreg(ht.high_coverage_proportion_observed, [1, ht.mu_snp],
                                                      weight=ht.possible_variants if weighted else None)
                                       ).map_values(lambda x: x.beta))


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


def pLI(ht: hl.Table, obs: hl.expr.Int32Expression, exp: hl.expr.Float32Expression) -> hl.Table:
    last_pi = {'Null': 0, 'Rec': 0, 'LI': 0}
    pi = {'Null': 1 / 3, 'Rec': 1 / 3, 'LI': 1 / 3}
    expected_values = {'Null': 1, 'Rec': 0.463, 'LI': 0.089}
    ht = ht.annotate(_obs=obs, _exp=exp)

    while abs(pi['LI'] - last_pi['LI']) > 0.001:
        last_pi = copy.deepcopy(pi)
        ht = ht.annotate(
            **{k: v * hl.dpois(ht._obs, ht._exp * expected_values[k]) for k, v in pi.items()})
        ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
        ht = ht.annotate(**{k: ht[k] / ht.row_sum for k, v in pi.items()})
        pi = ht.aggregate({k: hl.agg.mean(ht[k]) for k in pi.keys()})

    ht = ht.annotate(
        **{k: v * hl.dpois(ht._obs, ht._exp * expected_values[k]) for k, v in pi.items()})
    ht = ht.annotate(row_sum=hl.sum([ht[k] for k in pi]))
    return ht.select(**{f'p{k}': ht[k] / ht.row_sum for k, v in pi.items()})


def oe_confidence_interval(ht: hl.Table, obs: hl.expr.Int32Expression, exp: hl.expr.Float32Expression,
                           prefix: str = 'oe', alpha: float = 0.05) -> hl.Table:
    ht = ht.annotate(_obs=obs, _exp=exp)
    oe_ht = ht.annotate(range=hl.range(0, 2000).map(lambda x: hl.float64(x) / 1000))
    oe_ht = oe_ht.annotate(range_dpois=oe_ht.range.map(lambda x: hl.dpois(oe_ht._obs, oe_ht._exp * x)))

    oe_ht = oe_ht.annotate(cumulative_dpois=hl.cumulative_sum(oe_ht.range_dpois))
    max_cumulative_dpois = oe_ht.cumulative_dpois[-1]
    oe_ht = oe_ht.annotate(norm_dpois=oe_ht.cumulative_dpois.map(lambda x: x / max_cumulative_dpois))
    oe_ht = oe_ht.annotate(
        lower_idx=hl.argmax(oe_ht.norm_dpois.map(lambda x: hl.or_missing(x < alpha, x))),
        upper_idx=hl.argmin(oe_ht.norm_dpois.map(lambda x: hl.or_missing(x > 1 - alpha, x)))
    )
    return oe_ht.select(**{
        f'{prefix}_lower': oe_ht.range[oe_ht.lower_idx],
        f'{prefix}_upper': oe_ht.range[oe_ht.upper_idx]
    })


def calculate_z(input_ht: hl.Table, observed: hl.expr.NumericExpression, expected: hl.expr.NumericExpression, output: str = 'z') -> hl.Table:
    input_ht = input_ht.annotate(_obs=observed, _exp=expected)
    ht = input_ht.filter(input_ht._obs > 0)
    ht = ht.select(_z=ht._obs / ht._exp)
    return input_ht.annotate(
        reasons=hl.case().when(input_ht._obs == 0, 'zero_variants').or_missing(),
        **{output: ht[input_ht.key]._z}
    ).drop('_obs', '_exp')
