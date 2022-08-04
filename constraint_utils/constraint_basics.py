__author__ = 'konradk'

from .generic import *
import pickle
import copy
import uuid

root = 'gs://gnomad-public/papers/2019-flagship-lof/v1.0'

# Unprocessed files
fasta_path = f'{root}/context/source/Homo_sapiens_assembly19.fasta'
gerp_annotations_path = f'{root}/annotations/gerp.scores.GRCh37.ht'  # Gerp is in here as S

# Input datasets
raw_context_txt_path = f'{root}/context/source/Homo_sapiens_assembly19.fasta.snps_only.vep.txt.bgz/*'
raw_context_ht_path = f'{root}/context/Homo_sapiens_assembly19.fasta.snps_only.unsplit.ht'
vep_context_ht_path = f'{root}/context/Homo_sapiens_assembly19.fasta.snps_only.unsplit.vep_20181129.ht'
context_ht_path = f'{root}/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht'
processed_genomes_ht_path = f'{root}/model/genomes_processed.ht'
processed_exomes_ht_path = f'{root}/model/exomes_processed.ht'

# Possible variant files
all_possible_summary_pickle = f'{root}/model/tmp/all_possible_counts_by_context.pckl'
all_possible_summary_unfiltered_pickle = f'{root}/model/tmp/all_possible_counts_by_context_unfiltered.pckl'

mutation_rate_ht_path = f'{root}/model/mutation_rate_methylation_bins.ht'
po_coverage_ht_path = f'{root}/model/prop_observed_by_coverage_no_common_pass_filtered_bins.ht'
po_ht_path = f'{root}/{{subdir}}/prop_observed_{{subdir}}.ht'
raw_constraint_ht_path = f'{root}/{{subdir}}/constraint_{{subdir}}.ht'
final_constraint_ht_path = f'{root}/{{subdir}}/constraint_final_{{subdir}}.ht'

HIGH_COVERAGE_CUTOFF = 40
VARIANT_TYPES_FOR_MODEL = ('ACG', 'TCG', 'CCG', 'GCG', 'non-CpG')
POPS = ('global', 'afr', 'amr', 'eas', 'nfe', 'sas')
MODEL_KEYS = {
    'worst_csq': ['gene'],
    'tx_annotation': ['gene', 'expressed'],
    'standard': ['gene', 'transcript', 'canonical']
}


# Data loads
def get_old_mu_data() -> hl.Table:
    old_mu_data = hl.import_table(f'{root}/old_exac_data/fordist_1KG_mutation_rate_table.txt',
                                  delimiter=' ', impute=True)
    return old_mu_data.transmute(context=old_mu_data['from'], ref=old_mu_data['from'][1],
                                 alt=old_mu_data.to[1]).key_by('context', 'ref', 'alt')


def load_all_possible_summary(filtered: bool = True) -> Dict[hl.Struct, int]:
    fname = all_possible_summary_pickle if filtered else all_possible_summary_unfiltered_pickle
    with hl.hadoop_open(fname, 'rb') as f:
        return pickle.load(f)


def load_tx_expression_data(context=True):
    tx_ht_file = 'gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.021819.ht' if context \
        else 'gs://gnomad-public/papers/2019-tx-annotation/gnomad.exomes.r2.1.1.sites.tx_annotated.021319.ht'
    tx_ht = hl.read_matrix_table(tx_ht_file).rows()

    def process_expression_data(csq_expression):
        exprs_to_drop = ['ensg', 'csq', 'symbol', 'lof', 'lof_flag', 'mean_proportion']
        expression_data = csq_expression.drop(*exprs_to_drop)
        all_tissues = list(expression_data.values())
        expression_data_list = list(zip(list(expression_data), all_tissues))
        brain_tissues = [x[1] for x in expression_data_list if 'Brain' in x[0]]
        return csq_expression.select('ensg', 'csq', 'symbol', 'lof',
                                     mean_expression=hl.mean(hl.filter(lambda e: ~hl.is_nan(e), all_tissues), filter_missing=True),
                                     max_expression=hl.max(hl.filter(lambda e: ~hl.is_nan(e), all_tissues), filter_missing=True),
                                     mean_brain_expression=hl.mean(hl.filter(lambda k: ~hl.is_nan(k), brain_tissues), filter_missing=True),
                                     Brain_Cortex=csq_expression.Brain_Cortex
                                     )

    return tx_ht.annotate(tx_annotation=tx_ht.tx_annotation.map(process_expression_data))


def annotate_distance_to_splice(input_ht):
    # Load GTF file
    tmp_path = f'/tmp_{uuid.uuid4()}.ht'
    ht = hl.experimental.import_gtf('gs://hail-common/references/gencode/gencode.v19.annotation.gtf.bgz', 'GRCh37', True, min_partitions=100)
    ht = ht.filter((ht.feature == 'exon') & (ht.transcript_type == 'protein_coding')).key_by()
    ht = ht.select(gene_id=ht.gene_id.split('\\.')[0], transcript_id=ht.transcript_id.split('\\.')[0],
                   loc=[hl.struct(acceptor=ht.strand == '+', locus=ht.interval.start),
                        hl.struct(acceptor=ht.strand != '+', locus=ht.interval.end)]).explode('loc')
    ht.transmute(**ht.loc).key_by('locus').collect_by_key().write(tmp_path, True)

    # Scan to get bounds, create intervals
    tmp_path2 = f'/tmp_{uuid.uuid4()}.ht'
    ht = hl.read_table(tmp_path)
    last_locus = hl.scan.take(ht.row, 1, ordering=-ht.locus.global_position())
    ht.key_by(
        interval=hl.or_missing((hl.len(last_locus) > 0) &
                               (last_locus[0].contig == ht.locus.contig),
                               hl.interval(last_locus[0], ht.locus))
    ).write(tmp_path2)

    ht = hl.read_table(tmp_path2)
    # join against intervals
    data = ht[input_ht.locus]
    # find nearer splice site
    return input_ht.transmute(
        nearest_splice=hl.cond(
            hl.abs(input_ht.locus.position - input_ht.data.locus.position) <= hl.abs(
                input_ht.locus.position - data.last.locus.position),
            hl.struct(dist=input_ht.locus.position - data.locus.position, sites=input_ht.data.values),
            hl.struct(dist=input_ht.locus.position - data.last.locus.position, sites=input_ht.data.last.values)
        )
    )


def annotate_tx_expression_data(ht, tx_ht, location):
    key = ht.key if isinstance(ht, hl.Table) else ht.row_key
    return hl.find(lambda csq: (csq.ensg == location.gene_id) &
                               (csq.csq == location.most_severe_consequence),
                   tx_ht[key].tx_annotation)


# Pre-process
def export_fasta(hc) -> None:
    # Done with an 0.1 jar
    hc.read('gs://gnomad-resources/Homo_sapiens_assembly19.fasta.snps_only.vep.vds').export_variants('gs://gnomad-resources/context/source/Homo_sapiens_assembly19.fasta.snps_only.vep.txt.bgz', 'v, va.context, va.vep', types=True, parallel=True)


def import_fasta(raw_context_txt_path: str, raw_context_ht_path: str, overwrite: bool = False) -> None:
    # Works with --jar gs://hail-common/builds/devel/jars/hail-devel-47de006f2c62-Spark-2.2.0.jar
    ht = hl.import_table(raw_context_txt_path, no_header=True, min_partitions=40000,
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
    ht.transmute(**hl.parse_variant(ht.v)).key_by('locus', 'alleles').write(raw_context_ht_path, overwrite)


def vep_context_ht(raw_context_ht_path: str, vep_context_ht_path: str, overwrite: bool = False):
    ht_full = hl.read_table(raw_context_ht_path)
    chunks = ('1-2', '3', '4', '5-6', '7-9', '10-12', '13-15', '16-18', '19-22', 'X-Y')
    for i in chunks:
        print(i)
        ht = hl.filter_intervals(ht_full, [hl.parse_locus_interval(i)])
        hl.vep(ht, vep_config).write(f'{root}/context/parts/Homo_sapiens_assembly19.fasta.snps_only.{i}.ht', overwrite=overwrite, stage_locally=True)
    hts = [hl.read_table(f'{root}/context/parts/Homo_sapiens_assembly19.fasta.snps_only.{i}.ht') for i in chunks]
    ht = hts[0].union(*hts[1:])
    ht.write(vep_context_ht_path, overwrite)


def split_context_mt(raw_context_ht_path: str, coverage_ht_paths: Dict[str, str], methylation_ht_path: str,
                     context_ht_path: str, overwrite: bool = False) -> None:
    raw_context_ht = hl.split_multi_hts(hl.read_table(raw_context_ht_path))
    raw_context_ht.write(f'{raw_context_ht_path}.temp_split.ht', overwrite)
    raw_context_ht = hl.read_table(f'{raw_context_ht_path}.temp_split.ht')

    coverage_hts = {loc: hl.read_table(coverage_ht_path) for loc, coverage_ht_path in coverage_ht_paths.items()}
    coverage_hts = {loc: coverage_ht.drop('#chrom', 'pos') if '#chrom' in list(coverage_ht.row) else coverage_ht
                    for loc, coverage_ht in coverage_hts.items()}
    methylation_ht = hl.read_table(methylation_ht_path)
    gerp_ht = hl.read_table(gerp_annotations_path)

    raw_context_ht = raw_context_ht.annotate(
        methylation=methylation_ht[raw_context_ht.locus],
        coverage=hl.struct(**{loc: coverage_ht[raw_context_ht.locus] for loc, coverage_ht in coverage_hts.items()}),
        gerp=gerp_ht[raw_context_ht.locus].S)

    raw_context_ht.write(context_ht_path, overwrite)


def pre_process_data(ht: hl.Table, split_context_ht_path: str,
                     output_ht_path: str, overwrite: bool = False) -> None:
    """ 
    Add annotations from VEP context Table to gnomAD data.
    
    Function adds the following annotations:
        - context
        - methylation
        - coverage
        - gerp
        - pass_filters
    
    Function drops `a_index`, `was_split`, and`colocated_variants` annotations from gnomAD data.
    
    .. note::
        Function expects that multiallelic variants in the VEP context Table have been split.
    
    :param ht: gnomAD exomes or genomes public Hail Table. 
    :param split_context_ht_path: Path to VEP context Table.
    :param output_ht_path: Path to output Table.
    :param overwrite: Whether to overwrite existing data. Defaults to False.
    """
    context_ht = hl.read_table(split_context_ht_path).drop('a_index', 'was_split')
    context_ht = context_ht.annotate(vep=context_ht.vep.drop('colocated_variants'))
    ht.annotate(**context_ht[ht.key], pass_filters=hl.len(ht.filters) == 0).write(output_ht_path, overwrite)


def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    """
    Filter input Table and add annotations used in constraint calculations.
 
    Function filters to SNPs, removes rows with undefined contexts, collapses strands
    to deduplicate trimer or heptamer contexts, and annotates the input Table.

    Function adds the following annotations:
        - ref
        - alt
        - was_flipped
        - methylation_level
        - cpg
        - transition
        - variant_type
        - variant_type_model
        - exome_coverage
   
    :param ht: Input Table to be annotated.
    :param trimer: Whether to use trimers or heptamers. Defaults to False.
    :param annotate_coverage: Whether to annotate the coverage of exome. Defaults to True.
    :return: Table with annotations.
    """
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
    }
    if annotate_coverage:
        annotation['exome_coverage'] = ht.coverage.exomes.median
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
        context_ht = context_ht.filter(hl.all(lambda x: x, [hl.is_defined(context_ht[x]) for x in grouping_variables]))
        genome_ht = genome_ht.filter(hl.all(lambda x: x, [hl.is_defined(genome_ht[x]) for x in grouping_variables]))
    else:
        context_ht = filter_for_mu(context_ht)
        genome_ht = filter_for_mu(genome_ht)

    context_ht = context_ht.select('context', 'ref', 'alt', 'methylation_level', *grouping_variables)
    genome_ht = genome_ht.select('context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping_variables)

    genome_join = genome_ht[context_ht.key]
    if remove_common_downsampled:
        ac_cutoff = 5
        downsampling_level = 1000
        freq_index = hl.eval(
            hl.zip_with_index(genome_ht.freq_meta).find(lambda f:
                                                        (f[1].get('downsampling') == str(downsampling_level)) &
                                                        (f[1].get('pop') == 'global') &
                                                        (f[1].get('group') == 'adj') &
                                                        (f[1].size() == 3)))[0]
        # In the denominator, only keep variants not in the genome dataset, or with AC <= ac_cutoff and passing filters
        context_ht = context_ht.filter(
            hl.is_missing(genome_join) | ((genome_join.freq[freq_index].AC <= ac_cutoff) & genome_join.pass_filters)
        )
        # Keep only AC <= ac_cutoff in numerator
        genome_ht = genome_ht.filter(genome_ht.freq[freq_index].AC <= ac_cutoff)
    elif remove_common_ordinary:
        af_cutoff = 0.001
        context_ht = context_ht.filter(
            hl.is_missing(genome_join) | (
                (genome_join.freq[0].AF <= af_cutoff) & genome_join.pass_filters
            )
        )
        # Keep only AF < af_cutoff in numerator
        genome_ht = filter_by_frequency(genome_ht, 'below', frequency=af_cutoff, adj=True)
    else:
        context_ht = context_ht.filter(
            hl.is_missing(genome_join) | genome_join.pass_filters
        )
    genome_ht = genome_ht.filter(genome_ht.pass_filters)

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
    downsamplings = list(map(lambda x: x[1], get_downsamplings(ht)))
    index_1kg = downsamplings.index(1000)
    return ht.annotate(proportion_observed_1kg=ht.downsampling_counts_global[index_1kg] / ht.possible_variants,
                       proportion_observed=ht.variant_count / ht.possible_variants,
                       mu_snp=ht.downsamplings_mu_snp[index_1kg],
                       mu_snp_nfe=ht.downsamplings_mu_nfe[index_1kg],
                       mu_snp_afr=ht.downsamplings_mu_afr[index_1kg])


def get_proportion_observed_by_coverage(exome_ht: hl.Table, context_ht: hl.Table, mutation_ht: hl.Table,
                                        recompute_possible: bool = False, dataset: str = 'gnomad',
                                        impose_high_af_cutoff_upfront: bool = True) -> hl.Table:

    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    context_ht = add_most_severe_csq_to_tc_within_ht(context_ht)

    context_ht = fast_filter_vep(context_ht).select('context', 'ref', 'alt', 'methylation_level', 'exome_coverage')
    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage))

    exome_ht = fast_filter_vep(exome_ht).select('context', 'ref', 'alt', 'methylation_level', 'exome_coverage', 'freq', 'pass_filters')

    grouping = ('exome_coverage',)
    af_cutoff = 0.001

    exome_join = exome_ht[context_ht.key]
    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]

    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit
    context_ht = context_ht.filter(hl.is_missing(exome_join) | keep_criteria(exome_join))

    exome_ht = exome_ht.filter(keep_criteria(exome_ht))

    possible_file = f'{root}/tmp/possible_coverage.ht'
    if recompute_possible:
        possible_ht = count_variants(context_ht, additional_grouping=grouping, force_grouping=True)
        possible_ht = annotate_with_mu(possible_ht, mutation_ht)
        possible_ht.transmute(possible_variants=possible_ht.variant_count).write(possible_file, True)

    possible_ht = hl.read_table(possible_file)
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=100, count_downsamplings=POPS,
                        impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront)
    ht = ht.join(possible_ht, 'outer')
    return ht


def build_models(coverage_ht: hl.Table, trimers: bool = False, weighted: bool = False, half_cutoff = False,
                 ) -> Tuple[Tuple[float, float], Dict[str, Tuple[float, float]]]:
    keys = ['context', 'ref', 'alt', 'methylation_level', 'mu_snp']

    cov_cutoff = (HIGH_COVERAGE_CUTOFF / half_cutoff) if half_cutoff else HIGH_COVERAGE_CUTOFF
    all_high_coverage_ht = coverage_ht.filter(coverage_ht.exome_coverage >= cov_cutoff)
    agg_expr = {
        'observed_variants': hl.agg.sum(all_high_coverage_ht.variant_count),
        'possible_variants': hl.agg.sum(all_high_coverage_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'observed_{pop}'] = hl.agg.array_sum(all_high_coverage_ht[f'downsampling_counts_{pop}'])
    high_coverage_ht = all_high_coverage_ht.group_by(*keys).aggregate(**agg_expr)

    high_coverage_ht = annotate_variant_types(high_coverage_ht, not trimers)
    plateau_models = build_plateau_models_pop(high_coverage_ht, weighted=weighted)

    high_coverage_scale_factor = all_high_coverage_ht.aggregate(
        hl.agg.sum(all_high_coverage_ht.variant_count) /
        hl.agg.sum(all_high_coverage_ht.possible_variants * all_high_coverage_ht.mu_snp))

    all_low_coverage_ht = coverage_ht.filter((coverage_ht.exome_coverage < cov_cutoff) &
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


def get_proportion_observed(exome_ht: hl.Table, context_ht: hl.Table, mutation_ht: hl.Table,
                            plateau_models: Dict[str, Tuple[float, float]], coverage_model: Tuple[float, float],
                            recompute_possible: bool = False, remove_from_denominator: bool = True,
                            custom_model: str = None, dataset: str = 'gnomad',
                            impose_high_af_cutoff_upfront: bool = True, half_cutoff = False) -> hl.Table:

    exome_ht = add_most_severe_csq_to_tc_within_ht(exome_ht)
    context_ht = add_most_severe_csq_to_tc_within_ht(context_ht)

    if custom_model == 'syn_canonical':
        context_ht = take_one_annotation_from_tc_within_ht(fast_filter_vep(context_ht))
        context_ht = context_ht.transmute(transcript_consequences=context_ht.vep.transcript_consequences)
        exome_ht = take_one_annotation_from_tc_within_ht(fast_filter_vep(exome_ht))
        exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)
    elif custom_model == 'worst_csq':
        context_ht = process_consequences(context_ht)
        context_ht = context_ht.transmute(worst_csq_by_gene=context_ht.vep.worst_csq_by_gene)
        context_ht = context_ht.explode(context_ht.worst_csq_by_gene)
        exome_ht = process_consequences(exome_ht)
        exome_ht = exome_ht.transmute(worst_csq_by_gene=exome_ht.vep.worst_csq_by_gene)
        exome_ht = exome_ht.explode(exome_ht.worst_csq_by_gene)
    elif custom_model == 'tx_annotation':
        tx_ht = load_tx_expression_data()
        context_ht = context_ht.annotate(**tx_ht[context_ht.key])
        context_ht = context_ht.explode(context_ht.tx_annotation)
        tx_ht = load_tx_expression_data(context=False)
        exome_ht = exome_ht.annotate(**tx_ht[exome_ht.key])
        exome_ht = exome_ht.explode(exome_ht.tx_annotation)
    elif custom_model == 'distance_to_splice':
        context_ht = annotate_distance_to_splice(context_ht)
        exome_ht = annotate_distance_to_splice(exome_ht)
        # TODO:
        # context_ht = context_ht.explode(context_ht.tx_annotation)
        # exome_ht = exome_ht.explode(exome_ht.tx_annotation)
    else:
        context_ht = context_ht.transmute(transcript_consequences=context_ht.vep.transcript_consequences)
        context_ht = context_ht.explode(context_ht.transcript_consequences)
        exome_ht = exome_ht.transmute(transcript_consequences=exome_ht.vep.transcript_consequences)
        exome_ht = exome_ht.explode(exome_ht.transcript_consequences)

    context_ht, _ = annotate_constraint_groupings(context_ht, custom_model=custom_model)
    exome_ht, grouping = annotate_constraint_groupings(exome_ht, custom_model=custom_model)

    context_ht = context_ht.filter(hl.is_defined(context_ht.exome_coverage)).select(
        'context', 'ref', 'alt', 'methylation_level', *grouping)
    exome_ht = exome_ht.select(
        'context', 'ref', 'alt', 'methylation_level', 'freq', 'pass_filters', *grouping)

    af_cutoff = 0.001

    freq_index = exome_ht.freq_index_dict.collect()[0][dataset]

    def keep_criteria(ht):
        crit = (ht.freq[freq_index].AC > 0) & ht.pass_filters & (ht.coverage > 0)
        if impose_high_af_cutoff_upfront:
            crit &= (ht.freq[freq_index].AF <= af_cutoff)
        return crit

    exome_join = exome_ht[context_ht.key]
    if remove_from_denominator:
        context_ht = context_ht.filter(hl.is_missing(exome_join) | keep_criteria(exome_join))

    exome_ht = exome_ht.filter(keep_criteria(exome_ht))

    possible_file = f'{root}/model/possible_data/possible_transcript_pop_{custom_model}.ht'
    if recompute_possible:
        ht = count_variants(context_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True)
        ht = annotate_with_mu(ht, mutation_ht)
        ht = ht.transmute(possible_variants=ht.variant_count)
        ht = annotate_variant_types(ht.annotate(mu_agg=ht.mu_snp * ht.possible_variants))
        model = hl.literal(plateau_models.total)[ht.cpg]
        cov_cutoff = (HIGH_COVERAGE_CUTOFF / half_cutoff) if half_cutoff else HIGH_COVERAGE_CUTOFF
        ann_expr = {
            'adjusted_mutation_rate': ht.mu_agg * model[1] + model[0],
            'coverage_correction': hl.case()
                .when(ht.coverage == 0, 0)
                .when(ht.coverage >= cov_cutoff, 1)
                .default(coverage_model[1] * hl.log10(ht.coverage) + coverage_model[0])
        }
        for pop in POPS:
            pop_model = hl.literal(plateau_models[pop])
            slopes = hl.map(lambda f: f[ht.cpg][1], pop_model)
            intercepts = hl.map(lambda f: f[ht.cpg][0], pop_model)
            ann_expr[f'adjusted_mutation_rate_{pop}'] = ht.mu_agg * slopes + intercepts
        ht = ht.annotate(**ann_expr)
        ann_expr = {
            'expected_variants': ht.adjusted_mutation_rate * ht.coverage_correction,
            'mu': ht.mu_agg * ht.coverage_correction
        }
        for pop in POPS:
            ann_expr[f'expected_variants_{pop}'] = ht[f'adjusted_mutation_rate_{pop}'] * ht.coverage_correction
        ht = ht.annotate(**ann_expr)
        ht.write(possible_file, True)

    possible_variants_ht = hl.read_table(possible_file)
    ht = count_variants(exome_ht, additional_grouping=grouping, partition_hint=2000, force_grouping=True,
                        count_downsamplings=POPS, impose_high_af_cutoff_here=not impose_high_af_cutoff_upfront)
    ht = ht.join(possible_variants_ht, 'outer')
    ht.write(f'{root}/model/possible_data/all_data_transcript_pop_{custom_model}.ht', True)
    ht = hl.read_table(f'{root}/model/possible_data/all_data_transcript_pop_{custom_model}.ht')

    grouping.remove('coverage')
    agg_expr = {
        'variant_count': hl.agg.sum(ht.variant_count),
        'adjusted_mutation_rate': hl.agg.sum(ht.adjusted_mutation_rate),
        'possible_variants': hl.agg.sum(ht.possible_variants),
        'expected_variants': hl.agg.sum(ht.expected_variants),
        'mu': hl.agg.sum(ht.mu)
    }
    for pop in POPS:
        agg_expr[f'adjusted_mutation_rate_{pop}'] = hl.agg.array_sum(ht[f'adjusted_mutation_rate_{pop}'])
        agg_expr[f'expected_variants_{pop}'] = hl.agg.array_sum(ht[f'expected_variants_{pop}'])
        agg_expr[f'downsampling_counts_{pop}'] = hl.agg.array_sum(ht[f'downsampling_counts_{pop}'])
    ht = ht.group_by(*grouping).partition_hint(1000).aggregate(**agg_expr)
    return ht.annotate(obs_exp=ht.variant_count / ht.expected_variants)


def finalize_dataset(po_ht: hl.Table, keys: Tuple[str] = ('gene', 'transcript', 'canonical'),
                     n_partitions: int = 1000) -> hl.Table:
    """
    Compute the pLI scores, 90% confidence interval around the observed:expected ratio, and z scores 
    for synomynous variants, missense variants, and LoF variants.
    
    .. note::
        The following annotations should be present in `po_ht`:
            - modifier
            - annotation
            - variant_count
            - mu
            - possible_variants
            - expected_variants
            - expected_variants_{pop} (pop defaults to `POPS`)
            - downsampling_counts_{pop} (pop defaults to `POPS`)
 
    :param po_ht: Input table with the number of expected variants (output of `get_proportion_observed`).
    :param keys: The keys of the output table, defaults to ('gene', 'transcript', 'canonical').
    :param n_partitions: Desired number of partitions for `Table.repartition()`, defaults to 1000.
    :return: Table with pLI scores, confidence inteval of oe ratio, and z scores.
    """
    # This function aggregates over genes in all cases, as XG spans PAR and non-PAR X
    po_ht = po_ht.repartition(n_partitions).persist()

    # Getting classic LoF annotations (no LOFTEE)
    classic_lof_annotations = hl.literal({'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant'})
    lof_ht_classic = po_ht.filter(classic_lof_annotations.contains(po_ht.annotation) &
                                  ((po_ht.modifier == 'HC') | (po_ht.modifier == 'LC')))
    lof_ht_classic = collapse_lof_ht(lof_ht_classic, keys, False)
    lof_ht_classic = lof_ht_classic.rename({x: f'{x}_classic' for x in list(lof_ht_classic.row_value)})

    # Getting all LoF annotations (LOFTEE HC + OS)
    lof_ht_classic_hc = po_ht.filter((po_ht.modifier == 'HC') | (po_ht.modifier == 'OS'))
    lof_ht_classic_hc = collapse_lof_ht(lof_ht_classic_hc, keys, False)
    lof_ht_classic_hc = lof_ht_classic_hc.rename({x: f'{x}_with_os' for x in list(lof_ht_classic_hc.row_value)})

    # Getting all LoF annotations (LOFTEE HC)
    lof_ht = po_ht.filter(po_ht.modifier == 'HC')
    lof_ht = collapse_lof_ht(lof_ht, keys, False)

    mis_ht = po_ht.filter(po_ht.annotation == 'missense_variant')
    agg_expr = {
        'obs_mis': hl.agg.sum(mis_ht.variant_count),
        'exp_mis': hl.agg.sum(mis_ht.expected_variants),
        'oe_mis': hl.agg.sum(mis_ht.variant_count) / hl.agg.sum(mis_ht.expected_variants),
        'mu_mis': hl.agg.sum(mis_ht.mu),
        'possible_mis': hl.agg.sum(mis_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_mis_{pop}'] = hl.agg.array_sum(mis_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_mis_{pop}'] = hl.agg.array_sum(mis_ht[f'downsampling_counts_{pop}'])
    mis_ht = mis_ht.group_by(*keys).aggregate(**agg_expr)

    pphen_mis_ht = po_ht.filter(po_ht.modifier == 'probably_damaging')
    pphen_mis_ht = pphen_mis_ht.group_by(*keys).aggregate(obs_mis_pphen=hl.agg.sum(pphen_mis_ht.variant_count),
                                                          exp_mis_pphen=hl.agg.sum(pphen_mis_ht.expected_variants),
                                                          oe_mis_pphen=hl.agg.sum(pphen_mis_ht.variant_count) / hl.agg.sum(pphen_mis_ht.expected_variants),
                                                          possible_mis_pphen=hl.agg.sum(pphen_mis_ht.possible_variants))
    syn_ht = po_ht.filter(po_ht.annotation == 'synonymous_variant').key_by(*keys)
    agg_expr = {
        'obs_syn': hl.agg.sum(syn_ht.variant_count),
        'exp_syn': hl.agg.sum(syn_ht.expected_variants),
        'oe_syn': hl.agg.sum(syn_ht.variant_count) / hl.agg.sum(syn_ht.expected_variants),
        'mu_syn': hl.agg.sum(syn_ht.mu),
        'possible_syn': hl.agg.sum(syn_ht.possible_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_syn_{pop}'] = hl.agg.array_sum(syn_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_syn_{pop}'] = hl.agg.array_sum(syn_ht[f'downsampling_counts_{pop}'])
    syn_ht = syn_ht.group_by(*keys).aggregate(**agg_expr)

    ht = lof_ht_classic.annotate(**mis_ht[lof_ht_classic.key], **pphen_mis_ht[lof_ht_classic.key],
                                 **syn_ht[lof_ht_classic.key], **lof_ht[lof_ht_classic.key],
                                 **lof_ht_classic_hc[lof_ht_classic.key])
    syn_cis = oe_confidence_interval(ht, ht.obs_syn, ht.exp_syn, prefix='oe_syn')
    mis_cis = oe_confidence_interval(ht, ht.obs_mis, ht.exp_mis, prefix='oe_mis')
    lof_cis = oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof, prefix='oe_lof')
    ht = ht.annotate(**syn_cis[ht.key], **mis_cis[ht.key], **lof_cis[ht.key])
    return calculate_all_z_scores(ht)  # .annotate(**oe_confidence_interval(ht, ht.obs_lof, ht.exp_lof)[ht.key])


def collapse_lof_ht(lof_ht: hl.Table, keys: Tuple[str], calculate_pop_pLI: bool = False) -> hl.Table:
    """
    Computes the pLI scores using observed and expected varinat counts.
    
    Function sums the number of observed variants, possible variants, and expected variants across all the `keys`,
    and use the expected variant counts and obverved variant counts to compute the pLI scores.
            
    Function will add the following annotations to output Table:
        - obs_lof - the sum of observed variants
        - mu_lof - the sum of mutation rate
        - possible_lof - possible number of LoF variants
        - exp_lof - expected number of LoF variants
        - exp_lof_{pop} (pop defaults to `POPS`) - expected number of LoF variants per population
        - obs_lof_{pop} (pop defaults to `POPS`) - observed number of LoF variants per population
        - oe_lof - lof observed:expected ratio
        - annotations added by function `pLI()`

    .. note::
        The following annotations should be present in `lof_ht`:
            - variant_count
            - mu
            - possible_variants
            - expected_variants
            
    :param lof_ht: Table with specific LoF annotations.
    :param keys: The keys of ouput Table.
    :param calculate_pop_pLI: Whether to calculate the pLI score for each population, defaults to False.
    :return: A collapsed Table with pLI scores.
    """
    agg_expr = {
        'obs_lof': hl.agg.sum(lof_ht.variant_count),
        'mu_lof': hl.agg.sum(lof_ht.mu),
        'possible_lof': hl.agg.sum(lof_ht.possible_variants),
        'exp_lof': hl.agg.sum(lof_ht.expected_variants)
    }
    for pop in POPS:
        agg_expr[f'exp_lof_{pop}'] = hl.agg.array_sum(lof_ht[f'expected_variants_{pop}'])
        agg_expr[f'obs_lof_{pop}'] = hl.agg.array_sum(lof_ht[f'downsampling_counts_{pop}'])
    lof_ht = lof_ht.group_by(*keys).aggregate(**agg_expr).persist()
    lof_ht = lof_ht.filter(lof_ht.exp_lof > 0)
    if calculate_pop_pLI:
        pop_lengths = get_all_pop_lengths(lof_ht, 'obs_lof_')
        print(pop_lengths)
        for pop_length, pop in pop_lengths:
            print(f'Calculating pLI for {pop}...')
            plis = []
            for i in range(8, pop_length):
                print(i)
                ht = lof_ht.filter(lof_ht[f'exp_lof_{pop}'][i] > 0)
                pli_ht = pLI(ht, ht[f'obs_lof_{pop}'][i], ht[f'exp_lof_{pop}'][i])
                plis.append(pli_ht[lof_ht.key])
            lof_ht = lof_ht.annotate(**{
                f'pLI_{pop}': [pli.pLI for pli in plis],
                f'pRec_{pop}': [pli.pRec for pli in plis],
                f'pNull_{pop}': [pli.pNull for pli in plis],
            })
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
            'annotation': ht.worst_csq_by_gene.most_severe_consequence,
            'modifier': hl.case()
                .when(hl.is_defined(ht.worst_csq_by_gene.lof),
                      ht.worst_csq_by_gene.lof)
                .when(hl.is_defined(ht.worst_csq_by_gene.polyphen_prediction),
                      ht.worst_csq_by_gene.polyphen_prediction)
                .default('None'),
            'gene': ht.worst_csq_by_gene.gene_symbol,
            'coverage': ht.exome_coverage
        }
    elif custom_model == 'tx_annotation':
        groupings = {
            'annotation': ht.tx_annotation.csq,
            'modifier': ht.tx_annotation.lof,
            'gene': ht.tx_annotation.symbol,
            'expressed': hl.case(missing_false=True).when(
                ht.tx_annotation.mean_expression >= 0.9, 'high').when(
                ht.tx_annotation.mean_expression > 0.1, 'medium').when(
                hl.is_defined(ht.tx_annotation.mean_expression), 'low').default('missing'),
            'coverage': ht.exome_coverage
        }
    else:
        groupings = {
            'annotation': ht.transcript_consequences.most_severe_consequence,
            'modifier': hl.case()
                .when(hl.is_defined(ht.transcript_consequences.lof),
                      ht.transcript_consequences.lof)
                .when(hl.is_defined(ht.transcript_consequences.polyphen_prediction),
                      ht.transcript_consequences.polyphen_prediction)
                .default('None'),
            'transcript': ht.transcript_consequences.transcript_id,
            'gene': ht.transcript_consequences.gene_symbol,
            'canonical': hl.or_else(ht.transcript_consequences.canonical == 1, False),
            'coverage': ht.exome_coverage
        }
        if custom_model == 'splice_region':
            groupings['distance_splice'] = ht.transcript_consequences
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
    ht = ht.annotate(high_coverage_proportion_observed=ht.observed_variants / ht.possible_variants)
    return ht.aggregate(hl.agg.group_by(ht.cpg,
                                        hl.agg.linreg(ht.high_coverage_proportion_observed, [1, ht.mu_snp],
                                                      weight=ht.possible_variants if weighted else None)
                                       ).map_values(lambda x: x.beta))


def build_plateau_models_pop(ht: hl.Table, weighted: bool = False) -> Dict[str, Tuple[float, float]]:
    """
    Calibrates high coverage model (returns intercept and slope)
    """
    pop_lengths = get_all_pop_lengths(ht)
    agg_expr = {
        pop: [hl.agg.group_by(ht.cpg,
                             hl.agg.linreg(ht[f'observed_{pop}'][i] / ht.possible_variants, [1, ht.mu_snp],
                                           weight=ht.possible_variants if weighted else None)
                             ).map_values(lambda x: x.beta) for i in range(length)]
        for length, pop in pop_lengths
    }
    agg_expr['total'] = hl.agg.group_by(ht.cpg,
                                        hl.agg.linreg(ht.observed_variants / ht.possible_variants, [1, ht.mu_snp],
                                                      weight=ht.possible_variants if weighted else None)
                                        ).map_values(lambda x: x.beta)
    return ht.aggregate(hl.struct(**agg_expr))


def get_all_pop_lengths(ht, prefix: str = 'observed_', pops: List[str] = POPS, skip_assertion: bool = False):
    ds_lengths = ht.aggregate([hl.agg.min(hl.len(ht[f'{prefix}{pop}'])) for pop in pops])
    # temp_ht = ht.take(1)[0]
    # ds_lengths = [len(temp_ht[f'{prefix}{pop}']) for pop in pops]
    pop_lengths = list(zip(ds_lengths, pops))
    print('Found: ', pop_lengths)
    if not skip_assertion:
        assert ht.all(hl.all(lambda f: f, [hl.len(ht[f'{prefix}{pop}']) == length for length, pop in pop_lengths]))
    return pop_lengths


def get_downsamplings(ht):
    freq_meta = ht.freq_meta.collect()[0]
    downsamplings = [(i, int(x.get('downsampling'))) for i, x in enumerate(freq_meta)
                     if x.get('group') == 'adj' and x.get('pop') == 'global'
                     and x.get('downsampling') is not None]
    return downsamplings


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
    """
    Compute pLI score using the number of observed and expected variants.
    
    The output Table will include the following annotations:
        - pLI
        - pNull
        - pRec

    :param ht: Input Table.
    :param obs: The number of observed variants on each gene or transcript.
    :param exp: The number of expected variants on each gene or transcript.
    :return: StructExpression with pLI score.
    """
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
                           prefix: str = 'oe', alpha: float = 0.05, select_only_ci_metrics: bool = True) -> hl.Table:
    """
    Indicate the bounds of confidence interval around the observed:expected ratio.
    
    For a given pair of observed and expected values, the function computes the density of the Poisson distribution
    with fixed k (the observed number of variants) over a range of lambda values, which are given by the expected 
    number of variants times a varying parameter ranging between 0 and 2. The cumulative density function of this 
    function is computed and the value of the varying parameter is extracted at points corresponding to 5% and 95% 
    to indicate the bounds of the confidence interval.
    
    Function will have following annotations in the output Table in addition to keys:
        - oe_lof_lower - the lower bound of confidence interval
        - oe_lof_upper - the upper bound of confidence interval

    :param ht: Input Table with the observed and expected variant counts for LoF, missense, and synonymous variants.
    :param obs: observed variant counts of LoF, missense, or synonymous variants.
    :param exp: expected variant counts of LoF, missense, or synonymous variants.
    :param prefix: Prefix of upper and lower bounds, defaults to 'oe'.
    :param alpha: The lower bound of , defaults to 0.05.
    :param select_only_ci_metrics: Whether to return only upper and lower bounds, defaults to True.
    :return: Table with the bounds of confidence interval.
    """
    ht = ht.annotate(_obs=obs, _exp=exp)
    oe_ht = ht.annotate(_range=hl.range(0, 2000).map(lambda x: hl.float64(x) / 1000))
    oe_ht = oe_ht.annotate(_range_dpois=oe_ht._range.map(lambda x: hl.dpois(oe_ht._obs, oe_ht._exp * x)))

    oe_ht = oe_ht.transmute(_cumulative_dpois=hl.cumulative_sum(oe_ht._range_dpois))
    max_cumulative_dpois = oe_ht._cumulative_dpois[-1]
    oe_ht = oe_ht.transmute(_norm_dpois=oe_ht._cumulative_dpois.map(lambda x: x / max_cumulative_dpois))
    oe_ht = oe_ht.transmute(
        _lower_idx=hl.argmax(oe_ht._norm_dpois.map(lambda x: hl.or_missing(x < alpha, x))),
        _upper_idx=hl.argmin(oe_ht._norm_dpois.map(lambda x: hl.or_missing(x > 1 - alpha, x)))
    )
    oe_ht = oe_ht.transmute(**{
        f'{prefix}_lower': hl.cond(oe_ht._obs > 0, oe_ht._range[oe_ht._lower_idx], 0),
        f'{prefix}_upper': oe_ht._range[oe_ht._upper_idx]
    })
    if select_only_ci_metrics:
        return oe_ht.select(f'{prefix}_lower', f'{prefix}_upper')
    else:
        return oe_ht.drop('_exp')


def calculate_z(input_ht: hl.Table, obs: hl.expr.NumericExpression, exp: hl.expr.NumericExpression, output: str = 'z_raw') -> hl.Table:
    """
    Compute the signed raw z score using observed and expected variant counts.
    
    The raw z scores are positive when the transcript had fewer variants than expected, and are negative when transcripts had more variants than expected.
    
    Function will have following annotations in the output Table in addition to keys:
        - `output` - the raw z score

    :param input_ht: Input Table.
    :param obs: Observed variant counts.
    :param exp: Expected variant counts.
    :param output: The column name for raw z score, defaults to 'z_raw'.
    :return: A Table with raw z score.
    """
    ht = input_ht.select(_obs=obs, _exp=exp)
    ht = ht.annotate(_chisq=(ht._obs - ht._exp) ** 2 / ht._exp)
    return ht.select(**{output: hl.sqrt(ht._chisq) * hl.cond(ht._obs > ht._exp, -1, 1)})


def calculate_all_z_scores(ht: hl.Table) -> hl.Table:
    """
    Calculate z scores for synomynous variants, missense variants, and LoF variants.
    
    z score = {variant_type}_z_raw / {variant_type}_sd (variant_type could be syn, mis, or lof)
    
    Function will add the following annotations to output Table:
        - syn_sd (global) - standard deviation of synonymous variants raw z score
        - mis_sd (global) - standard deviation of missense varinats raw z score
        - lof_sd (global) - standard deviation of LoF variants raw z score
        - constraint_flag - flags to determine if the variant is 'no_variants', 'no_exp_syn', 'no_exp_mis', 'no_exp_lof', 'syn_outlier', 'mis_too_many', or 'lof_too_many'
        - syn_z - z score of synonymous variants
        - mis_z - z score of missense varinats
        - lof_z - z score of LoF variants

    :param ht: Input table with observed and expected variant counts for synomynous variants, missense variants, and LoF variants.
    :return: Table with z scores.
    """
    ht = ht.annotate(**calculate_z(ht, ht.obs_syn, ht.exp_syn, 'syn_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_mis, ht.exp_mis, 'mis_z_raw')[ht.key])
    ht = ht.annotate(**calculate_z(ht, ht.obs_lof, ht.exp_lof, 'lof_z_raw')[ht.key])
    reasons = hl.empty_set(hl.tstr)
    reasons = hl.cond(hl.or_else(ht.obs_syn, 0) + hl.or_else(ht.obs_mis, 0) + hl.or_else(ht.obs_lof, 0) == 0, reasons.add('no_variants'), reasons)
    reasons = hl.cond(ht.exp_syn > 0, reasons, reasons.add('no_exp_syn'), missing_false=True)
    reasons = hl.cond(ht.exp_mis > 0, reasons, reasons.add('no_exp_mis'), missing_false=True)
    reasons = hl.cond(ht.exp_lof > 0, reasons, reasons.add('no_exp_lof'), missing_false=True)
    reasons = hl.cond(hl.abs(ht.syn_z_raw) > 5, reasons.add('syn_outlier'), reasons, missing_false=True)
    reasons = hl.cond(ht.mis_z_raw < -5, reasons.add('mis_too_many'), reasons, missing_false=True)
    reasons = hl.cond(ht.lof_z_raw < -5, reasons.add('lof_too_many'), reasons, missing_false=True)
    ht = ht.annotate(constraint_flag=reasons)
    sds = ht.aggregate(hl.struct(
        syn_sd=hl.agg.filter(
            ~ht.constraint_flag.contains('no_variants') &
            ~ht.constraint_flag.contains('syn_outlier') &
            ~ht.constraint_flag.contains('no_exp_syn') &
            hl.is_defined(ht.syn_z_raw),
            hl.agg.stats(ht.syn_z_raw)).stdev,
        mis_sd=hl.agg.filter(
            ~ht.constraint_flag.contains('no_variants') &
            ~ht.constraint_flag.contains('mis_outlier') &
            ~ht.constraint_flag.contains('no_exp_mis') &
            hl.is_defined(ht.mis_z_raw) & (ht.mis_z_raw < 0),
            hl.agg.explode(lambda x: hl.agg.stats(x), [ht.mis_z_raw, -ht.mis_z_raw])
        ).stdev,
        lof_sd=hl.agg.filter(
            ~ht.constraint_flag.contains('no_variants') &
            ~ht.constraint_flag.contains('lof_outlier') &
            ~ht.constraint_flag.contains('no_exp_lof') &
            hl.is_defined(ht.lof_z_raw) & (ht.lof_z_raw < 0),
            hl.agg.explode(lambda x: hl.agg.stats(x), [ht.lof_z_raw, -ht.lof_z_raw])
        ).stdev
    ))
    print(sds)
    ht = ht.annotate_globals(**sds)
    return ht.transmute(syn_z=ht.syn_z_raw / sds.syn_sd,
                        mis_z=ht.mis_z_raw / sds.mis_sd,
                        lof_z=ht.lof_z_raw / sds.lof_sd)
