#!/usr/bin/env python3

__author__ = 'konradk'

from pprint import pprint

from constraint_utils import *

root = 'gs://gnomad-resources/lof_paper'
gene_lof_matrix_path = f'{root}/individual_level/full_gene_lof_matrix_{{}}.mt'
all_lof_metrics_path = f'{root}/full_lof_metrics_{{}}.ht'
homozygous_lof_mt_path = f'{root}/individual_level/all_homozygous_lofs.mt'
lofs_by_gene_ht_path = f'{root}/homozygous_lof_summary.ht'


def load_gtf_data():
    ht = hl.experimental.import_gtf('gs://konradk/gencode.v19.annotation.gtf.bgz', 'GRCh37', True, min_partitions=12)
    ht = ht.annotate(gene_id=ht.gene_id.split('\\.')[0],
                     transcript_id=ht.transcript_id.split('\\.')[0],
                     length=ht.interval.end.position - ht.interval.start.position + 1)
    genes = ht.filter(ht.feature == 'gene').select(
        'gene_id', 'gene_type', 'gene_name', 'length').rename({'length': 'gene_length'}).key_by('gene_id')
    coding_regions = ht.filter(ht.feature == 'CDS').select('gene_id', 'transcript_id', 'transcript_type', 'length', 'level')
    transcripts = coding_regions.group_by('transcript_id', 'transcript_type', 'gene_id',
                                          transcript_level=coding_regions.level).aggregate(
        cds_length=hl.agg.sum(coding_regions.length),
        num_coding_exons=hl.agg.count()
    ).key_by('transcript_id')
    return transcripts.annotate(**genes[transcripts.gene_id])


def load_gene_expression_data():
    ht = hl.import_table('gs://konradk/GTEx.v7.median_expression_per_tx_per_tissue.021018.tsv.bgz', impute=True, min_partitions=12)
    # Filtering out tissues with < 100 samples
    brain_tissues = [x[1] for x in list(ht.row.items()) if 'Brain' in x[0] and x[0] not in ('Brain_Spinalcord_cervicalc_1_', 'Brain_Substantianigra')]
    return ht.select(
        'transcript_id', 'gene_id',
        brain_expression=hl.mean(hl.filter(lambda k: ~hl.is_nan(k), brain_tissues), filter_missing=True)
    ).key_by('transcript_id')


def load_exac_pli_data():
    exac_pli = hl.import_table('gs://konradk/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt.gz', impute=True)
    return exac_pli.key_by(transcript=exac_pli.transcript.split('\\.')[0])


def load_gene_list_data():
    gene_lists = hl.import_table('gs://konradk/all_gene_lists.txt').key_by('gene').collect_by_key()
    return gene_lists.transmute(lists=gene_lists.values.map(lambda x: x.list))


def load_omim_data():
    omim = hl.import_table('gs://konradk/omim.use.tsv')
    omim = omim.annotate(gene=omim.genes.split('\|'))
    return omim.explode('gene').key_by('gene')


def add_rank(ht, field, ascending=True, total_genes=None, bins=10, defined_only=False):
    if total_genes is None:
        if defined_only:
            total_genes = ht.aggregate(hl.agg.count_where(hl.is_defined(ht[field])))
        else:
            total_genes = ht.count()
    rank_field = ht[field] if ascending else -ht[field]
    ht = ht.key_by(_rank=rank_field).add_index(f'{field}_rank').key_by()
    ht = ht.annotate(**{f'{field}_rank': hl.or_missing(
        hl.is_defined(ht._rank), ht[f'{field}_rank']
    )}).drop('_rank')
    return ht.annotate(**{
        f'{field}_bin': hl.int(ht[f'{field}_rank'] * bins / total_genes),
        f'{field}_bin_6': hl.int(ht[f'{field}_rank'] * 6 / total_genes)
    })


def select_primitives_from_ht(ht):
    return ht.select(**{x: v for x, v in ht.row_value.items() if
                        v.dtype in {hl.tstr, hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64, hl.tbool}})


def generate_gene_lof_matrix(mt: hl.MatrixTable, tx_ht: hl.Table, by_transcript: bool = False,
                             filter_an_adj: bool = False, common_var_filt: bool = False, pre_loftee: bool = False
                             ) -> hl.MatrixTable:
    filt_criteria = hl.len(mt.filters) == 0
    if filter_an_adj:
        filt_criteria &= get_an_adj_criteria(mt)
    if common_var_filt:
        filt_criteria &= mt.freq[0].AF < 0.05
    mt = mt.filter_rows(filt_criteria)
    if by_transcript:
        explode_field = 'transcript_consequences'
    else:
        mt = process_consequences(mt)
        explode_field = 'worst_csq_by_gene'
    if pre_loftee:
        lof_cats = hl.literal({"splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"})
        criteria = lambda x: lof_cats.contains(add_most_severe_consequence_to_consequence(x).most_severe_consequence)
    else:
        criteria = lambda x: (x.lof == 'HC') & hl.is_missing(x.lof_flags)
    lof_csqs = mt.vep[explode_field].filter(criteria)
    mt = mt.select_rows(mt.freq, lof_csqs=lof_csqs)
    mt = mt.explode_rows(mt.lof_csqs)
    annotation_expr = {
        'gene_id': mt.lof_csqs.gene_id,
        'gene': mt.lof_csqs.gene_symbol,
        'indel': hl.is_indel(mt.alleles[0], mt.alleles[1])
    }
    if by_transcript:
        annotation_expr['transcript_id'] = mt.lof_csqs.transcript_id
        annotation_expr['canonical'] = hl.is_defined(mt.lof_csqs.canonical)
    else:
        tx_annotation = annotate_tx_expression_data(mt, tx_ht, mt.lof_csqs).mean_expression
        annotation_expr['expressed'] = hl.case().when(
            tx_annotation >= 0.9, 'high').when(
            tx_annotation > 0.1, 'medium').when(
            hl.is_defined(tx_annotation), 'low').default('missing')
    mt = mt.annotate_rows(**annotation_expr)
    mt.describe()
    return mt.group_rows_by(*list(annotation_expr.keys())).aggregate_rows(
        n_sites=hl.agg.count(),
        n_sites_array=hl.agg.array_sum(mt.freq.map(lambda x: hl.int(x.AC > 0))),
        classic_caf=hl.agg.sum(mt.freq[0].AF),
        max_af=hl.agg.max(mt.freq[0].AF),
        classic_caf_array=hl.agg.array_sum(mt.freq.map(lambda x: x.AF)),
    ).aggregate_entries(
        num_homs=hl.agg.count_where(mt.GT.is_hom_var()),
        num_hets=hl.agg.count_where(mt.GT.is_het()),
        defined_sites=hl.agg.count_where(hl.is_defined(mt.GT))
    ).result()


def generate_gene_lof_summary(mt, collapse_indels: bool = False, by_transcript: bool = False):
    if collapse_indels:
        grouping = ['gene_id', 'gene']
        if by_transcript:
            grouping.extend(['transcript_id', 'canonical'])
        else:
            grouping.append('expressed')
        mt = mt.group_rows_by(*grouping).aggregate_rows(
            n_sites=hl.agg.sum(mt.n_sites),
            n_sites_array=hl.agg.array_sum(mt.n_sites_array),
            classic_caf=hl.agg.sum(mt.classic_caf),
            max_af=hl.agg.max(mt.max_af),
            classic_caf_array=hl.agg.array_sum(mt.classic_caf_array)
        ).aggregate_entries(
            num_homs=hl.agg.sum(mt.num_homs),
            num_hets=hl.agg.sum(mt.num_hets),
            defined_sites=hl.agg.sum(mt.defined_sites)
        ).result()
    ht = mt.annotate_rows(no_lofs=hl.agg.count_where((mt.defined_sites > 0) & (mt.num_homs + mt.num_hets == 0)),
                          obs_het_lof=hl.agg.count_where((mt.num_hets > 0) & (mt.num_homs == 0)),
                          obs_hom_lof=hl.agg.count_where(mt.num_homs > 0),
                          defined=hl.agg.count_where(mt.defined_sites > 0),
                          pop_no_lofs=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where((mt.defined_sites > 0) & (mt.num_homs + mt.num_hets == 0))),
                          pop_obs_het_lof=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where((mt.num_hets > 0) & (mt.num_homs == 0))),
                          pop_obs_hom_lof=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where(mt.num_homs > 0)),
                          pop_defined=hl.agg.group_by(
                              mt.meta.pop, hl.agg.count_where(mt.defined_sites > 0)),
                          ).rows()
    ht = ht.annotate(p=1 - hl.sqrt(hl.float64(ht.no_lofs) / ht.defined),
                     pop_p=hl.dict(hl.array(ht.pop_defined).map(
                         lambda x: (x[0], 1 - hl.sqrt(hl.float64(ht.pop_no_lofs.get(x[0])) / x[1])))))
    ht = ht.annotate(exp_hom_lof=ht.defined * ht.p * ht.p)
    return ht.annotate(oe=ht.obs_hom_lof / ht.exp_hom_lof)


def combine_lof_metrics(gene_lof_matrix, by_transcript: bool = False, pop_specific: bool = False):
    keys = ['gene']
    caf_keys = ['gene']
    caf_drop = ['gene_id', 'oe']
    constraint_path = raw_constraint_ht_path.format(subdir="standard" if by_transcript else "tx_annotation")
    if pop_specific: constraint_path = constraint_path.replace('standard/', 'pop_specific/standard/')
    ht = hl.read_table(constraint_path)
    if by_transcript:
        keys.append('transcript')
        caf_keys.append('transcript_id')
        caf_drop.append('canonical')
    else:
        keys.append('expressed')
        caf_keys.append('expressed')
    constraint_ht = ht.drop(*[x for x in ht.row if x.endswith('_classic') or x.endswith('_with_os')])
    constraint_ht = add_rank(constraint_ht, 'oe_lof_upper', defined_only=True).key_by(*keys)
    caf_ht = hl.read_table(gene_lof_matrix.replace('.mt', '.summary.ht'))
    mapping = caf_ht.freq_index_dict.collect()[0]
    caf_dict = {
        f'classic_caf_{pop}': caf_ht.classic_caf_array[mapping[f'gnomad_{pop}']] for pop in
        map(lambda x: x.lower(), EXOME_POPS)
    }
    caf_dict.update({
        f'p_{pop}': caf_ht.pop_p[pop] for pop in map(lambda x: x.lower(), EXOME_POPS)
    })
    caf_ht = caf_ht.annotate(**caf_dict)
    caf_ht = caf_ht.key_by(*caf_keys).drop(*caf_drop)
    exac_pli = load_exac_pli_data()
    gene_ht = load_gtf_data()
    expr_ht = load_gene_expression_data()
    if by_transcript:
        exac = exac_pli[constraint_ht.transcript]
        gene = gene_ht.drop('gene_name')[constraint_ht.transcript]
        expr_ht = expr_ht.drop('gene_id')
    else:
        exac = exac_pli.key_by('gene')[constraint_ht.gene]
        gene = gene_ht.key_by('gene_name').drop('transcript_id')[constraint_ht.gene]
        expr_ht = expr_ht.key_by('gene_id').drop('transcript_id')
    constraint_ht = constraint_ht.annotate(
        **caf_ht[constraint_ht.key], **gene,
        exac_pLI=exac.pLI, exac_obs_lof=exac.n_lof, exac_exp_lof=exac.exp_lof, exac_oe_lof=exac.n_lof / exac.exp_lof
    )
    # If CAF data is missing and LoFs are possible in the gene, set all CAF metrics to zero
    constraint_ht = constraint_ht.annotate(
        **{x: hl.cond(hl.is_missing(constraint_ht[x]) & (constraint_ht.possible_lof > 0),
                      0,
                      constraint_ht[x])
           for x in list(caf_ht.row)
           if x in constraint_ht.row and constraint_ht[x].dtype in {hl.tint32, hl.tfloat32, hl.tint64, hl.tfloat64}}
    )
    return constraint_ht.annotate(
        **expr_ht[constraint_ht.gene_id]
    ).annotate_globals(freq_meta=caf_ht.index_globals().freq_meta, freq_index_dict=caf_ht.index_globals().freq_index_dict)


def find_in_meta(downsampling, pop, freq_meta):
    for i, x in enumerate(freq_meta):
        if x.get('group') == 'adj' and \
                x.get('pop') == pop and \
                x.get('downsampling') == str(downsampling):
            return i
    print(f'failed on {downsampling} and {pop}')
    return None


def explode_downsamplings(ht, select_canonical=True):
    pop_lengths = get_all_pop_lengths(ht, prefix='exp_lof_')
    downsamplings = list(map(lambda x: x[1], get_downsamplings(ht)))
    freq_meta = ht.freq_meta.collect()[0]
    meta_locs = [(i, pop, downsamplings[i], find_in_meta(downsamplings[i], pop, freq_meta)) for count, pop in pop_lengths for i in range(count)]
    fields = ['data']
    if select_canonical: fields.insert(0, 'canonical')
    ht = ht.annotate(
        data=[
            hl.struct(
                pop=pop, downsampling=downsampling,
                exp_syn=ht[f'exp_syn_{pop}'][i], obs_syn=ht[f'obs_syn_{pop}'][i],
                exp_mis=ht[f'exp_mis_{pop}'][i], obs_mis=ht[f'obs_mis_{pop}'][i],
                exp_lof=ht[f'exp_lof_{pop}'][i], obs_lof=ht[f'obs_lof_{pop}'][i],
                n_sites=ht.n_sites_array[meta_loc],
                caf=ht.classic_caf_array[meta_loc]
            )
            for i, pop, downsampling, meta_loc in meta_locs] + [
            hl.struct(pop='global', downsampling=125748,
                      exp_syn=ht.exp_syn, obs_syn=ht.obs_syn,
                      exp_mis=ht.exp_mis, obs_mis=ht.obs_mis,
                      exp_lof=ht.exp_lof, obs_lof=ht.obs_lof,
                      n_sites=ht.n_sites, caf=ht.classic_caf)
        ]
    ).select(*fields)
    ht = ht.explode('data')
    return ht.transmute(**ht.data)


def compute_homozygous_lof(mt, tx_ht):
    mt = mt.filter_rows((hl.len(mt.filters) == 0) & (mt.freq[0].homozygote_count >= 1))
    mt = process_consequences(mt)
    mt = mt.explode_rows(mt.vep.worst_csq_by_gene)
    mt = mt.filter_rows((mt.vep.worst_csq_by_gene.lof == 'HC') & hl.is_missing(mt.vep.worst_csq_by_gene.lof_flags))
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_hom_var() &
                                   (mt.AD[1] / hl.sum(mt.AD) >= 0.8)
                                   ))
    mt = mt.annotate_rows(indel=hl.is_indel(mt.alleles[0], mt.alleles[1]),
                          tx_annotation=annotate_tx_expression_data(mt, tx_ht, mt.vep.worst_csq_by_gene))
    return mt.repartition(100)


def main(args):
    hl.init(log='/constraint.log')

    extension = 'by_transcript' if args.by_transcript else 'worst'
    if args.filter_an_adj:
        extension += '_an_adj'
    if args.no_loftee:
        extension += '_no_loftee'
    gene_lof_matrix = gene_lof_matrix_path.format(extension)
    if args.pop_specific:
        extension += '_pop_specific'
    all_lof_metrics = all_lof_metrics_path.format(extension)

    tx_ht = load_tx_expression_data(context=False)
    if args.calculate_gene_lof_matrix:
        mt = get_gnomad_data('exomes', adj=True, release_samples=True, release_annotations=True)
        mt = generate_gene_lof_matrix(mt, tx_ht, by_transcript=args.by_transcript,
                                      filter_an_adj=args.filter_an_adj, pre_loftee=args.no_loftee)
        mt.write(gene_lof_matrix, args.overwrite)
        send_message(args.slack_channel, 'Gene LoF matrix computed!')

    if args.export_gene_lof_matrix:
        mt = hl.read_matrix_table(gene_lof_matrix)
        ht = generate_gene_lof_summary(mt, not args.dont_collapse_indels, args.by_transcript)
        if args.dont_collapse_indels:
            extension += '_by_indel'
            gene_lof_matrix = gene_lof_matrix_path.format(extension)
        ht.write(gene_lof_matrix.replace('.mt', '.summary.ht'), args.overwrite)
        hl.read_table(gene_lof_matrix.replace('.mt', '.summary.ht')).drop('classic_caf_array').export(
            gene_lof_matrix.replace('.mt', '.summary.txt.bgz'))

    if args.combine_lof_metrics:
        ht = combine_lof_metrics(gene_lof_matrix, args.by_transcript, pop_specific=args.pop_specific)
        ht.write(all_lof_metrics, args.overwrite)
        # This file has been spot-checked. Of the canonical transcripts,
        # # 10.5% are missing CAF data due to having 0 HC no-flag LoFs
        # # # Many are actually 0, some e.g. single exon flag remove all LoFs
        # # # Some additional ones (SETD1B, TCF4, etc) are DQ'ed bc of AN_Adj
        # # 7.5% are missing ExAC pLI data (AC genes, low coverage, etc)
        # # 2.7% are missing LoF data (only zero LoF possible genes), confirmed by:
        # # # cht.filter(hl.is_missing(cht.obs_lof) & ~cht.gene_issues.contains('no_exp_lof')).show(10)
        # # 1 is missing oe_syn CIs: tiny Y chromosome gene with negative expected
        # Checked that number of exons is similar to ExAC (minor differences = UTRs)

    if args.export_combined_metrics:
        ht = hl.read_table(all_lof_metrics)
        ht = ht.annotate(constraint_flag=hl.delimit(ht.constraint_flag, '|'),
                         chromosome=ht.interval.start.contig,
                         start_position=ht.interval.start.position,
                         end_position=ht.interval.end.position)
        ht = select_primitives_from_ht(ht)
        ht.export(all_lof_metrics.replace('.ht', '.txt.bgz'))
        if args.by_transcript:
            ht = ht.filter(ht.canonical)
            ht = add_rank(ht, 'oe_lof_upper', defined_only=True).drop('canonical')
            ht.export(all_lof_metrics.replace('.ht', '_by_gene.txt.bgz'))

        ht = hl.read_table(all_lof_metrics)
        ht = explode_downsamplings(ht, args.by_transcript)
        ht.export(all_lof_metrics.replace('.ht', '.downsamplings.txt.bgz'))

    if args.compute_homozygous_lof:
        mt = get_gnomad_data('exomes', adj=True, non_refs_only=True,
                             release_samples=True, release_annotations=True)
        ht = compute_homozygous_lof(mt, tx_ht)
        ht.write(homozygous_lof_mt_path, args.overwrite)
        send_message(args.slack_channel, 'Homozygous LoFs computed!')

    if args.export_homozygous_lof:
        mt = hl.read_matrix_table(homozygous_lof_mt_path)
        vep_ht = hl.read_table(annotations_ht_path('exomes', 'vep_csq'))
        mt = mt.filter_rows(mt.tx_annotation.mean_expression > 0.1)
        mt = mt.annotate_rows(info=hl.struct(
            AC=mt.freq[0].AC, AN=mt.freq[0].AN, AF=mt.freq[0].AF, n_hom=mt.freq[0].homozygote_count,
            CSQ=vep_ht[mt.row_key].vep
        )).drop('is_missing')
        hl.export_vcf(mt, homozygous_lof_mt_path.replace('.mt', '.vcf.bgz'), metadata={
            'info': {'CSQ': {'Description': vep_ht.vep_csq_header.collect()[0]}}
        })
        ht = mt.rows()
        # ht.count()  # 3385 variants
        ht.select(
            'indel', AC=ht.info.AC, AN=ht.info.AN, AF=ht.info.AF, n_hom=ht.info.n_hom,
            variant_id=hl.delimit([ht.locus.contig, hl.str(ht.locus.position), ht.alleles[0], ht.alleles[1]], '-'),
            csq=ht.vep.worst_csq_by_gene
        ).flatten().export(homozygous_lof_mt_path.replace('.mt', '.txt.bgz'))

        ht = ht.group_by(ht.vep.worst_csq_by_gene.gene_symbol, ht.indel).aggregate(
            total=hl.agg.count(),
            total_twohit=hl.agg.count_where(ht.freq[0].homozygote_count >= 2)
        )
        ht.write(lofs_by_gene_ht_path, args.overwrite)
        hl.read_table(lofs_by_gene_ht_path).export(lofs_by_gene_ht_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--by_transcript', help='Use all transcripts instead of worst', action='store_true')
    parser.add_argument('--filter_an_adj', help='Filter variants with less than 80% callrate', action='store_true')
    parser.add_argument('--no_loftee', help='Do not filter with LOFTEE (only for comparison)', action='store_true')
    parser.add_argument('--dont_collapse_indels', help='Do not collapse indels when exporting Gene LoF matrix', action='store_true')
    parser.add_argument('--pop_specific', help='Use pop-specific matrix instead', action='store_true')
    parser.add_argument('--calculate_gene_lof_matrix', help='Calculate Gene LoF Matrix', action='store_true')
    parser.add_argument('--export_gene_lof_matrix', help='Export Gene LoF Matrix', action='store_true')
    parser.add_argument('--combine_lof_metrics', help='Combine all LoF metrics', action='store_true')
    parser.add_argument('--export_combined_metrics', help='Export combined LoF metrics', action='store_true')
    parser.add_argument('--compute_homozygous_lof', help='Compute Homozygous LoFs', action='store_true')
    parser.add_argument('--export_homozygous_lof', help='Export Homozygous LoF summary', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
