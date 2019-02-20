from constraint_utils import *

subdir = 'summary_results'
maps_ht_path = f'{root}/{subdir}/maps_plain_{{data_type}}.ht'
loftee_maps_ht_path = f'{root}/{subdir}/maps_loftee_{{data_type}}.ht'
fifty_bp_maps_ht_path = f'{root}/{subdir}/maps_end_trunc_50bp_{{data_type}}.ht'
end_trunc_maps_ht_path = f'{root}/{subdir}/maps_end_trunc_gerp_{{data_type}}.ht'
loftee_assess_ht_path = f'{root}/{subdir}/freq_loftee_{{data_type}}.ht'
observed_possible_ht_path = f'{root}/{subdir}/observed_possible_expanded_{{data_type}}.ht'
indels_summary_ht_path = f'{root}/{subdir}/indels_summary_{{data_type}}.ht'
methylation_hist_file = f'{root}/{subdir}/methylation_hist.txt.bgz'


def get_worst_consequence_with_non_coding(ht):
    def get_worst_csq(csq_list, check_biotype):
        all_csq_terms = csq_list.flatmap(lambda x: x.consequence_terms)
        worst_csq = hl.literal(CSQ_ORDER).find(lambda x: all_csq_terms.contains(x))
        biotype = csq_list.any(lambda x: (x.biotype == 'protein_coding') &
                                         x.consequence_terms.contains(worst_csq)) if check_biotype else False
        return hl.struct(worst_csq=worst_csq, protein_coding=biotype)

    return ht.annotate(**hl.case(missing_false=True)
                       .when(hl.len(ht.vep.transcript_consequences) > 0, get_worst_csq(ht.vep.transcript_consequences, True))
                       .when(hl.len(ht.vep.regulatory_feature_consequences) > 0, get_worst_csq(ht.vep.regulatory_feature_consequences, False))
                       .when(hl.len(ht.vep.motif_feature_consequences) > 0, get_worst_csq(ht.vep.motif_feature_consequences, False))
                       .default(get_worst_csq(ht.vep.intergenic_consequences, False)))


def explode_downsamplings(ht, full_sample_size):
    downsamplings = get_downsamplings(ht)

    ht = ht.transmute(
        data=[
            hl.struct(
                downsampling=downsamplings[i][1],
                singletons=ht.singleton_downsampling_array[downsamplings[i][0]],
                observed=ht.downsampling_array[downsamplings[i][0]],
                possible=ht.possible
            )
            for i in range(len(downsamplings))] + [
            hl.struct(downsampling=full_sample_size, singletons=ht.singletons,
                      observed=ht.observed, possible=ht.possible)
        ]
    )
    ht = ht.explode('data')
    return ht.transmute(**ht.data)


def parse_lof_info(ht: hl.Table, location: hl.expr.StringExpression = None):
    """
    Must be exploded
    """
    if location is None:
        location = ht.vep.transcript_consequences.lof_info
    ht = ht.annotate(lof_data=location.split(',').map(lambda x: x.split(':')))
    return ht.transmute(lof_data=hl.dict(ht.lof_data.map(lambda x: (x[0], hl.cond(hl.len(x) == 1, x[0], x[1])))))


def get_bin_boundaries(ht: hl.Table, feature: str = 'GERP_DIST', n_bins: int = 20, chrom: str = '10'):
    temp_ht = hl.filter_intervals(ht, [hl.parse_locus_interval(chrom)])
    temp_ht = temp_ht.filter(hl.is_defined(temp_ht.lof_data.get(feature)))
    data = temp_ht.aggregate(hl.agg.collect(hl.float(temp_ht.lof_data.get(feature))))
    data = sorted(data)
    l = int(len(data) / n_bins)
    bin_boundaries = [data[x - 1] for x in range(len(data), -1, -l)]
    bin_boundaries.append(min(data))
    return bin_boundaries


def main(args):
    context_ht = hl.read_table(context_ht_path)
    mutation_ht = hl.read_table(mutation_rate_ht_path)

    for data_type, sample_sizes in data_type_sex_counts.items():
        print(f'Running {data_type}...')
        ht = get_gnomad_public_data(data_type)
        coverage_ht = hl.read_table(coverage_ht_path(data_type))

        ht = ht.annotate(coverage=coverage_ht[ht.locus].median)
        ht = ht.filter((hl.len(ht.filters) == 0) & get_an_adj_criteria(ht, sample_sizes))
        ht = filter_vep_to_canonical_transcripts(ht)
        ht = get_worst_consequence_with_non_coding(ht)

        if args.run_indels:
            print(f'Running indels for {data_type}...')
            indel_ht = ht.filter(hl.is_indel(ht.alleles[0], ht.alleles[1]))

            indels = indel_ht.group_by(indel_ht.worst_csq, indel_ht.coverage).partition_hint(100).aggregate(
                singletons=hl.agg.count_where(indel_ht.freq[0].AC == 1),
                observed=hl.agg.count())
            total_real_estate = coverage_ht.group_by(coverage_ht.median).partition_hint(100).aggregate(real_estate=hl.agg.count())
            indels.annotate(coverage_real_estate=total_real_estate[indels.coverage].real_estate).write(
                indels_summary_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(indels_summary_ht_path.format(data_type=data_type)).export(
                indels_summary_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
            )

        context = context_ht[ht.key]
        snp_ht = prepare_ht(ht.annotate(context=context.context, methylation=context.methylation), True, False).persist()

        if args.run_maps:
            print(f'Running MAPS for {data_type}...')
            maps_ht = maps(snp_ht, mutation_ht, additional_grouping=['protein_coding', ])
            maps_ht.write(maps_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(maps_ht_path.format(data_type=data_type)).export(
                maps_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
            )

        if args.run_loftee_maps:
            exploded_snp_ht = snp_ht.transmute(transcript_consequences=snp_ht.vep.transcript_consequences)

            exploded_snp_ht = exploded_snp_ht.explode(exploded_snp_ht.transcript_consequences)
            exploded_snp_ht = parse_lof_info(exploded_snp_ht, exploded_snp_ht.transcript_consequences.lof_info)
            exploded_snp_ht = exploded_snp_ht.annotate(
                worst_csq=hl.literal(CSQ_ORDER).find(
                    lambda x: exploded_snp_ht.transcript_consequences.consequence_terms.contains(x)),
                lof=exploded_snp_ht.transcript_consequences.lof,
                lof_filter=exploded_snp_ht.transcript_consequences.lof_filter,
                lof_flags=exploded_snp_ht.transcript_consequences.lof_flags,
                lof_filter_simplified=hl.cond(exploded_snp_ht.transcript_consequences.lof_filter.contains(','),
                                              'MULTIPLE',
                                              exploded_snp_ht.transcript_consequences.lof_filter))

            maps_ht = maps(exploded_snp_ht, mutation_ht, ['lof', 'lof_filter', 'lof_flags'])
            maps_ht.write(loftee_maps_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(loftee_maps_ht_path.format(data_type=data_type)).export(
                loftee_maps_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
            )

        if args.run_end_trunc_maps:
            exploded_snp_ht = snp_ht.transmute(transcript_consequences=snp_ht.vep.transcript_consequences)

            exploded_snp_ht = exploded_snp_ht.explode(exploded_snp_ht.transcript_consequences)
            exploded_snp_ht = parse_lof_info(exploded_snp_ht, exploded_snp_ht.transcript_consequences.lof_info)

            gerp_bin_boundaries = get_bin_boundaries(exploded_snp_ht, 'GERP_DIST', n_bins=40)
            bp_bin_boundaries = get_bin_boundaries(exploded_snp_ht, 'BP_DIST', n_bins=40)
            percentile_breaks = list(
                map(lambda x: x / 10, range(9, 0, -1))) + list(
                map(lambda x: x / 200, range(18, -1, -1)))

            exploded_snp_ht = exploded_snp_ht.annotate(
                worst_csq=hl.literal(CSQ_ORDER).find(
                    lambda x: exploded_snp_ht.transcript_consequences.consequence_terms.contains(x)),
                fifty_bp_rule=exploded_snp_ht.lof_data.get('50_BP_RULE'),
                gerp=hl.literal(gerp_bin_boundaries).find(
                    lambda x: hl.float(exploded_snp_ht.lof_data.get('GERP_DIST')) >= x),
                bp=hl.literal(bp_bin_boundaries).find(lambda x: hl.float(
                    exploded_snp_ht.lof_data.get('BP_DIST')) >= x),
                percentile=hl.literal(percentile_breaks).find(lambda x: hl.float(
                    exploded_snp_ht.lof_data.get('PERCENTILE')) >= x)
            )
            # maps_ht = maps(exploded_snp_ht, mutation_ht, ['bp', ])
            # maps_ht = maps(exploded_snp_ht, mutation_ht, ['percentile', ])

            maps_ht = maps(exploded_snp_ht, mutation_ht, ['fifty_bp_rule', ])
            maps_ht.write(fifty_bp_maps_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(fifty_bp_maps_ht_path.format(data_type=data_type)).export(
                fifty_bp_maps_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
            )
            exploded_snp_ht = exploded_snp_ht.filter((exploded_snp_ht.worst_csq == 'synonymous_variant') |
                                                     (exploded_snp_ht.fifty_bp_rule == 'FAIL'))
            maps_ht = maps(exploded_snp_ht, mutation_ht, ['gerp', ])
            maps_ht.write(end_trunc_maps_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(end_trunc_maps_ht_path.format(data_type=data_type)).export(
                end_trunc_maps_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
            )

        if args.run_obs_poss:
            print(f'Running observed/possible for {data_type}...')

            context_ht = prepare_ht(context_ht, True, False)
            context_ht = filter_vep_to_canonical_transcripts(context_ht)
            context_ht = get_worst_consequence_with_non_coding(context_ht)

            observed = snp_ht.group_by(snp_ht.worst_csq, snp_ht.context, snp_ht.ref, snp_ht.alt,
                                       snp_ht.methylation_level, snp_ht.coverage).partition_hint(1000).aggregate(
                observed=hl.agg.count(),
                singletons=hl.agg.count_where(snp_ht.freq[0].AC == 1),
                downsampling_array=hl.agg.array_sum(snp_ht.freq.map(lambda x: x.AC > 0)),
                singleton_downsampling_array=hl.agg.array_sum(snp_ht.freq.map(lambda x: x.AC == 1))
            )
            possible = context_ht.group_by(context_ht.worst_csq, context_ht.context, context_ht.ref, context_ht.alt,
                                           context_ht.methylation_level, context_ht.coverage[data_type].median).partition_hint(1000).aggregate(
                possible=hl.agg.count())
            annotate_variant_types(observed.join(possible), False).write(
                observed_possible_ht_path.format(data_type=data_type), args.overwrite)

            explode_downsamplings(hl.read_table(observed_possible_ht_path.format(data_type=data_type)),
                                  sum(sample_sizes.values())).export(
                observed_possible_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz'))

        if data_type == 'exomes' and args.assess_loftee:
            ht = ht.annotate(vep=ht.vep.annotate(
                transcript_consequences=ht.vep.transcript_consequences.filter(lambda x: hl.is_defined(x.lof)))
            )
    
            ht = ht.annotate(
                freq_bin=hl.case(missing_false=True)
                    .when(ht.freq[0].AC == 1, 'Singleton')
                    .when(ht.freq[0].AC == 2, 'Doubleton')
                    .when(ht.freq[0].AF < 1e-4, '< 0.01%')
                    .when(ht.freq[0].AF < 1e-3, '0.01% - 0.1%')
                    .when(ht.freq[0].AF < 1e-2, '0.1% - 1%')
                    .when(ht.freq[0].AF < 1e-1, '1% - 10%')
                    .default('>10%'),
                fail_loftee=hl.any(lambda tc: (tc.lof == 'LC'), ht.vep.transcript_consequences),
                loftee_os=hl.any(lambda tc: (tc.lof == 'OS'), ht.vep.transcript_consequences),
                pass_loftee=hl.any(lambda tc: (tc.lof == 'HC'), ht.vep.transcript_consequences),
                pass_loftee_with_flags=hl.any(lambda tc: (tc.lof == 'HC') & hl.is_missing(tc.lof_flags),
                                              ht.vep.transcript_consequences)
            )
            ht = ht.group_by(ht.freq_bin, ht.worst_csq, ht.fail_loftee,
                             ht.loftee_os, ht.pass_loftee, ht.pass_loftee_with_flags
                             ).aggregate(n=hl.agg.count())
    
            clinvar_ht = hl.read_table(clinvar_ht_path)
            path_labels = hl.literal({'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'})
            clinvar_ht = clinvar_ht.filter(
                hl.or_missing(hl.len(clinvar_ht.info.CLNSIG) > 0,
                              path_labels.contains(clinvar_ht.info.CLNSIG[0])) &
                hl.is_missing(clinvar_ht.info.CLNSIGCONF)
            )
            clinvar_ht = filter_vep_to_canonical_transcripts(clinvar_ht)
            clinvar_ht = get_worst_consequence_with_non_coding(clinvar_ht)
    
            clinvar_ht = clinvar_ht.annotate(vep=clinvar_ht.vep.annotate(
                transcript_consequences=clinvar_ht.vep.transcript_consequences.filter(lambda x: hl.is_defined(x.lof)))
            )
            clinvar_ht = clinvar_ht.annotate(freq_bin=hl.delimit(clinvar_ht.info.CLNREVSTAT[:2], ',').replace('_', ' '),
                                             fail_loftee=hl.any(lambda tc: (tc.lof == 'LC'), clinvar_ht.vep.transcript_consequences),
                                             loftee_os=hl.any(lambda tc: (tc.lof == 'OS'), clinvar_ht.vep.transcript_consequences),
                                             pass_loftee=hl.any(lambda tc: (tc.lof == 'HC'), clinvar_ht.vep.transcript_consequences),
                                             pass_loftee_with_flags=hl.any(lambda tc: (tc.lof == 'HC') & hl.is_missing(tc.lof_flags),
                                                                           clinvar_ht.vep.transcript_consequences))
            clinvar_ht = clinvar_ht.group_by(
                clinvar_ht.freq_bin, clinvar_ht.worst_csq, clinvar_ht.fail_loftee,
                clinvar_ht.loftee_os, clinvar_ht.pass_loftee, clinvar_ht.pass_loftee_with_flags
            ).aggregate(n=hl.agg.count())
    
            ht = ht.union(clinvar_ht)
    
            ht.write(loftee_assess_ht_path.format(data_type=data_type), args.overwrite)
            hl.read_table(loftee_assess_ht_path.format(data_type=data_type)).export(
                loftee_assess_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
            )

    if args.methylation_hist:
        methylation_hist = context_ht.aggregate(hl.agg.hist(context_ht.methylation.MEAN, 0, 1, 40))
        data = list(zip(methylation_hist.bin_edges, methylation_hist.bin_freq))
        with hl.hadoop_open(methylation_hist_file, 'w') as f:
            f.write('edge\tfreq\n')
            for edge, freq in data:
                f.write(f'{edge}\t{freq}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_indels', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_maps', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_loftee_maps', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_end_trunc_maps', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_obs_poss', help='Overwrite everything', action='store_true')
    parser.add_argument('--assess_loftee', help='Overwrite everything', action='store_true')
    parser.add_argument('--methylation_hist', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


