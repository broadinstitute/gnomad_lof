from constraint_utils import *

root = 'gs://gnomad-public/papers/2019-flagship-lof/v1.1'
subdir = 'summary_results'
maps_ht_path = f'{root}/{subdir}/maps_plain_{{data_type}}.ht'
sfs_ht_path = f'{root}/{subdir}/sfs_{{data_type}}.ht'
loftee_maps_ht_path = f'{root}/{subdir}/maps_loftee_{{data_type}}.ht'
fifty_bp_maps_ht_path = f'{root}/{subdir}/maps_end_trunc_50bp_{{data_type}}.ht'
end_trunc_maps_ht_path = f'{root}/{subdir}/maps_end_trunc_gerp_{{data_type}}.ht'
loftee_assess_ht_path = f'{root}/{subdir}/freq_loftee_{{data_type}}.ht'
variants_per_sample_ht_path = f'{root}/{subdir}/variants_per_sample_{{data_type}}.ht'
observed_possible_ht_path = f'{root}/{subdir}/observed_possible_expanded_{{data_type}}.ht'
observed_possible_sites_ht_path = f'{root}/{subdir}/observed_possible_sites_{{data_type}}.txt'
indels_summary_ht_path = f'{root}/{subdir}/indels_summary_{{data_type}}.ht'
methylation_hist_file = f'{root}/{subdir}/methylation_hist.txt.bgz'


def get_worst_consequence_with_non_coding(ht):
    def get_worst_csq(csq_list: hl.expr.ArrayExpression, protein_coding: bool) -> hl.struct:
        lof = hl.null(hl.tstr)
        no_lof_flags = hl.null(hl.tbool)
        # lof_filters = hl.null(hl.tstr)
        # lof_flags = hl.null(hl.tstr)
        if protein_coding:
            all_lofs = csq_list.map(lambda x: x.lof)
            lof = hl.literal(['HC', 'OS', 'LC']).find(lambda x: all_lofs.contains(x))
            csq_list = hl.cond(hl.is_defined(lof), csq_list.filter(lambda x: x.lof == lof), csq_list)
            no_lof_flags = hl.or_missing(hl.is_defined(lof),
                                         csq_list.any(lambda x: (x.lof == lof) & hl.is_missing(x.lof_flags)))
            # lof_filters = hl.delimit(hl.set(csq_list.map(lambda x: x.lof_filter).filter(lambda x: hl.is_defined(x))), '|')
            # lof_flags = hl.delimit(hl.set(csq_list.map(lambda x: x.lof_flags).filter(lambda x: hl.is_defined(x))), '|')
        all_csq_terms = csq_list.flatmap(lambda x: x.consequence_terms)
        worst_csq = hl.literal(CSQ_ORDER).find(lambda x: all_csq_terms.contains(x))
        return hl.struct(worst_csq=worst_csq, protein_coding=protein_coding, lof=lof, no_lof_flags=no_lof_flags,
                         # lof_filters=lof_filters, lof_flags=lof_flags
                         )

    protein_coding = ht.vep.transcript_consequences.filter(lambda x: x.biotype == 'protein_coding')
    return ht.annotate(**hl.case(missing_false=True)
                       .when(hl.len(protein_coding) > 0, get_worst_csq(protein_coding, True))
                       .when(hl.len(ht.vep.transcript_consequences) > 0, get_worst_csq(ht.vep.transcript_consequences, False))
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


def get_curation_data():
    ht = hl.import_table('gs://gnomad-public/papers/2019-flagship-lof/v1.1/Final_gnomad_LOF_curation.csv',
                         delimiter=',', impute=True, quote='"')
    ht = ht.annotate(curated=(ht.verdict == 'LoF') | (ht.verdict == 'likely_LoF'))
    var = ht.variant_id.split('-')
    return ht.key_by(locus=hl.locus(var[0], hl.int(var[1])), alleles=var[2:4])


def main(args):
    context_ht = hl.read_table(context_ht_path())
    mutation_ht = hl.read_table(mutation_rate_ht_path)

    for data_type, sample_sizes in data_type_sex_counts.items():
        print(f'Running {data_type}...')
        ht = get_gnomad_public_data(data_type)
        coverage_ht = hl.read_table(coverage_ht_path(data_type))

        ht = ht.annotate(coverage=coverage_ht[ht.locus].median)
        ht = ht.filter((hl.len(ht.filters) == 0))  # & get_an_adj_criteria(ht, sample_sizes))
        ht = filter_vep_to_canonical_transcripts(ht)
        ht = get_worst_consequence_with_non_coding(ht)

        if args.run_indels:
            print(f'Running indels for {data_type}...')
            indel_ht = ht.filter(hl.is_indel(ht.alleles[0], ht.alleles[1]))

            indels = indel_ht.group_by(indel_ht.worst_csq, indel_ht.coverage,
                                       indel_length=hl.len(indel_ht.alleles[1]) - hl.len(indel_ht.alleles[0])).partition_hint(100).aggregate(
                singletons=hl.agg.count_where(indel_ht.freq[0].AC == 1),
                observed=hl.agg.count(),
                downsampling_array=hl.agg.array_sum(indel_ht.freq.map(lambda x: x.AC > 0)),
                singleton_downsampling_array=hl.agg.array_sum(indel_ht.freq.map(lambda x: x.AC == 1)))
            total_real_estate = coverage_ht.group_by(coverage_ht.median).partition_hint(100).aggregate(real_estate=hl.agg.count())
            indels.annotate(coverage_real_estate=total_real_estate[indels.coverage].real_estate).write(
                indels_summary_ht_path.format(data_type=data_type), args.overwrite)
            indels = hl.read_table(indels_summary_ht_path.format(data_type=data_type))
            indels.export(
                indels_summary_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
            )
            explode_downsamplings(indels.annotate(possible=indels.coverage_real_estate),
                                  sum(sample_sizes.values())).export(
                indels_summary_ht_path.format(data_type=data_type).replace('.ht', '.downsampling.txt.bgz'))

        if args.run_sfs:
            criteria = hl.case(missing_false=True).when(ht.freq[0].AC == 1, 'Singleton')
            if data_type == 'genomes':
                criteria = criteria.when(ht.freq[0].AF < 1e-3, '< 0.1%')
            else:
                criteria = (criteria.when(ht.freq[0].AC == 2, 'Doubleton')
                            .when(ht.freq[0].AC <= 5, 'AC 3-5')
                            .when(ht.freq[0].AF < 1e-4, '< 0.01%')
                            .when(ht.freq[0].AF < 1e-3, '0.01% - 0.1%'))
            sfs_ht = ht.annotate(
                freq_bin=criteria
                    .when(ht.freq[0].AF < 1e-2, '0.1% - 1%')
                    .when(ht.freq[0].AF < 1e-1, '1% - 10%')
                    .default('>10%')
            )
            sfs_ht.group_by(
                sfs_ht.freq_bin, sfs_ht.worst_csq, snp=hl.is_snp(sfs_ht.alleles[0], sfs_ht.alleles[1])
            ).aggregate(total=hl.agg.count()).write(sfs_ht_path.format(data_type=data_type), overwrite=args.overwrite)
            hl.read_table(sfs_ht_path.format(data_type=data_type)).export(
                sfs_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
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

        if args.run_obs_poss_sites:
            print(f'Running observed/possible sites for {data_type}...')

            initial_alleles = hl.read_table(annotations_ht_path(data_type, 'allele_data'))
            initial_alleles = initial_alleles.filter(initial_alleles.a_index == 1)

            sites_snp_ht = snp_ht.filter(hl.is_defined(initial_alleles[snp_ht.key]))

            observed = sites_snp_ht.aggregate(hl.agg.counter(sites_snp_ht.coverage))
            possible = context_ht.aggregate(hl.agg.counter(context_ht.coverage[data_type].median))
            for coverage, n_sites in observed.items():
                if coverage not in possible:
                    raise ValueError(f'{coverage} not found in possible dict')

            temp_f = hl.utils.new_local_temp_file()
            print(f'Writing to {temp_f}...')
            with open(temp_f, 'w') as f:
                f.write('coverage\tn_sites_observed\tn_sites_possible\n')
                for coverage, n_sites_possible in possible.items():
                    n_sites_observed = observed.get(coverage, 0)
                    f.write(f'{coverage}\t{n_sites_observed}\t{n_sites_possible / 3}\n')
            hl.hadoop_copy(f'file://{temp_f}', observed_possible_sites_ht_path.format(data_type=data_type))

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

        if args.run_variants_per_sample:
            ht = get_gnomad_public_data(data_type)

            ht = ht.filter((hl.len(ht.filters) == 0))  # & get_an_adj_criteria(ht, sample_sizes))
            ht = filter_vep_to_canonical_transcripts(ht)
            ht = get_worst_consequence_with_non_coding(ht)
            ht = filter_low_conf_regions(ht, annotate_instead_of_filter=True)
            ht = ht.annotate(curated=get_curation_data()[ht.key].curated)

            def build_criteria(ht: hl.Table, data_type: str, index: int = 0):
                criteria = hl.case(missing_false=True).when(
                    ht.freq[index].AC == 0, 'Not found').when(
                    ht.freq[index].AC == 1, 'Singleton').when(
                    ht.freq[index].AC == 2, 'Doubleton')
                if data_type == 'genomes':
                    criteria = criteria.when(ht.freq[index].AF < 1e-3, '< 0.1%')
                else:
                    criteria = (criteria
                                .when(ht.freq[index].AC <= 5, 'AC 3 - 5')
                                .when(ht.freq[index].AF < 1e-4, 'AC 6 - 0.01%')
                                .when(ht.freq[index].AF < 1e-3, '0.01% - 0.1%'))
                return (criteria
                        .when(ht.freq[index].AF < 1e-2, '0.1% - 1%')
                        .when(ht.freq[index].AF < 1e-1, '1% - 10%')
                        .when(ht.freq[index].AF > 0.95, '>95%')
                        .default('10% - 95%'))

            pops = list(map(lambda x: x.lower(), EXOME_POPS if data_type == 'exomes' else GENOME_POPS))
            pops = [(pop, hl.eval(ht.freq_index_dict[f'gnomad_{pop}'])) for pop in pops] + [('global', 0)]
            ht = ht.group_by('worst_csq', 'lof', 'no_lof_flags', 'protein_coding', 'curated', **ht.regions).aggregate(
                pop_bin_sums=[(pop, hl.agg.group_by(build_criteria(ht, data_type, index),
                                                    [hl.agg.sum(ht.freq[index].AC),
                                                     hl.agg.sum(ht.freq[index].homozygote_count)]))
                              for pop, index in pops]
            )
            ht = ht.explode('pop_bin_sums')
            ht = ht.transmute(pop=ht.pop_bin_sums[0], bin_sums=hl.array(ht.pop_bin_sums[1]))
            ht = ht.explode('bin_sums')
            ht = ht.transmute(bin=ht.bin_sums[0], total=ht.bin_sums[1][0], total_hom=ht.bin_sums[1][1])
            # # Alternative to above code:
            # pop_freq_mapping = {f'bin_{pop}': build_criteria(ht, data_type, index) for pop, index in pops}
            # ht = ht.annotate(**pop_freq_mapping)
            # ht = ht.group_by(*list(pop_freq_mapping), 'worst_csq', 'lof', 'no_lof_flags', 'protein_coding').aggregate(
            #     **{f'count_{pop}': hl.agg.sum(ht.freq[index].AC) for pop, index in pops},
            #     ** {f'hom_{pop}': hl.agg.sum(ht.freq[index].homozygote_count) for pop, index in pops}
            # )
            ht.write(variants_per_sample_ht_path.format(data_type=data_type), overwrite=args.overwrite)
            hl.read_table(variants_per_sample_ht_path.format(data_type=data_type)).export(
                variants_per_sample_ht_path.format(data_type=data_type).replace('.ht', '.txt.bgz')
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
    parser.add_argument('--run_sfs', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_maps', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_loftee_maps', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_end_trunc_maps', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_obs_poss', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_obs_poss_sites', help='Overwrite everything', action='store_true')
    parser.add_argument('--run_variants_per_sample', help='Overwrite everything', action='store_true')
    parser.add_argument('--assess_loftee', help='Overwrite everything', action='store_true')
    parser.add_argument('--methylation_hist', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)


