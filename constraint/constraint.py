#!/usr/bin/env python3

__author__ = 'konradk'

from constraint_utils import *


def main(args):
    hl.init(log='/constraint.log')

    if args.pre_process_data:
        # import_fasta()
        # vep_context_mt()
        # split_context_mt(unsplit_context_mt_path, {'exomes': coverage_ht_path('exomes'), 'genomes': coverage_ht_path('genomes')},
        #                  methylation_sites_mt_path(), context_mt_path, args.overwrite)
        pre_process_data(annotations_ht_path('genomes', 'frequencies'), annotations_ht_path('genomes', 'rf'),
                         context_mt_path, processed_genomes_ht_path, args.overwrite)
        pre_process_data(annotations_ht_path('exomes', 'frequencies'), annotations_ht_path('exomes', 'rf'),
                         context_mt_path, processed_exomes_ht_path, args.overwrite)

    full_context_ht = prepare_ht(hl.read_matrix_table(context_mt_path), args.trimers)  # .rows()
    full_genome_ht = prepare_ht(hl.read_matrix_table(processed_genomes_ht_path), args.trimers)
    full_exome_ht = prepare_ht(hl.read_matrix_table(processed_exomes_ht_path), args.trimers)

    context_ht = full_context_ht.filter_rows(full_context_ht.locus.in_autosome_or_par())
    genome_ht = full_genome_ht.filter_rows(full_genome_ht.locus.in_autosome_or_par())
    exome_ht = full_exome_ht.filter_rows(full_exome_ht.locus.in_autosome_or_par())

    if args.calculate_mutation_rate:
        raw_mu_ht = calculate_mu_by_downsampling(genome_ht, full_context_ht, recalculate_all_possible_summary=True,
                                                 recalculate_all_possible_summary_unfiltered=False)
        raw_mu_ht.write(mutation_rate_ht_path, overwrite=args.overwrite)
        # This command was run for Nicky for genome mutation rates
        # raw_mu_ht = calculate_mu_by_downsampling(genome_ht, full_context_ht,
        #                                          summary_file=f'gs://konradk/tmp/all_possible_counts_by_context_remove_common_ordinary_trimers_{args.trimers}.pckl',
        #                                          recalculate_all_possible_summary=True, remove_common_downsampled=False, remove_common_ordinary=True)
        # raw_mu_ht.write(f'{root}/exploratory/mutation_rate_methylation_remove_common_ordinary_trimers_{args.trimers}.ht', overwrite=args.overwrite)
        hl.read_table(mutation_rate_ht_path).select('mu_snp').export(mutation_rate_ht_path.replace('.ht', '.txt.bgz'))
        send_message(args.slack_channel, 'Mutation rate calculated!')

    mutation_ht = hl.read_table(mutation_rate_ht_path).select('mu_snp')

    context_x_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('X')])
    context_x_ht = context_x_ht.filter_rows(context_x_ht.locus.in_x_nonpar())
    context_y_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('Y')])
    context_y_ht = context_y_ht.filter_rows(context_y_ht.locus.in_y_nonpar())

    exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
    exome_x_ht = exome_x_ht.filter_rows(exome_x_ht.locus.in_x_nonpar())
    exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
    exome_y_ht = exome_y_ht.filter_rows(exome_y_ht.locus.in_y_nonpar())

    if args.build_model:
        coverage_ht = get_proportion_observed_by_coverage(exome_ht, context_ht, mutation_ht, True)
        annotate_variant_types(coverage_ht).write(po_coverage_ht_path, overwrite=args.overwrite)

        coverage_x_ht = get_proportion_observed_by_coverage(exome_x_ht, context_x_ht, mutation_ht, True)
        annotate_variant_types(coverage_x_ht).write(po_coverage_ht_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)

        # TODO: consider 20X cutoff for Y
        coverage_y_ht = get_proportion_observed_by_coverage(exome_y_ht, context_y_ht, mutation_ht, True)
        annotate_variant_types(coverage_y_ht).write(po_coverage_ht_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)

        send_message(args.slack_channel, 'Coverage data calculated!')

    coverage_ht = hl.read_table(po_coverage_ht_path)
    coverage_x_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_x.ht'))
    coverage_y_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_y.ht'))

    coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True)
    _, plateau_x_models = build_models(coverage_x_ht, args.trimers, True)
    _, plateau_y_models = build_models(coverage_y_ht, args.trimers, True)

    po_output_path = po_ht_path.replace('.ht', f'_{args.model}.ht')
    output_path = constraint_ht_path.replace('.ht', f'_{args.model}.ht')
    if args.apply_model:
        get_proportion_observed(exome_ht, context_ht, mutation_ht, plateau_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model).write(po_output_path, overwrite=args.overwrite)
        hl.read_table(po_output_path).export(po_output_path.replace('.ht', '.txt.bgz'))

        get_proportion_observed(exome_x_ht, context_x_ht, mutation_ht, plateau_x_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model).write(po_output_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        hl.read_table(po_output_path.replace('.ht', '_x.ht')).export(po_output_path.replace('.ht', '_x.txt.bgz'))

        get_proportion_observed(exome_y_ht, context_y_ht, mutation_ht, plateau_y_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model).write(po_output_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)
        hl.read_table(po_output_path.replace('.ht', '_y.ht')).export(po_output_path.replace('.ht', '_y.txt.bgz'))

    if args.finalize:
        ht = hl.read_table(po_output_path).union(
            hl.read_table(po_output_path.replace('.ht', '_x.ht'))
        ).union(
            hl.read_table(po_output_path.replace('.ht', '_y.ht'))
        )
        if args.model != 'syn_canonical':
            ht = finalize_dataset(ht, keys=MODEL_KEYS[args.model])
        ht.write(output_path, args.overwrite)
        hl.read_table(output_path).export(output_path.replace('.ht', '.txt.bgz'))

        # Prepare release
        ht = hl.read_table(output_path)
        ht.select(
            obs_lof=ht.obs_lof_classic_hc, exp_lof=ht.exp_lof_classic_hc, oe_lof=ht.oe_lof_classic_hc,
            oe_lof_lower=ht.oe_classic_hc_lower, oe_lof_upper=ht.oe_classic_hc_upper,
            obs_mis=ht.obs_mis, exp_mis=ht.exp_mis, oe_mis=ht.oe_mis,
            oe_mis_lower=ht.oe_mis_lower, oe_mis_upper=ht.oe_mis_upper,
            obs_syn=ht.obs_syn, exp_syn=ht.exp_syn, oe_syn=ht.oe_syn,
            oe_syn_lower=ht.oe_syn_lower, oe_syn_upper=ht.oe_syn_upper,
            lof_z=ht.lof_z, mis_z=ht.mis_z, syn_z=ht.syn_z,
            pLI=ht.pLI_classic_hc, pRec=ht.pRec_classic_hc, pNull=ht.pNull_classic_hc, gene_issues=ht.constraint_flag
        ).select_globals().write(constraint_ht_path, overwrite=args.overwrite)
        ht = hl.read_table(constraint_ht_path)
        ht.export(constraint_ht_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
    parser.add_argument('--calculate_mutation_rate', help='Calculate mutation rate', action='store_true')
    parser.add_argument('--build_model', help='Calculate proportion observed by mu by coverage', action='store_true')
    parser.add_argument('--apply_model', help='Apply constraint model', action='store_true')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now)')
    parser.add_argument('--finalize', help='Combine autosomes, X, Y, and finalize', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel.split(','), main, args)
    else:
        main(args)
