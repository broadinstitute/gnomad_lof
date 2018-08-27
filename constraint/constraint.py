#!/usr/bin/env python3

__author__ = 'konradk'

from constraint_utils import *


def main(args):
    hl.init(log='/constraint.log')

    if args.pre_process_data:
        # import_fasta()
        # split_context_mt(context_mt_path, {'exomes': coverage_ht_path('exomes'), 'genomes': coverage_ht_path('genomes')},
        #                  methylation_sites_mt_path(), split_context_mt_path, args.overwrite)
        # mt = hl.read_matrix_table(split_context_mt_path)
        # ht = hl.read_table(coverage_ht_path('exomes'))
        # mt.annotate_rows(coverage=mt.coverage.annotate(exomes=ht[mt.locus])).write(new_split_context_mt_path)
        pre_process_data(annotations_ht_path('genomes', 'frequencies'), annotations_ht_path('genomes', 'rf'),
                         split_context_mt_path, processed_genomes_ht_path, args.overwrite)
        pre_process_data(annotations_ht_path('exomes', 'frequencies'), annotations_ht_path('exomes', 'rf'),
                         split_context_mt_path, processed_exomes_ht_path, args.overwrite)

    full_context_ht = prepare_ht(hl.read_matrix_table(split_context_mt_path), args.trimers)  # .rows()
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
        # raw_mu_ht = calculate_mu_by_downsampling(autosomes_genome_ht, full_context_ht,
        #                                          summary_file='gs://konradk/tmp/all_possible_counts_by_context_with_common_var.pckl',
        #                                          recalculate_all_possible_summary=True, remove_common=False,
        #                                          recalculate_all_possible_summary_unfiltered=False)
        # raw_mu_ht.write(f'{root}/exploratory/mutation_rate_methylation_with_common.ht', overwrite=args.overwrite)
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

    if args.get_mu_coverage:
        coverage_ht = get_proportion_observed_by_coverage(exome_ht, context_ht, mutation_ht, True)
        annotate_variant_types(coverage_ht).write(po_coverage_ht_path, overwrite=args.overwrite)

        coverage_x_ht = get_proportion_observed_by_coverage(exome_x_ht, context_x_ht, mutation_ht, False)
        annotate_variant_types(coverage_x_ht).write(po_coverage_ht_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)

        coverage_y_ht = get_proportion_observed_by_coverage(exome_y_ht, context_y_ht, mutation_ht, True)
        annotate_variant_types(coverage_y_ht).write(po_coverage_ht_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)

        send_message(args.slack_channel, 'Coverage data calculated!')

    coverage_ht = hl.read_table(po_coverage_ht_path)
    coverage_x_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_x.ht'))
    coverage_y_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_y.ht'))

    coverage_model, plateau_models = build_models(coverage_ht, args.trimers)
    _, plateau_x_models = build_models(coverage_x_ht, args.trimers)
    _, plateau_y_models = build_models(coverage_y_ht, args.trimers)

    if args.confirm_model:
        get_proportion_observed(exome_ht, context_ht, mutation_ht,
                                plateau_models, coverage_model, recompute_possible=True,
                                confirm_model_only=True).write(po_syn_ht_path, overwrite=args.overwrite)
        hl.read_table(po_syn_ht_path).export(po_syn_ht_path.replace('.ht', '.txt.bgz'))

        get_proportion_observed(exome_x_ht, context_x_ht, mutation_ht,
                                plateau_x_models, coverage_model, recompute_possible=True,
                                confirm_model_only=True).write(po_syn_ht_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        hl.read_table(po_syn_ht_path.replace('.ht', '_x.ht')).export(po_syn_ht_path.replace('.ht', '_x.txt.bgz'))

        get_proportion_observed(exome_y_ht, context_y_ht, mutation_ht,
                                plateau_y_models, coverage_model, recompute_possible=True,
                                confirm_model_only=True).write(po_syn_ht_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)
        hl.read_table(po_syn_ht_path.replace('.ht', '_y.ht')).export(po_syn_ht_path.replace('.ht', '_y.txt.bgz'))

    if args.build_full_model:
        get_proportion_observed(exome_ht, context_ht, mutation_ht, plateau_models,
                                coverage_model, recompute_possible=True).write(po_ht_path, overwrite=args.overwrite)
        hl.read_table(po_ht_path).export(po_ht_path.replace('.ht', '.txt.bgz'))

        get_proportion_observed(exome_x_ht, context_x_ht, mutation_ht, plateau_x_models,
                                coverage_model, recompute_possible=True).write(po_ht_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        hl.read_table(po_ht_path.replace('.ht', '_x.ht')).export(po_ht_path.replace('.ht', '_x.txt.bgz'))

        get_proportion_observed(exome_y_ht, context_y_ht, mutation_ht, plateau_y_models,
                                coverage_model, recompute_possible=True).write(po_ht_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)
        hl.read_table(po_ht_path.replace('.ht', '_y.ht')).export(po_ht_path.replace('.ht', '_y.txt.bgz'))

        ht = hl.read_table(po_ht_path).union(
            hl.read_table(po_ht_path.replace('.ht', '_x.ht'))
        ).union(
            hl.read_table(po_ht_path.replace('.ht', '_y.ht'))
        )
        ht = finalize_dataset(ht)
        # ht = ht.annotate(**oe_confidence_interval(ht)[ht.key])
        ht.write(constraint_ht_path, args.overwrite)
        hl.read_table(constraint_ht_path).export(constraint_ht_path.replace('.ht', '.txt.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
    parser.add_argument('--calculate_mutation_rate', help='Calculate mutation rate', action='store_true')
    parser.add_argument('--get_mu_coverage', help='Calculate proportion observed by mu by coverage', action='store_true')
    parser.add_argument('--confirm_model', help='Confirm model on synonymous variants on canonical transcripts', action='store_true')
    parser.add_argument('--build_full_model', help='Apply model to all transcripts and variant types', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
