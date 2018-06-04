#!/usr/bin/env python3

__author__ = 'konradk'

from constraint_utils import *


def main(args):
    hl.init(log='/constraint.log')

    if args.pre_process_data:
        import_fasta()
        split_context_mt(context_mt_path, {'exomes': coverage_ht_path('exomes'), 'genomes': coverage_ht_path('genomes')},
                         methylation_sites_mt_path(), split_context_mt_path, args.overwrite)
        pre_process_data(annotations_ht_path('genomes', 'frequencies'),
                         split_context_mt_path, processed_genomes_ht_path, args.overwrite)
        pre_process_data(annotations_ht_path('exomes', 'frequencies'),
                         split_context_mt_path, processed_exomes_ht_path, args.overwrite)

    full_context_ht = prepare_ht(hl.read_matrix_table(split_context_mt_path), args.trimers)  # .rows()
    full_genome_ht = prepare_ht(hl.read_matrix_table(processed_genomes_ht_path), args.trimers)
    full_exome_ht = prepare_ht(hl.read_table(processed_exomes_ht_path), args.trimers)

    autosomes_context_ht = filter_to_autosomes(full_context_ht)
    autosomes_genome_ht = filter_to_autosomes(full_genome_ht)
    autosomes_exome_ht = filter_to_autosomes(full_exome_ht)

    if args.calculate_mutation_rate:
        raw_mu_ht = calculate_mu_by_downsampling(autosomes_genome_ht, full_context_ht, recalculate_all_possible_summary=False,
                                                 recalculate_all_possible_summary_unfiltered=False)
        raw_mu_ht.write(mutation_rate_ht_path, overwrite=args.overwrite)
        hl.read_table(mutation_rate_ht_path).select('mu_snp').export(mutation_rate_ht_path.replace('.ht', '.txt.bgz'))
        send_message(args.slack_channel, 'Mutation rate calculated!')

    mutation_ht = hl.read_table(mutation_rate_ht_path).select('mu_snp')

    # pass_autosomes_exome_ht = filter_to_pass(autosomes_exome_ht)  # TODO
    rare_autosomes_exome_ht = filter_by_frequency(autosomes_exome_ht, 0.001, direction='below')

    if args.get_mu_coverage:
        coverage_ht = get_proportion_observed_by_coverage(rare_autosomes_exome_ht, autosomes_context_ht, mutation_ht, True)
        annotate_variant_types(coverage_ht).write(po_coverage_ht_path, overwrite=args.overwrite)
        send_message(args.slack_channel, 'Coverage data calculated!')

    if args.confirm_model or args.build_full_model:
        coverage_ht = hl.read_table(po_coverage_ht_path)

        keys = ['context', 'ref', 'alt', 'methylation_level', 'mutation_rate', 'cpg', 'transition', 'variant_type', 'variant_type_model']
        # po_coverage_rounded_kt = round_coverage(po_coverage_kt, keys + ['coverage']).key_by(keys)
        #
        # plateau_models = build_plateau_models(get_high_coverage_ht(po_coverage_rounded_kt, keys))
        # coverage_model = build_coverage_model(po_coverage_rounded_kt, keys)
        #
        # if args.confirm_model:
        #     get_proportion_observed(exome_vds, context_vds, mutation_kt, plateau_models, coverage_model, canonical=True, synonymous=True).write(po_syn_kt_path, overwrite=args.overwrite)
        #     hc.read_table(po_syn_kt_path).export(po_syn_kt_path.replace('.kt', '.txt.bgz'))
        #
        # if args.build_full_model:
        #     get_proportion_observed(exome_vds, context_vds, mutation_kt, plateau_models, coverage_model).write(po_kt_path, overwrite=args.overwrite)
        #     hc.read_table(po_kt_path).export(po_kt_path.replace('.kt', '.txt.bgz'))


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
