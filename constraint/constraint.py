#!/usr/bin/env python3

__author__ = 'konradk'

from constraint_utils import *


def main(args):
    hl.init(log='/constraint.log')
    hl.set_upload_email('konrad')
    hl.enable_pipeline_upload()
    hl._set_flags(cpp='1')

    if args.pre_process_data:
        import_fasta(raw_context_txt_path, raw_context_ht_path, args.overwrite)
        vep_context_ht(raw_context_ht_path, vep_context_ht_path, args.overwrite)
        split_context_mt(vep_context_ht_path, {'exomes': coverage_ht_path('exomes'), 'genomes': coverage_ht_path('genomes')},
                         methylation_sites_mt_path(), context_ht_path, args.overwrite)
        pre_process_data(get_gnomad_public_data('genomes'), context_ht_path, processed_genomes_ht_path, args.overwrite)
        pre_process_data(get_gnomad_public_data('exomes'), context_ht_path, processed_exomes_ht_path, args.overwrite)

    full_context_ht = prepare_ht(hl.read_table(context_ht_path), args.trimers)
    full_genome_ht = prepare_ht(hl.read_table(processed_genomes_ht_path), args.trimers)
    full_exome_ht = prepare_ht(hl.read_table(processed_exomes_ht_path), args.trimers)

    context_ht = full_context_ht.filter(full_context_ht.locus.in_autosome_or_par())
    genome_ht = full_genome_ht.filter(full_genome_ht.locus.in_autosome_or_par())
    exome_ht = full_exome_ht.filter(full_exome_ht.locus.in_autosome_or_par())

    if args.calculate_mutation_rate:
        raw_mu_ht = calculate_mu_by_downsampling(genome_ht, full_context_ht, recalculate_all_possible_summary=True,
                                                 recalculate_all_possible_summary_unfiltered=False)
        raw_mu_ht.write(mutation_rate_ht_path, overwrite=args.overwrite)
        hl.read_table(mutation_rate_ht_path).select('mu_snp').export(mutation_rate_ht_path.replace('.ht', '.txt.bgz'))
        send_message(args.slack_channel, 'Mutation rate calculated!')

    mutation_ht = hl.read_table(mutation_rate_ht_path).select('mu_snp')

    context_x_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('X')])
    context_x_ht = context_x_ht.filter(context_x_ht.locus.in_x_nonpar())
    context_y_ht = hl.filter_intervals(full_context_ht, [hl.parse_locus_interval('Y')])
    context_y_ht = context_y_ht.filter(context_y_ht.locus.in_y_nonpar())

    exome_x_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('X')])
    exome_x_ht = exome_x_ht.filter(exome_x_ht.locus.in_x_nonpar())
    exome_y_ht = hl.filter_intervals(full_exome_ht, [hl.parse_locus_interval('Y')])
    exome_y_ht = exome_y_ht.filter(exome_y_ht.locus.in_y_nonpar())

    global po_coverage_ht_path
    if args.dataset != 'gnomad':
        po_coverage_ht_path = po_coverage_ht_path.replace(root, root + f'/{args.dataset}')
    if args.skip_af_filter_upfront:
        po_coverage_ht_path = po_coverage_ht_path.replace(root, root + '/pop_specific')

    if args.build_model:
        coverage_ht = get_proportion_observed_by_coverage(exome_ht, context_ht, mutation_ht, True, args.dataset, not args.skip_af_filter_upfront)
        annotate_variant_types(coverage_ht).write(po_coverage_ht_path, overwrite=args.overwrite)
        hl.read_table(po_coverage_ht_path).export(po_coverage_ht_path.replace('.ht', '.txt.bgz'))

        coverage_x_ht = get_proportion_observed_by_coverage(exome_x_ht, context_x_ht, mutation_ht, True, args.dataset, not args.skip_af_filter_upfront)
        annotate_variant_types(coverage_x_ht).write(po_coverage_ht_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)

        # TODO: consider 20X cutoff for Y
        coverage_y_ht = get_proportion_observed_by_coverage(exome_y_ht, context_y_ht, mutation_ht, True, args.dataset, not args.skip_af_filter_upfront)
        annotate_variant_types(coverage_y_ht).write(po_coverage_ht_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)

        send_message(args.slack_channel, 'Coverage data calculated!')

    coverage_ht = hl.read_table(po_coverage_ht_path)
    coverage_x_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_x.ht'))
    coverage_y_ht = hl.read_table(po_coverage_ht_path.replace('.ht', '_y.ht'))

    coverage_model, plateau_models = build_models(coverage_ht, args.trimers, True)
    _, plateau_x_models = build_models(coverage_x_ht, args.trimers, True)
    _, plateau_y_models = build_models(coverage_y_ht, args.trimers, True)

    po_output_path = po_ht_path.format(subdir=args.model)
    output_path = raw_constraint_ht_path.format(subdir=args.model)
    if args.dataset != 'gnomad':
        po_output_path = po_output_path.replace(root, root + f'/{args.dataset}')
        output_path = output_path.replace(root, root + f'/{args.dataset}')
    if args.skip_af_filter_upfront:
        po_output_path = po_output_path.replace(root, root + '/pop_specific')
        output_path = output_path.replace(root, root + '/pop_specific')

    if args.apply_model:
        get_proportion_observed(exome_ht, context_ht, mutation_ht, plateau_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model, dataset=args.dataset,
                                impose_high_af_cutoff_upfront=not args.skip_af_filter_upfront
                                ).write(po_output_path, overwrite=args.overwrite)
        hl.read_table(po_output_path).export(po_output_path.replace('.ht', '.txt.bgz'))

        get_proportion_observed(exome_x_ht, context_x_ht, mutation_ht, plateau_x_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model, dataset=args.dataset,
                                impose_high_af_cutoff_upfront=not args.skip_af_filter_upfront
                                ).write(po_output_path.replace('.ht', '_x.ht'), overwrite=args.overwrite)
        hl.read_table(po_output_path.replace('.ht', '_x.ht')).export(po_output_path.replace('.ht', '_x.txt.bgz'))

        get_proportion_observed(exome_y_ht, context_y_ht, mutation_ht, plateau_y_models,
                                coverage_model, recompute_possible=True,
                                custom_model=args.model, dataset=args.dataset,
                                impose_high_af_cutoff_upfront=not args.skip_af_filter_upfront
                                ).write(po_output_path.replace('.ht', '_y.ht'), overwrite=args.overwrite)
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

    if args.prepare_release:
        ht = hl.read_table(output_path)
        var_types = ('lof', 'mis', 'syn')
        ht.select(
            *[f'{t}_{v}{ci}' for v in var_types
              for t, ci in zip(('obs', 'exp', 'oe', 'mu', 'oe', 'oe'),
                               ('', '', '', '', '_lower', '_upper'))],
            *[f'{v}_z' for v in var_types], 'pLI', 'pRec', 'pNull', gene_issues=ht.constraint_flag
        ).select_globals().write(f'{root}/constraint_final_{args.model}.ht', overwrite=args.overwrite)
        ht = hl.read_table(f'{root}/constraint_final_{args.model}.ht')
        ht.export(f'{root}/constraint_final_{args.model}.txt.bgz')

    hl.upload_log()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--trimers', help='Use trimers instead of heptamers', action='store_true')
    parser.add_argument('--pre_process_data', help='Pre-process data', action='store_true')
    parser.add_argument('--calculate_mutation_rate', help='Calculate mutation rate', action='store_true')
    parser.add_argument('--build_model', help='Calculate proportion observed by mu by coverage', action='store_true')
    parser.add_argument('--apply_model', help='Apply constraint model', action='store_true')
    parser.add_argument('--dataset', help='Which dataset to use (one of gnomad, non_neuro, non_cancer, controls)', default='gnomad')
    parser.add_argument('--model', help='Which model to apply (one of "standard", "syn_canonical", or "worst_csq" for now)', default='standard')
    parser.add_argument('--skip_af_filter_upfront', help='Skip AF filter up front (to be applied later to ensure that it is not affecting population-specific constraint): not generally recommended', action='store_true')
    parser.add_argument('--finalize', help='Combine autosomes, X, Y, and finalize', action='store_true')
    parser.add_argument('--prepare_release', help='Prepare release file', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack_and_upload_log(args.slack_channel.split(','), hl.upload_log, main, args)
        # try_slack(args.slack_channel.split(','), main, args)
    else:
        main(args)
