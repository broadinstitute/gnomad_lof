#!/usr/bin/env python3

__author__ = 'konradk'

from constraint_utils import *

gene_lof_matrix = f'{root}/recessive/gene_matrix/full_gene_lof_matrix.mt'


def generate_gene_lof_matrix(mt: hl.MatrixTable) -> hl.MatrixTable:
    mt = process_consequences(mt)
    mt = mt.select_rows(gene_lofs=mt.vep.worst_csq_by_gene.values().filter(
        lambda x: (x.lof == 'HC') & (x.lof_flags == '')).map(
        lambda x: x.gene_symbol
    ))
    mt = mt.explode_rows(mt.gene_lofs)
    # indel=hl.is_indel(mt.alleles[0], mt.alleles[1])
    return mt.group_rows_by(mt.gene_lofs).aggregate(
        num_homs=hl.agg.count_where(mt.GT.is_hom_var()),
        num_hets=hl.agg.count_where(mt.GT.is_het()),
        defined_sites=hl.agg.count_where(hl.is_defined(mt.GT))
    )


def main(args):
    hl.init(log='/constraint.log')

    mt = get_gnomad_data('exomes', adj=True, release_samples=True)
    # mt = hl.filter_intervals(mt, [hl.parse_locus_interval('22')])
    annotations_mt = hl.read_matrix_table(processed_exomes_ht_path)
    exome_mt = mt.annotate_rows(**annotations_mt[mt.row_key, :])
    autosomes_exome_mt = filter_to_autosomes(exome_mt)

    # pass_autosomes_exome_mt = filter_to_pass(autosomes_exome_mt)
    pass_autosomes_exome_mt = autosomes_exome_mt.filter_rows(~autosomes_exome_mt.fail_filters)
    rare_autosomes_exome_mt = filter_by_frequency(pass_autosomes_exome_mt, 'below', 0.05)
    rare_autosomes_exome_mt = filter_by_frequency(rare_autosomes_exome_mt, 'above', 0)

    mt = generate_gene_lof_matrix(rare_autosomes_exome_mt)
    mt.write(gene_lof_matrix, args.overwrite)
    send_message(args.slack_channel, 'Gene LoF matrix computed!')

    mt = hl.read_matrix_table(gene_lof_matrix)
    # mt = mt.group_rows_by('gene_lofs').aggregate(defined_sites=hl.agg.sum(mt.defined_sites),
    #                                              num_homs=hl.agg.sum(mt.num_homs),
    #                                              num_hets=hl.agg.sum(mt.num_hets))
    ht = mt.annotate_rows(no_lofs=hl.agg.count_where((mt.defined_sites > 0) & (mt.num_homs + mt.num_hets == 0)),
                          obs_lof=hl.agg.count_where(mt.num_homs > 0),
                          defined=hl.agg.count_where(mt.defined_sites > 0)).rows()
    ht = ht.annotate(p=1 - hl.sqrt(ht.no_lofs / ht.defined))
    ht = ht.annotate(exp_lof=ht.defined * ht.p * ht.p)
    ht2 = ht.annotate(oe=ht.obs_lof / ht.exp_lof)
    ht2 = ht2.filter(ht2.exp_lof > 10)
    ht2.order_by(hl.asc(ht2.oe)).select('obs_lof', 'exp_lof', 'oe', 'p').show(20)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
