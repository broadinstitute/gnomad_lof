#!/usr/bin/env python3

__author__ = 'konradk'

from constraint_utils import *

gene_lof_matrix = f'{root}/recessive/gene_matrix/full_gene_lof_matrix.mt'


def generate_gene_lof_matrix(mt: hl.MatrixTable) -> hl.MatrixTable:
    mt = process_consequences(mt)
    mt = mt.select_rows(gene_lofs=mt.vep.worst_csq_by_gene.filter(
        lambda x: (x.lof == 'HC') & (x.lof_flags == '')).map(
        lambda x: x.gene_symbol
    ))
    mt = mt.explode_rows(mt.gene_lofs)
    return mt.group_rows_by(mt.gene_lofs, indel=hl.is_indel(mt.alleles[0], mt.alleles[1])).aggregate(
        num_homs=hl.agg.count_where(mt.GT.is_hom_var()),
        num_hets=hl.agg.count_where(mt.GT.is_het()),
        defined_sites=hl.agg.count_where(hl.is_defined(mt.GT))
    )


def main(args):
    hl.init(log='/constraint.log')

    full_exome_ht = prepare_ht(hl.read_table(processed_exomes_ht_path), args.trimers)
    autosomes_exome_ht = filter_to_autosomes(full_exome_ht)

    # pass_autosomes_exome_ht = filter_to_pass(autosomes_exome_ht)  # TODO
    rare_autosomes_exome_ht = filter_by_frequency(autosomes_exome_ht, 0.001, direction='below')

    mt = generate_gene_lof_matrix(rare_autosomes_exome_ht)
    mt.write(gene_lof_matrix)
    send_message(args.slack_channel, 'Gene LoF matrix computed!')

    mt = hl.read_matrix_table(gene_lof_matrix)
    mt.annotate_rows(no_lofs=hl.agg.count_where(mt.num_homs + mt.num_hets == 0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite everything', action='store_true')
    parser.add_argument('--slack_channel', help='Send message to Slack channel/user', default='@konradjk')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
