from gnomad_hail import *
from gnomad_hail.utils.plotting import *
import statsmodels.formula.api as smf


def reverse_complement_bases(bases: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return hl.delimit(hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), '')
    # return bases[::-1].map(lambda x: flip_base(x))


def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
    return (hl.switch(base)
            .when('A', 'T')
            .when('T', 'A')
            .when('G', 'C')
            .when('C', 'G')
            .default(base))


def collapse_strand(ht: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    collapse_expr = {
        'ref': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.ref), ht.ref),
        'alt': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.alt), ht.alt),
        'context': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                           reverse_complement_bases(ht.context), ht.context),
        'was_flipped': (ht.ref == 'G') | (ht.ref == 'T')
    }
    return ht.annotate(**collapse_expr) if isinstance(ht, hl.Table) else ht.annotate_rows(**collapse_expr)


def downsampling_counts_expr(ht: Union[hl.Table, hl.MatrixTable], pop: str = 'global', variant_quality: str = 'adj',
                             singleton: bool = False) -> hl.expr.ArrayExpression:
    return hl.agg.array_sum(
        hl.map(lambda f: hl.int(f.AC[1] == 1) if singleton else hl.int(f.AC[1] > 0), hl.sorted(
            hl.filter(
                lambda f: (f.meta.size() == 3) & (f.meta.get('group') == variant_quality) &
                          (f.meta.get('pop') == pop) & f.meta.contains('downsampling'),
                ht.freq),
            key=lambda f: hl.int(f.meta['downsampling'])
        )))


def count_variants(ht: hl.Table,
                   count_singletons: bool = False, count_downsamplings: Optional[List[str]] = None,
                   additional_grouping: Optional[List[str]] = (), partition_hint: int = 1,
                   omit_methylation: bool = False, return_type_only: bool = False) -> Union[hl.Table, Any]:
    """
    Count variants by context, ref, alt, methylation_level
    """

    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    if not omit_methylation:
        grouping = grouping.annotate(methylation_level=ht.methylation_level)
    for group in additional_grouping:
        grouping = grouping.annotate(**{group: ht[group]})

    if count_singletons:
        singleton = hl.any(lambda f: (f.meta.size() == 1) & (f.meta.get('group') == 'adj') & (f.AC[1] == 1), ht.freq)

    if count_downsamplings:
        # Slower, but more flexible (allows for downsampling agg's)
        output = {'variant_count': hl.agg.count()}
        for pop in count_downsamplings:
            output[f'downsampling_counts_{pop}'] = downsampling_counts_expr(ht, pop)
        if count_singletons:
            output['singleton_count'] = hl.agg.count_where(singleton)
            for pop in count_downsamplings:
                output[f'singleton_downsampling_counts_{pop}'] = downsampling_counts_expr(ht, pop, singleton=True)
        return ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**output)
    else:
        agg = {'variant_count': hl.agg.counter(grouping)}
        if count_singletons:
            agg['singleton_count'] = hl.agg.counter(hl.agg.filter(singleton, grouping))

        if return_type_only:
            return agg['variant_count'].dtype
        else:
            return ht.aggregate(hl.struct(**agg))


def annotate_variant_types(t: Union[hl.MatrixTable, hl.Table],
                           heptamers: bool = False) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (((t.ref == "A") & (t.alt == "G")) | ((t.ref == "G") & (t.alt == "A")) |
                       ((t.ref == "T") & (t.alt == "C")) | ((t.ref == "C") & (t.alt == "T")))
    cpg_expr = (((t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1:mid_index] == 'C')) |
                ((t.ref == "C") & (t.alt == "T") & (t.context[mid_index + 1:mid_index + 2] == 'G')))
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (hl.case()
                         .when(t.cpg, 'CpG')
                         .when(t.transition, 'non-CpG transition')
                         .default('transversion'))
    variant_type_model_expr = hl.cond(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)
    else:
        return t.annotate(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)


def round_fields(ht: hl.Table, field: str, bin_size: int) -> hl.Table:
    """
    Round field in table to bin_size
    """
    return ht.annotate(**{field: hl.int(ht[field] / bin_size) * bin_size})


def set_kt_cols_to_zero(ht: hl.Table, fields: List[hl.expr.NumericExpression]) -> hl.Table:
    """
    Sets values of fields in ht to zero if missing
    """
    return ht.annotate(**{field: hl.or_else(field, 0) for field in fields})


def rebin_methylation(t: Union[hl.MatrixTable, hl.Table], bins: int=20) -> Union[hl.MatrixTable, hl.Table]:
    """
    Rebins methylation.level based on methylation.value to `bins` (assumes bi-allelic)
    """
    methylation_expr = t.methylation.annotate(level=hl.or_missing(
        hl.is_transition(hl.alleles[0], hl.alleles[1]),
        hl.range(bins - 1, -1, -1).find(lambda e: t.methylation.value * bins >= e)))
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(methylation=methylation_expr)
    else:
        return t.annotate(methylation=methylation_expr)


def trimer_from_heptamer(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    trimer_expr = hl.cond(hl.len(t.context) == 7, t.context[2:5], t.context)
    return t.annotate_rows(context=trimer_expr) if isinstance(t, hl.MatrixTable) else t.annotate(context=trimer_expr)


def filter_by_frequency(t: Union[hl.MatrixTable, hl.Table], frequency: float = None, allele_count: int = None,
                        population: str = None, subpop: str = None, downsampling: int = None,
                        direction: str = 'equal', keep: bool = True, adj: bool = True) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter MatrixTable or Table with gnomAD-format frequency data (assumed bi-allelic/split)
    (i.e. Array[Struct(Array[AC], Array[AF], AN, homozygote_count, meta)])
    At least one of frequency or allele_count is required.
    Subpop can be specified without a population if desired.

    :param MatrixTable t: Input MatrixTable or Table
    :param float frequency: Frequency to filter by
    :param int allele_count: Allele count to filter by
    :param str population:
    :param str subpop:
    :param int downsampling:
    :param str direction: One of "above", "below", and "equal" (how to apply the filter, default equal)
    :param bool keep: Whether to keep rows passing this frequency (passed to filter_rows)
    :param bool adj: Whether to use adj frequency
    :return: Filtered MatrixTable or Table
    :rtype: MatrixTable or Table
    """
    if frequency is None and allele_count is None:
        raise ValueError('At least one of frequency or allele_count must be specified')
    if direction not in ('above', 'below', 'equal'):
        raise ValueError('direction needs to be one of "above", "below", or "equal"')
    group = 'adj' if adj else 'raw'
    criteria = [lambda f: f.meta.get('group') == group]
    if frequency is not None:
        if direction == 'above':
            criteria.append(lambda f: f.AF[1] > frequency)
        elif direction == 'below':
            criteria.append(lambda f: f.AF[1] < frequency)
        else:
            criteria.append(lambda f: f.AF[1] == frequency)
    if allele_count is not None:
        if direction == 'above':
            criteria.append(lambda f: f.AC[1] > allele_count)
        elif direction == 'below':
            criteria.append(lambda f: f.AC[1] < allele_count)
        else:
            criteria.append(lambda f: f.AC[1] == allele_count)
    size = 1
    if population:
        criteria.append(lambda f: f.meta.get('pop') == population)
        size += 1
    if subpop:
        criteria.append(lambda f: f.meta.get('subpop') == subpop)
        size += 1
        # If one supplies a subpop but not a population, this will ensure this gets it right
        if not population: size += 1
    if downsampling:
        criteria.append(lambda f: f.meta.get('downsampling') == str(downsampling))
        size += 1
        if not population:
            size += 1
            criteria.append(lambda f: f.meta.get('pop') == 'global')
        if subpop:
            raise Exception('No downsampling data for subpopulations implemented')
    criteria.append(lambda f: f.meta.size() == size)

    def combine_functions(func_list, x):
        cond = func_list[0](x)
        for c in func_list[1:]:
            cond &= c(x)
        return cond

    filt = lambda x: combine_functions(criteria, x)
    criteria = hl.any(filt, t.freq)
    return t.filter_rows(criteria, keep=keep) if isinstance(t, hl.MatrixTable) else t.filter(criteria, keep=keep)


def filter_for_mu(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    """
    Filter to non-coding annotations, remove GERP outliers
    GERP cutoffs determined by finding 5% and 95% percentiles on list generated by:
    ```
    gerp_data = ht.aggregate(gerp=hl.agg.hist(context_ht.gerp, -12.3, 6.17, 100))
    cumulative_data = np.cumsum(summary_hist.gerp.bin_freq) + summary_hist.gerp.n_smaller
    np.append(cumulative_data, [cumulative_data[-1] + summary_hist.gerp.n_larger])
    list(zip(summary_hist.gerp.bin_edges, cumulative_data / max(cumulative_data)))
    ```
    """
    # criteria = ((t.gerp > -3.9885) & (t.gerp < 2.6607) &
    #             ((t.vep.most_severe_consequence == 'intron_variant') |
    #              (t.vep.most_severe_consequence == 'intergenic_variant')))
    criteria = ((t.gerp > -3.9885) & (t.gerp < 2.6607) &
                ((t.csq == 'intron_variant') |
                 (t.csq == 'intergenic_variant')))
    return t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)


def filter_to_pass(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    return t
    # criteria = hl.is_missing(t.filters)
    # return t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)


def filter_vep(t: Union[hl.MatrixTable, hl.Table],
               canonical: bool = True, synonymous: bool = True) -> Union[hl.MatrixTable, hl.Table]:
    if synonymous: t = filter_vep_to_synonymous_variants(t)
    if canonical: t = filter_vep_to_canonical_transcripts(t)
    criteria = hl.is_defined(t.vep.transcript_consequences)
    return t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)


def fast_filter_vep(t: Union[hl.Table, hl.MatrixTable], vep_root: str = 'vep') -> Union[hl.Table, hl.MatrixTable]:
    from gnomad_hail.utils.constants import CSQ_ORDER

    csqs = hl.literal(CSQ_ORDER)

    def add_most_severe_consequence(tc: hl.expr.StructExpression) -> hl.expr.StructExpression:
        """
        Add most_severe_consequence annotation to transcript consequences
        This is for a given transcript, as there are often multiple annotations for a single transcript:
        e.g. splice_region_variant&intron_variant -> splice_region_variant
        """
        csq_term_array = tc.consequence_terms.split('&')
        return tc.annotate(
            most_severe_consequence=csqs.find(lambda c: csq_term_array.contains(c))
        )

    transcript_csqs = t[vep_root].transcript_consequences.map(add_most_severe_consequence)
    transcript_csqs = transcript_csqs.filter(lambda csq: (csq.most_severe_consequence == "synonymous_variant") &
                                                         (csq.canonical == 1))
    vep_data = t[vep_root].annotate(transcript_consequences=transcript_csqs)
    t = t.annotate_rows(**{vep_root: vep_data}) if isinstance(t, hl.MatrixTable) else t.annotate(**{vep_root: vep_data})
    criteria = hl.is_defined(t.vep.transcript_consequences)
    return t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)


def remove_coverage_outliers(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
    """
    Keep only loci where genome coverage was between 15 and 60
    """
    criteria = (t.coverage.genomes.mean >= 15) & (t.coverage.genomes.mean <= 60)
    return t.filter_rows(criteria) if isinstance(t, hl.MatrixTable) else t.filter(criteria)


# Misc
def maps(ht: hl.Table, mutation_ht: hl.Table, vep_root: str = 'vep') -> hl.Table:
    ht = ht.annotate(csq=ht[vep_root]['most_severe_consequence'])
    ht = count_variants(ht, count_singletons=True, additional_grouping=('csq', ))
    ht = ht.annotate(mu=mutation_ht[ht.key].mu_snp,
                     ps=ht.singleton_count / ht.variant_count)
    syn_ps_pd = ht.filter(ht.csq == 'synonymous_variant').to_pandas()

    lm = smf.ols(formula='ps ~ mu', data=syn_ps_pd).fit()
    slope = lm.params['mu']
    intercept = lm.params['Intercept']
    ht = ht.annotate(expected_singletons=(ht.mu * slope + intercept) * ht.variant_count)

    agg_ht = (ht.group_by('csq')
              .aggregate(singleton_count=hl.agg.sum(ht.singleton_count),
                         expected_singletons=hl.agg.sum(ht.expected_singletons),
                         variant_count=hl.agg.sum(ht.variant_count)))
    agg_ht = agg_ht.annotate(ps=agg_ht.singleton_count / agg_ht.variant_count,
                             maps=(agg_ht.singleton_count - agg_ht.expected_singletons) / agg_ht.variant_count)
    agg_ht = agg_ht.annotate(sem_ps=(agg_ht.ps * (1 - agg_ht.ps) / agg_ht.variant_count) ** 0.5)
    return agg_ht
