### Human genes description:
- list_CEGv2.tsv: list of 684 genes deemed essential in multiple cultured cell lines based on CRISPR/Cas screen data [Hart 2017](http://www.g3journal.org/content/7/8/2719)
- list_NEGv1.tsv: list of 927 genes deemed non-essential in multiple cultured cell lines based on CRISPR/Cas screen data [Hart 2017](http://www.g3journal.org/content/7/8/2719)

### Mouse genes description:

- list_mouse_het_lethal_genes.tsv: list of 404 mouse genes that has at least one report of heterozygous knockout being lethal. Extracted using [Mousemine](http://www.mousemine.org/)
- list_mouse_lethal_genes.tsv: list of 4071 mouse genes that has at least one report of knockout being lethal (any strain / any zygousity). Extracted using [MGI](http://www.informatics.jax.org)
- list_mouse_viable_genes.tsv: list of 282 mouse genes that has no any phenotype report other than "no phenotype associated". Extracted using [Mousemine](http://www.mousemine.org/)


For a given human gene, query against a list is defined to be true if at least one mouse homolog of that human gene are in the list. 
(e.g. if 1 of 10 homolog of human gene A is in the het_lethal list, then that human gene is defined as "mouse het lethal".)
The homolog correspondence is defined based on 
HMD_HumanPhenotype.rpt, which is extracted from [MGI](http://www.informatics.jax.org/downloads/reports/).

- HGene_MPhenotype_all.tsv: is the raw data for extracting viable genes, which is the result of querying "MP:0010769" in [Mousemine](http://www.mousemine.org/)
- MPhenotype_MGenotype_heterozygotes (1).tsv: is the raw data for mouse het-lethal genes, which is the result of querying "MP:0010769" in [Mousemine](http://www.mousemine.org/)
