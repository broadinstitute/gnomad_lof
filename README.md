# Code for gnomAD LoF flagship manuscript

Manuscript at: [Karczewski et al., 2019](https://www.biorxiv.org/content/10.1101/531210v2)

## Details

This repo serves as a home for two main purposes:
- The constraint computation pipeline, written in [Hail 0.2](https://hail.is).
- The figure-generating code for the manuscript

### Constraint

The constraint pipeline, as initially described [here](https://www.ncbi.nlm.nih.gov/pubmed/25086666) and [here](https://www.nature.com/articles/nature19057), has been updated with a number of improvements as described [in the supplement of Karczewski et al](https://www.biorxiv.org/content/10.1101/531210v2.supplementary-material). Notably, the pipeline is now written in [Hail](https://hail.is), which enables scalability to large datasets like gnomAD, and can compute constraint against arbitrary sets of variants.

The main components of the pipeline can be found in `constraint/constraint.py`, which uses the public gnomAD data and a dataset of every possible variant (~9B variants) to compute the observed and expected number of variants per transcript/gene. This script is provided primarily for reference (it could be run with modifications, but cannot be run as-is, as it has paths to buckets on Google cloud hard-coded). Additionally, we combine the constraint data with aggregate LoF frequencies in `constraint/gene_lof_matrix.py` - this script cannot be run outside of the gnomAD team, as it requires access to the individual level data, but it is also provided for reference.

### Figure-generating code

The code to generate all the figures can be found in `R/`. These scripts use only aggregated data, and thus, can be run by anyone. They look for data files in `./data`, and if they are not found, downloads them as needed from the [public data repository](https://storage.googleapis.com/gnomad-public/papers/2019-flagship-lof/v1.0/).

Each figure panel can be generated individually, or figures as a whole. For instance, in `R/fig3_spectrum.R`, we provide a function called `figure3` which can generate the entirety of figure 3. Alternatively, running the code inside the function can generate each figure panel separately. Note that for some figures, on some R setups, attempting to generate the full figure by calling the function directly can crash R: we are uncertain of the cause of the issue, but it can be resolved by running the code inside the function step-wise.
