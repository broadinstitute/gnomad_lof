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

The code was run using R 3.5.1 using the following packages:

```
> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin18.2.0 (64-bit)
Running under: macOS  10.14.3

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /opt/local/Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] hexbin_1.27.2        bindrcpp_0.2.2       cowplot_0.9.4        RMySQL_0.10.16       DBI_1.0.0           
 [6] ggrepel_0.8.0        pbapply_1.3-4        rlang_0.3.1          tidygraph_1.1.1      STRINGdb_1.22.0     
[11] meta_4.9-4           ggrastr_0.1.7        ggpubr_0.2           ggridges_0.5.1       readxl_1.2.0        
[16] corrr_0.3.0          corrplot_0.84        patchwork_0.0.1      naniar_0.4.1         plotROC_2.2.1       
[21] gghighlight_0.1.0    skimr_1.0.4          gapminder_0.3.0      trelliscopejs_0.1.18 scales_1.0.0        
[26] magrittr_1.5         slackr_1.4.2         plotly_4.8.0         broom_0.5.1          forcats_0.3.0       
[31] stringr_1.3.1        dplyr_0.7.8          purrr_0.2.5          readr_1.3.1          tidyr_0.8.2         
[36] tibble_2.0.1         tidyverse_1.2.1      Hmisc_4.1-1          ggplot2_3.1.0        Formula_1.2-3       
[41] survival_2.43-3      lattice_0.20-38     

loaded via a namespace (and not attached):
 [1] colorspace_1.4-0        visdat_0.5.2            mclust_5.4.2            htmlTable_1.13.1       
 [5] base64enc_0.1-3         rstudioapi_0.9.0        hash_2.2.6              bit64_0.9-7            
 [9] fansi_0.4.0             lubridate_1.7.4         sqldf_0.4-11            xml2_1.2.0             
[13] codetools_0.2-16        splines_3.5.1           knitr_1.21              jsonlite_1.6           
[17] Cairo_1.5-9             cluster_2.0.7-1         png_0.1-7               compiler_3.5.1         
[21] httr_1.4.0              backports_1.1.3         assertthat_0.2.0        Matrix_1.2-15          
[25] lazyeval_0.2.1          cli_1.0.1               acepack_1.4.1           htmltools_0.3.6        
[29] prettyunits_1.0.2       tools_3.5.1             igraph_1.2.2            gtable_0.2.0           
[33] glue_1.3.0              Rcpp_1.0.0              cellranger_1.1.0        gdata_2.18.0           
[37] nlme_3.1-137            autocogs_0.1.1          xfun_0.4                proto_1.0.0            
[41] rvest_0.3.2             gtools_3.8.1            DistributionUtils_0.6-0 hms_0.4.2              
[45] parallel_3.5.1          RColorBrewer_1.1-2      yaml_2.2.0              memoise_1.1.0          
[49] gridExtra_2.3           rpart_4.1-13            latticeExtra_0.6-28     stringi_1.2.4          
[53] RSQLite_2.1.1           plotrix_3.7-4           checkmate_1.9.1         caTools_1.17.1.1       
[57] chron_2.3-53            pkgconfig_2.0.2         bitops_1.0-6            bindr_0.1.1            
[61] labeling_0.3            htmlwidgets_1.3         bit_1.1-14              tidyselect_0.2.5       
[65] plyr_1.8.4              R6_2.3.0                gplots_3.0.1            generics_0.0.2         
[69] gsubfn_0.7              pillar_1.3.1            haven_2.0.0             foreign_0.8-71         
[73] withr_2.1.2             RCurl_1.95-4.11         nnet_7.3-12             modelr_0.1.2           
[77] crayon_1.3.4            utf8_1.1.4              KernSmooth_2.23-15      progress_1.2.0         
[81] data.table_1.12.0       blob_1.1.1              digest_0.6.18           webshot_0.5.1          
[85] munsell_0.5.0           viridisLite_0.3.0       egg_0.4.2              
```
