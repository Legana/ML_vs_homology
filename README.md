
# Machine learning AMP prediction models vs. homology (BLAST) to find AMPs in proteomes

-   [Find species represented with the most AMPs in a constructed AMP
    database and construct positive
    datasets](01_amp_training_data_preparations.md) : Rmd file
    [01_amp_training_data_preparations.Rmd](01_amp_training_data_preparations.Rmd)
-   [Test BLAST and classification
    methods](02_method_evaluation_on_proteomes.md) : Rmd file
    [02_method_evaluation_on_proteomes.Rmd](02_method_evaluation_on_proteomes.Rmd)
-   [Use phylogenetic data to assess performance over a taxonomic
    scale](03_taxonomic_distance_vs_performance.md) : Rmd file
    [03_taxonomic_distance_vs_performance.Rmd](03_taxonomic_distance_vs_performance.Rmd)

The files required to run the code in these Rmd files can be obtained by
clicking [here](https://cloudstor.aarnet.edu.au/plus/s/yXYa5zVk5rrvRpz)
or by using the command:

``` bash
wget 'https://cloudstor.aarnet.edu.au/plus/s/yXYa5zVk5rrvRpz/download' -O data.tgz
tar -zxvf data.tgz 
```

### `sessionInfo()`

    R version 4.1.2 (2021-11-01)
    Platform: x86_64-apple-darwin17.0 (64-bit)
    Running under: macOS Monterey 12.3.1

    Matrix products: default
    LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

    locale:
    [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] ggtree_3.0.4        randomcoloR_1.1.0.1 broom_0.7.11        ggtext_0.1.1        pals_1.7           
     [6] treeio_1.16.2       ape_5.6-1           precrec_0.12.7      patchwork_1.1.1     ampir_1.1.0        
    [11] forcats_0.5.1       stringr_1.4.0       dplyr_1.0.7         purrr_0.3.4         readr_2.1.1        
    [16] tidyr_1.1.4         tibble_3.1.6        ggplot2_3.3.5       tidyverse_1.3.1    

    loaded via a namespace (and not attached):
     [1] Rtsne_0.15           colorspace_2.0-2     ellipsis_0.3.2       class_7.3-19         Peptides_2.4.4      
     [6] fs_1.5.2             aplot_0.1.2          gridtext_0.1.4       dichromat_2.0-0      rstudioapi_0.13     
    [11] listenv_0.8.0        prodlim_2019.11.13   fansi_1.0.2          lubridate_1.8.0      xml2_1.3.3          
    [16] codetools_0.2-18     splines_4.1.2        knitr_1.37           jsonlite_1.7.2       pROC_1.18.0         
    [21] caret_6.0-90         cluster_2.1.2        dbplyr_2.1.1         mapproj_1.2.8        compiler_4.1.2      
    [26] httr_1.4.2           backports_1.4.1      assertthat_0.2.1     Matrix_1.3-4         fastmap_1.1.0       
    [31] lazyeval_0.2.2       cli_3.1.0            htmltools_0.5.2      tools_4.1.2          gtable_0.3.0        
    [36] glue_1.6.0           reshape2_1.4.4       maps_3.4.0           V8_4.0.0             Rcpp_1.0.8          
    [41] cellranger_1.1.0     vctrs_0.3.8          nlme_3.1-153         iterators_1.0.13     timeDate_3043.102   
    [46] gower_0.2.2          xfun_0.30            globals_0.14.0       rvest_1.0.2          lifecycle_1.0.1     
    [51] future_1.23.0        MASS_7.3-54          scales_1.1.1         ipred_0.9-12         hms_1.1.1           
    [56] parallel_4.1.2       yaml_2.2.1           curl_4.3.2           ggfun_0.0.4          yulab.utils_0.0.4   
    [61] rpart_4.1-15         stringi_1.7.6        foreach_1.5.1        tidytree_0.3.7       lava_1.6.10         
    [66] rlang_0.4.12         pkgconfig_2.0.3      evaluate_0.14        lattice_0.20-45      recipes_0.1.17      
    [71] tidyselect_1.1.1     parallelly_1.30.0    plyr_1.8.6           magrittr_2.0.1       R6_2.5.1            
    [76] generics_0.1.1       DBI_1.1.2            pillar_1.6.4         haven_2.4.3          withr_2.4.3         
    [81] survival_3.2-13      nnet_7.3-16          future.apply_1.8.1   modelr_0.1.8         crayon_1.4.2        
    [86] utf8_1.2.2           tzdb_0.2.0           rmarkdown_2.13       grid_4.1.2           readxl_1.3.1        
    [91] data.table_1.14.2    ModelMetrics_1.2.2.2 reprex_2.0.1         digest_0.6.29        gridGraphics_0.5-1  
    [96] stats4_4.1.2         munsell_0.5.0        ggplotify_0.1.0  
