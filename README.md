# Overview

Public GitHub repository for downloading Swedish freshwater aquatic plant and water chemistry data, combining them into a single dataset, and performing analysis. Companion script to manuscript under review.

Main files are
* `1-s2.0-S0304377019300300-mmc2.xlsx`: List from [Murphy et al., 2019](https://doi.org/10.1016/j.aquabot.2019.06.006) listing main vascular plants considered as "aquatic". Downloaded from [supplementary material](https://ars.els-cdn.com/content/image/1-s2.0-S0304377019300300-mmc2.xlsx).
* `seChemAPI.R`: Raw companion chemistry download and cleaning into long format.
* `seDataCombine.R`: Combining chemistry and plant data together into a single wide dataset.
    * Run after `seChemAPI.R` and `seMacroAPI.R`.
* `seMacroAPI.R`: Raw plant data download and cleaning into wide format.
* `seSpatialModels.R`: Running several [`spCF`](https://cran.r-project.org/package=spCF) models for spatial analysis on the data.
    * Run after `seDataCombine.R`.

### Session Information
Code was most recently run in the below environment.
```
R version 4.5.2 (2025-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] here_1.0.2         readxl_1.4.5       jsonlite_2.0.0    
 [4] httr_1.4.7         ggpattern_1.3.1    sf_1.0-22         
 [7] viridis_0.6.5      viridisLite_0.4.2  ggcorrplot_0.1.4.1
[10] spCF_0.1.1         lubridate_1.9.4    forcats_1.0.1     
[13] stringr_1.6.0      dplyr_1.2.1        purrr_1.2.0       
[16] readr_2.1.6        tidyr_1.3.1        tibble_3.3.0      
[19] ggplot2_4.0.3      tidyverse_2.0.0    devtools_2.4.6    
[22] usethis_3.2.1     

loaded via a namespace (and not attached):
 [1] dotCall64_1.2      gtable_0.3.6       spam_2.11-1       
 [4] remotes_2.5.0      lattice_0.22-7     tzdb_0.5.0        
 [7] vctrs_0.7.3        tools_4.5.2        generics_0.1.4    
[10] proxy_0.4-27       pkgconfig_2.0.3    Matrix_1.7-4      
[13] KernSmooth_2.23-26 RColorBrewer_1.1-3 S7_0.2.1          
[16] lifecycle_1.0.5    compiler_4.5.2     farver_2.1.2      
[19] FNN_1.1.4.1        fields_17.1        maps_3.4.3        
[22] class_7.3-23       pillar_1.11.1      nloptr_2.2.1      
[25] ellipsis_0.3.2     classInt_0.4-11    cachem_1.1.0      
[28] dbscan_1.2.4       sessioninfo_1.2.3  tidyselect_1.2.1  
[31] stringi_1.8.7      rprojroot_2.1.1    fastmap_1.2.0     
[34] grid_4.5.2         cli_3.6.5          magrittr_2.0.4    
[37] pkgbuild_1.4.8     e1071_1.7-16       withr_3.0.2       
[40] scales_1.4.0       timechange_0.3.0   gridExtra_2.3     
[43] cellranger_1.1.0   ranger_0.17.0      hms_1.1.4         
[46] memoise_2.0.1      rlang_1.2.0        Rcpp_1.1.0        
[49] glue_1.8.0         DBI_1.2.3          pkgload_1.4.1     
[52] rstudioapi_0.17.1  R6_2.6.1           fs_1.6.6          
[55] units_1.0-0
```