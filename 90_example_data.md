# (PART) Appendix {-}

# Example Data {#example-data}

<script>
document.addEventListener("click", function (event) {
    if (event.target.classList.contains("rebook-collapse")) {
        event.target.classList.toggle("active");
        var content = event.target.nextElementSibling;
        if (content.style.display === "block") {
            content.style.display = "none";
        } else {
            content.style.display = "block";
        }
    }
})
</script>

<style>
.rebook-collapse {
  background-color: #eee;
  color: #444;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

.rebook-content {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: #f1f1f1;
}
</style>

## Package data {#package-data}

The datasets in `mia` are conversions of the `phyloseq` datasets 
`GlobalPatterns` `enterotype`, `esophagus` and `soilrep`.

### GlobalPatterns


```r
library(mia)
# Example how to load data
data("GlobalPatterns", package="mia")

GlobalPatterns
```

```
## class: TreeSummarizedExperiment 
## dim: 19216 26 
## metadata(0):
## assays(1): counts
## rownames(19216): 549322 522457 ... 200359 271582
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```




## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=C             
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] mia_1.1.11                     TreeSummarizedExperiment_2.1.4
 [3] Biostrings_2.61.2              XVector_0.33.0                
 [5] SingleCellExperiment_1.15.1    SummarizedExperiment_1.23.1   
 [7] Biobase_2.53.0                 GenomicRanges_1.45.0          
 [9] GenomeInfoDb_1.29.3            IRanges_2.27.2                
[11] S4Vectors_0.31.1               BiocGenerics_0.39.2           
[13] MatrixGenerics_1.5.3           matrixStats_0.60.0            
[15] BiocStyle_2.21.3               rebook_1.3.0                  

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0            colorspace_2.0-2           
 [3] ellipsis_0.3.2              scuttle_1.3.1              
 [5] BiocNeighbors_1.11.0        ggrepel_0.9.1              
 [7] bit64_4.0.5                 fansi_0.5.0                
 [9] decontam_1.13.0             splines_4.1.0              
[11] codetools_0.2-18            sparseMatrixStats_1.5.2    
[13] cachem_1.0.5                knitr_1.33                 
[15] scater_1.21.3               jsonlite_1.7.2             
[17] cluster_2.1.2               graph_1.71.2               
[19] BiocManager_1.30.16         compiler_4.1.0             
[21] assertthat_0.2.1            Matrix_1.3-4               
[23] fastmap_1.1.0               lazyeval_0.2.2             
[25] BiocSingular_1.9.1          htmltools_0.5.1.1          
[27] tools_4.1.0                 rsvd_1.0.5                 
[29] gtable_0.3.0                glue_1.4.2                 
[31] GenomeInfoDbData_1.2.6      reshape2_1.4.4             
[33] dplyr_1.0.7                 Rcpp_1.0.7                 
[35] jquerylib_0.1.4             vctrs_0.3.8                
[37] ape_5.5                     nlme_3.1-152               
[39] DECIPHER_2.21.0             DelayedMatrixStats_1.15.2  
[41] xfun_0.25                   stringr_1.4.0              
[43] beachmat_2.9.1              lifecycle_1.0.0            
[45] irlba_2.3.3                 XML_3.99-0.6               
[47] zlibbioc_1.39.0             MASS_7.3-54                
[49] scales_1.1.1                parallel_4.1.0             
[51] yaml_2.2.1                  memoise_2.0.0              
[53] gridExtra_2.3               ggplot2_3.3.5              
[55] sass_0.4.0                  stringi_1.7.3              
[57] RSQLite_2.2.7               ScaledMatrix_1.1.0         
[59] tidytree_0.3.4              permute_0.9-5              
[61] filelock_1.0.2              BiocParallel_1.27.4        
[63] rlang_0.4.11                pkgconfig_2.0.3            
[65] bitops_1.0-7                evaluate_0.14              
[67] lattice_0.20-44             purrr_0.3.4                
[69] treeio_1.17.2               CodeDepends_0.6.5          
[71] bit_4.0.4                   tidyselect_1.1.1           
[73] plyr_1.8.6                  magrittr_2.0.1             
[75] bookdown_0.22               R6_2.5.0                   
[77] generics_0.1.0              DelayedArray_0.19.1        
[79] DBI_1.1.1                   mgcv_1.8-36                
[81] pillar_1.6.2                RCurl_1.98-1.3             
[83] tibble_3.1.3                dir.expiry_1.1.0           
[85] crayon_1.4.1                utf8_1.2.2                 
[87] rmarkdown_2.10              viridis_0.6.1              
[89] grid_4.1.0                  blob_1.2.2                 
[91] vegan_2.5-7                 digest_0.6.27              
[93] tidyr_1.1.3                 munsell_0.5.0              
[95] DirichletMultinomial_1.35.0 beeswarm_0.4.0             
[97] viridisLite_0.4.0           vipor_0.4.5                
[99] bslib_0.2.5.1              
```
</div>
