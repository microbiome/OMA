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
data("GlobalPatterns")

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

### Enterotype


### Esophagus


### Soilrep


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
 [1] mia_1.1.7                      TreeSummarizedExperiment_2.1.3
 [3] Biostrings_2.61.1              XVector_0.33.0                
 [5] SingleCellExperiment_1.15.1    SummarizedExperiment_1.23.1   
 [7] Biobase_2.53.0                 GenomicRanges_1.45.0          
 [9] GenomeInfoDb_1.29.3            IRanges_2.27.0                
[11] S4Vectors_0.31.0               BiocGenerics_0.39.1           
[13] MatrixGenerics_1.5.1           matrixStats_0.60.0            
[15] BiocStyle_2.21.3               rebook_1.3.0                  

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0            colorspace_2.0-2           
 [3] ellipsis_0.3.2              scuttle_1.3.0              
 [5] BiocNeighbors_1.11.0        bit64_4.0.5                
 [7] fansi_0.5.0                 decontam_1.13.0            
 [9] splines_4.1.0               codetools_0.2-18           
[11] sparseMatrixStats_1.5.0     cachem_1.0.5               
[13] knitr_1.33                  scater_1.21.2              
[15] jsonlite_1.7.2              cluster_2.1.2              
[17] graph_1.71.2                BiocManager_1.30.16        
[19] compiler_4.1.0              assertthat_0.2.1           
[21] Matrix_1.3-4                fastmap_1.1.0              
[23] lazyeval_0.2.2              BiocSingular_1.9.1         
[25] htmltools_0.5.1.1           tools_4.1.0                
[27] rsvd_1.0.5                  gtable_0.3.0               
[29] glue_1.4.2                  GenomeInfoDbData_1.2.6     
[31] reshape2_1.4.4              dplyr_1.0.7                
[33] Rcpp_1.0.7                  jquerylib_0.1.4            
[35] vctrs_0.3.8                 ape_5.5                    
[37] nlme_3.1-152                DECIPHER_2.21.0            
[39] DelayedMatrixStats_1.15.0   xfun_0.24                  
[41] stringr_1.4.0               beachmat_2.9.0             
[43] lifecycle_1.0.0             irlba_2.3.3                
[45] XML_3.99-0.6                zlibbioc_1.39.0            
[47] MASS_7.3-54                 scales_1.1.1               
[49] parallel_4.1.0              yaml_2.2.1                 
[51] memoise_2.0.0               gridExtra_2.3              
[53] ggplot2_3.3.5               sass_0.4.0                 
[55] stringi_1.7.3               RSQLite_2.2.7              
[57] ScaledMatrix_1.1.0          tidytree_0.3.4             
[59] permute_0.9-5               filelock_1.0.2             
[61] BiocParallel_1.27.2         rlang_0.4.11               
[63] pkgconfig_2.0.3             bitops_1.0-7               
[65] evaluate_0.14               lattice_0.20-44            
[67] purrr_0.3.4                 treeio_1.17.2              
[69] CodeDepends_0.6.5           bit_4.0.4                  
[71] tidyselect_1.1.1            plyr_1.8.6                 
[73] magrittr_2.0.1              bookdown_0.22              
[75] R6_2.5.0                    generics_0.1.0             
[77] DelayedArray_0.19.1         DBI_1.1.1                  
[79] mgcv_1.8-36                 pillar_1.6.1               
[81] RCurl_1.98-1.3              tibble_3.1.3               
[83] dir.expiry_1.1.0            crayon_1.4.1               
[85] utf8_1.2.2                  rmarkdown_2.9              
[87] viridis_0.6.1               grid_4.1.0                 
[89] blob_1.2.2                  vegan_2.5-7                
[91] digest_0.6.27               tidyr_1.1.3                
[93] munsell_0.5.0               DirichletMultinomial_1.35.0
[95] beeswarm_0.4.0              viridisLite_0.4.0          
[97] vipor_0.4.5                 bslib_0.2.5.1              
```
</div>
