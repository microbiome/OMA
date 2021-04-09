# (PART) Appendix {-}

# Example data {#example-data}

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
R Under development (unstable) (2021-04-05 r80145)
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
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] mia_0.99.10                      TreeSummarizedExperiment_1.99.11
 [3] Biostrings_2.59.2                XVector_0.31.1                  
 [5] SingleCellExperiment_1.13.14     SummarizedExperiment_1.21.3     
 [7] Biobase_2.51.0                   GenomicRanges_1.43.4            
 [9] GenomeInfoDb_1.27.10             IRanges_2.25.7                  
[11] S4Vectors_0.29.15                BiocGenerics_0.37.1             
[13] MatrixGenerics_1.3.1             matrixStats_0.58.0              
[15] BiocStyle_2.19.2                 rebook_1.1.19                   
[17] BiocManager_1.30.12             

loaded via a namespace (and not attached):
 [1] nlme_3.1-152                bitops_1.0-6               
 [3] DirichletMultinomial_1.33.2 bit64_4.0.5                
 [5] filelock_1.0.2              tools_4.1.0                
 [7] bslib_0.2.4                 vegan_2.5-7                
 [9] utf8_1.2.1                  R6_2.5.0                   
[11] irlba_2.3.3                 vipor_0.4.5                
[13] mgcv_1.8-34                 DBI_1.1.1                  
[15] colorspace_2.0-0            permute_0.9-5              
[17] gridExtra_2.3               tidyselect_1.1.0           
[19] bit_4.0.4                   compiler_4.1.0             
[21] graph_1.69.0                BiocNeighbors_1.9.4        
[23] DelayedArray_0.17.10        bookdown_0.21              
[25] sass_0.3.1                  scales_1.1.1               
[27] stringr_1.4.0               digest_0.6.27              
[29] rmarkdown_2.7               scater_1.19.11             
[31] pkgconfig_2.0.3             htmltools_0.5.1.1          
[33] sparseMatrixStats_1.3.7     fastmap_1.1.0              
[35] rlang_0.4.10                RSQLite_2.2.5              
[37] DelayedMatrixStats_1.13.5   jquerylib_0.1.3            
[39] generics_0.1.0              jsonlite_1.7.2             
[41] BiocParallel_1.25.5         dplyr_1.0.5                
[43] RCurl_1.98-1.3              magrittr_2.0.1             
[45] BiocSingular_1.7.2          GenomeInfoDbData_1.2.4     
[47] scuttle_1.1.18              Matrix_1.3-2               
[49] Rcpp_1.0.6                  ggbeeswarm_0.6.0           
[51] munsell_0.5.0               fansi_0.4.2                
[53] DECIPHER_2.19.2             viridis_0.5.1              
[55] ape_5.4-1                   lifecycle_1.0.0            
[57] stringi_1.5.3               yaml_2.2.1                 
[59] MASS_7.3-53.1               debugme_1.1.0              
[61] zlibbioc_1.37.0             blob_1.2.1                 
[63] grid_4.1.0                  crayon_1.4.1               
[65] dir.expiry_0.99.4           lattice_0.20-41            
[67] splines_4.1.0               beachmat_2.7.7             
[69] CodeDepends_0.6.5           knitr_1.31                 
[71] pillar_1.5.1                codetools_0.2-18           
[73] ScaledMatrix_0.99.2         XML_3.99-0.6               
[75] glue_1.4.2                  evaluate_0.14              
[77] vctrs_0.3.7                 tidyr_1.1.3                
[79] gtable_0.3.0                purrr_0.3.4                
[81] assertthat_0.2.1            cachem_1.0.4               
[83] ggplot2_3.3.3               xfun_0.22                  
[85] rsvd_1.0.3                  viridisLite_0.3.0          
[87] tibble_3.1.0                memoise_2.0.0              
[89] beeswarm_0.3.1              cluster_2.1.1              
[91] ellipsis_0.3.1             
```
</div>
