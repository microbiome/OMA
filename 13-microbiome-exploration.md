# Microbiome Exploration {#microbiome-exploration}

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

This chapter focuses on the exploration of microbiome data and establish 
commonly used descriptors of a microbiome. The main difference to quality
control is that the exploration assumes the technical aspects of the dataset
have been investigated to your satisfaction. Generally speaking at this point 
you should be quite certain, that the dataset doesn't suffer from severe 
technical biases or you should at least be aware of potential problems.

In reality you might need to go back and forth between QC and exploration, 
since you discover through exploration of your dataset technical aspects you 
need to check.


```r
library(mia)
data("GlobalPatterns")
se <- GlobalPatterns 
```
## Prevalence

Prevalence is a measurements, which describes in how many samples certain
microbes were detected.

Investigating the prevalence of microbes allows you either to focus on changes,
which pertain to most of the samples, or to focus on less often found microbes,
which are nonetheless abundantly found in a small number of samples.

Population prevalence (frequency) at a 1% relative abundance threshold:


```r
head(getPrevalence(se, detection = 1/100, sort = TRUE, as_relative = TRUE))
```

```
## 331820 158660  98605 326977 145149 114821 
## 0.2308 0.2308 0.1923 0.1923 0.1538 0.1538
```

Population prevalence (frequency) at the absolute abundance threshold at read count 1:


```r
head(getPrevalence(se, detection = 1, sort = TRUE, abund_values = "counts",
                   as_relative = FALSE))
```

```
## 145149 114821 108747 526804  98605 180658 
##      1      1      1      1      1      1
```

### Prevalent microbiota analysis

If you only need the names of the prevalent taxa, do as follows. This
returns the taxa that exceed the given prevalence and detection
thresholds.


```r
prev <- getPrevalentTaxa(se, detection = 0, prevalence = 50/100)
```

See also related functions for the analysis of rare and variable taxa
(rareMembers; rareAbundance; lowAbundance). 

### Plotting prevalence

TODO


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04 LTS

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
 [1] mia_0.0.0.9009                   MicrobiomeExperiment_0.99.0.9014
 [3] Biostrings_2.58.0                XVector_0.30.0                  
 [5] TreeSummarizedExperiment_1.6.0   SingleCellExperiment_1.12.0     
 [7] SummarizedExperiment_1.20.0      Biobase_2.50.0                  
 [9] GenomicRanges_1.42.0             GenomeInfoDb_1.26.1             
[11] IRanges_2.24.0                   S4Vectors_0.28.0                
[13] BiocGenerics_0.36.0              MatrixGenerics_1.2.0            
[15] matrixStats_0.57.0               BiocStyle_2.18.1                
[17] rebook_1.0.0                     BiocManager_1.30.10             

loaded via a namespace (and not attached):
 [1] tidyselect_1.1.0          xfun_0.19                
 [3] beachmat_2.6.2            purrr_0.3.4              
 [5] lattice_0.20-41           vctrs_0.3.5              
 [7] generics_0.1.0            htmltools_0.5.0          
 [9] yaml_2.2.1                XML_3.99-0.5             
[11] rlang_0.4.9               pillar_1.4.7             
[13] scuttle_1.0.3             glue_1.4.2               
[15] BiocParallel_1.24.1       CodeDepends_0.6.5        
[17] GenomeInfoDbData_1.2.4    lifecycle_0.2.0          
[19] stringr_1.4.0             zlibbioc_1.36.0          
[21] codetools_0.2-18          evaluate_0.14            
[23] knitr_1.30                callr_3.5.1              
[25] ps_1.4.0                  Rcpp_1.0.5               
[27] DelayedArray_0.16.0       graph_1.68.0             
[29] digest_0.6.27             stringi_1.5.3            
[31] bookdown_0.21             processx_3.4.5           
[33] dplyr_1.0.2               grid_4.0.3               
[35] tools_4.0.3               bitops_1.0-6             
[37] magrittr_2.0.1            RCurl_1.98-1.2           
[39] tibble_3.0.4              tidyr_1.1.2              
[41] pkgconfig_2.0.3           crayon_1.3.4             
[43] ape_5.4-1                 ellipsis_0.3.1           
[45] Matrix_1.2-18             DelayedMatrixStats_1.12.1
[47] sparseMatrixStats_1.2.0   rmarkdown_2.5            
[49] R6_2.5.0                  nlme_3.1-150             
[51] compiler_4.0.3           
```
</div>
