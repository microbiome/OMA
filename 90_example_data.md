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


Open demonstration data for testing and benchmarking purposes is
available from multiple locations. This chapter introduces some
options. The other chapters of this book provide ample examples about
the use of the data.


## Package data {#package-data}

The `mia` R package contains example data sets that are direct
conversions from the alternative `phyloseq` container to the
`TreeSummarizedExperiment` container.

List the [available datasets](https://microbiome.github.io/mia/reference/index.html) in the `mia` package:



```r
library(mia)
data(package="mia")
```

Check the documentation for one of these data sets:


```
## No documentation for 'GlobalPatterns' in specified packages and libraries:
## you could try '??GlobalPatterns'
```

Load the `GlobalPatterns` data from the `mia` package:


```r
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


## ExperimentHub data

The [ExperimentHub](https://bioconductor.org/packages/release/bioc/vignettes/ExperimentHub/inst/doc/ExperimentHub.html) provides a variety of data resources, including the [microbiomeDataSets](https://bioconductor.org/packages/devel/data/experiment/html/microbiomeDataSets.html) package.

A table of the available data sets is available through the `availableDataSets`
function.


```r
library(microbiomeDataSets)
availableDataSets()
```

```
##             Dataset
## 1  GrieneisenTSData
## 2    HintikkaXOData
## 3       LahtiMLData
## 4        LahtiMData
## 5       LahtiWAData
## 6      OKeefeDSData
## 7 SilvermanAGutData
## 8        SongQAData
## 9   SprockettTHData
```

All data are downloaded from ExperimentHub and cached for local
re-use. Check the [man pages of each
function](https://microbiome.github.io/microbiomeDataSets/reference/index.html)
for a detailed documentation of the data contents and references. The
following example data set is a *[MultiAssayExperiment](https://bioconductor.org/packages/3.14/MultiAssayExperiment)*:


```r
#mae <- HintikkaXOData()
```

Example data sets are available in *[SummarizedExperiment](https://bioconductor.org/packages/3.14/SummarizedExperiment)*, *[TreeSummarizedExperiment](https://bioconductor.org/packages/3.14/TreeSummarizedExperiment)*, and *[MultiAssayExperiment](https://bioconductor.org/packages/3.14/MultiAssayExperiment)* data containers; see the separate page on [alternative containers](https://microbiome.github.io/OMA/multitable.html) for more details.



## Other data sources

The [curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData) is an independent source that provides various example data sets as `(Tree)SummarizedExperiment` objects. This resource provides curated human microbiome data including gene families, marker abundance, marker presence, pathway abundance, pathway coverage, and relative abundance for samples from different body sites. See the package homepage for more details on data availability and access.

As one example, let us retrieve the AsnicarF data set:



```r
library(curatedMetagenomicData)
tse <- curatedMetagenomicData("AsnicarF_20.+.relative_abundance", dryrun = FALSE, counts = TRUE)
```

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

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
 [1] curatedMetagenomicData_3.1.4   microbiomeDataSets_1.1.5      
 [3] MultiAssayExperiment_1.19.11   TreeSummarizedExperiment_2.1.4
 [5] Biostrings_2.61.2              XVector_0.33.0                
 [7] SingleCellExperiment_1.15.2    SummarizedExperiment_1.23.4   
 [9] Biobase_2.53.0                 GenomicRanges_1.45.0          
[11] GenomeInfoDb_1.29.8            IRanges_2.27.2                
[13] S4Vectors_0.31.3               BiocGenerics_0.39.2           
[15] MatrixGenerics_1.5.4           matrixStats_0.60.1-9001       
[17] BiocStyle_2.21.3               rebook_1.3.1                  

loaded via a namespace (and not attached):
 [1] nlme_3.1-153                  bitops_1.0-7                 
 [3] bit64_4.0.5                   httr_1.4.2                   
 [5] filelock_1.0.2                tools_4.1.1                  
 [7] bslib_0.3.0                   utf8_1.2.2                   
 [9] R6_2.5.1                      DBI_1.1.1                    
[11] lazyeval_0.2.2                withr_2.4.2                  
[13] tidyselect_1.1.1              curl_4.3.2                   
[15] bit_4.0.4                     compiler_4.1.1               
[17] graph_1.71.2                  DelayedArray_0.19.2          
[19] bookdown_0.24                 sass_0.4.0                   
[21] rappdirs_0.3.3                stringr_1.4.0                
[23] digest_0.6.27                 rmarkdown_2.10               
[25] pkgconfig_2.0.3               htmltools_0.5.2              
[27] dbplyr_2.1.1                  fastmap_1.1.0                
[29] rlang_0.4.11                  RSQLite_2.2.8                
[31] shiny_1.6.0                   jquerylib_0.1.4              
[33] generics_0.1.0                jsonlite_1.7.2               
[35] BiocParallel_1.27.6           dplyr_1.0.7                  
[37] RCurl_1.98-1.4                magrittr_2.0.1               
[39] GenomeInfoDbData_1.2.6        Matrix_1.3-4                 
[41] Rcpp_1.0.7                    fansi_0.5.0                  
[43] ape_5.5                       lifecycle_1.0.0              
[45] stringi_1.7.4                 yaml_2.2.1                   
[47] zlibbioc_1.39.0               BiocFileCache_2.1.1          
[49] AnnotationHub_3.1.5           grid_4.1.1                   
[51] blob_1.2.2                    promises_1.2.0.1             
[53] parallel_4.1.1                ExperimentHub_2.1.4          
[55] crayon_1.4.1                  dir.expiry_1.1.0             
[57] lattice_0.20-44               KEGGREST_1.33.0              
[59] CodeDepends_0.6.5             knitr_1.33                   
[61] pillar_1.6.2                  codetools_0.2-18             
[63] XML_3.99-0.7                  glue_1.4.2                   
[65] BiocVersion_3.14.0            evaluate_0.14                
[67] BiocManager_1.30.16           png_0.1-7                    
[69] httpuv_1.6.2                  vctrs_0.3.8                  
[71] treeio_1.17.2                 purrr_0.3.4                  
[73] tidyr_1.1.3                   assertthat_0.2.1             
[75] cachem_1.0.6                  xfun_0.25                    
[77] mime_0.11                     xtable_1.8-4                 
[79] tidytree_0.3.5                later_1.3.0                  
[81] tibble_3.1.4                  AnnotationDbi_1.55.1         
[83] memoise_2.0.0                 interactiveDisplayBase_1.31.2
[85] ellipsis_0.3.2               
```
</div>