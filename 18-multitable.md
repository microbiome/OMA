# Representation of multiple data tables {#multitable}

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

Microbiome data can be part of multiomics experiments and analysis strategies
and we want to outline the understanding in which we think the packages 
explained and used in this book relate to these experiment layouts
using the `TreeSummarizedExperiment` and classes beyond.

Many microbiome experiments include multiple versions and types of
data generated independently or derived from each other through transformation
or agglomeration. We start by providing recommendations on how to represent
different varieties of multi-table data within the
`TreeSummarizedExperiment` class.

The options and recommendations are summarized in Table \@ref(tab:options).


## Assay data

The original count-based taxonomic abundance tables may have different 
transformations, such as logarithmic, Centered Log-Ratio (CLR), or relative 
abundance. These are typically stored in _**assays**_.


```r
library(mia)
data(GlobalPatterns)
se <- GlobalPatterns
assays(se)
```

```
## List of length 1
## names(1): counts
```

As an example the relative abundance is calculated.


```r
se <- relAbundanceCounts(se)
assays(se)
```

```
## List of length 2
## names(2): counts relabundance
```

Here the dimension of the count data remains unchanged. This is
actually a requirement for any `SummarizedExperiment` object.


## Alternative experiments

_**Alternative experiments**_ differ from transformations as they can
contain complementary data, which is no longer tied to the same
dimensions as the assay data. However, the number of samples (columns)
must be the same, however.

This can come into play for instance when one has taxonomic abundance
profiles quantified with different measurement technologies, such as
phylogenetic microarrays, amplicon sequencing, or metagenomic
sequencing. Such alternative experiments that concern the same samples
can be stored as

1. Separate _assays_ assuming that the taxonomic information can be mapped 
between feature directly 1:1; or 
2. data in the _altExp_ slot of the `TreeSummarizedExperiment`, if the feature 
dimensions differ. Each element of the _altExp_ slot is a `SummarizedExperiment`
or an object from a derived class with independent feature data.


As an example, we show how to store taxonomic abundance tables
agglomerated at different taxonomic levels. However, the data could as
well originate from entirely different measurement sources as long as
the samples are matched.


```r
# Agglomerate the data to Phylym level
se.phylum <- agglomerateByRank(se, "Phylum")
# both have the same number of columns (samples)
dim(se)
```

```
## [1] 19216    26
```

```r
dim(se.phylum)
```

```
## [1] 67 26
```

```r
# Add the new table as an alternative experiment
altExp(se, "Phylum") <- se.phylum
altExpNames(se)
```

```
## [1] "Phylum"
```

```r
# Pick a sample subset: this acts on both altExp and assay data
se[,1:10]
```

```
## class: TreeSummarizedExperiment 
## dim: 19216 10 
## metadata(0):
## assays(2): counts relabundance
## rownames(19216): 549322 522457 ... 200359 271582
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(10): CL3 CC1 ... M31Tong M11Tong
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(1): Phylum
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

```r
dim(altExp(se[,1:10],"Phylum"))
```

```
## [1] 67 10
```

For more details of altExp have a look at the [Intro vignette](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html) of the 
`SingleCellExperiment` package [@R-SingleCellExperiment].



## MultiAssayExperiments

_**Multiple experiments**_ relate to complementary measurement types,
such as transcriptomic or metabolomic profiling of the microbiome or
the host. Multiple experiments can be represented using the same
options as alternative experiments, or by using the
`MultiAssayExperiment` class [@R-MultiAssayExperiment]. Depending on how the 
datasets relate to each other the data can be stored as:

1. Separate _altExp_ if the samples can be matched directly 1:1; or
2. As `MultiAssayExperiment` objects, in which the connections between
sample are defined through a `sampleMap`. Each element on the
`experimentsList` of an `MultiAssayExperiment` is `matrix` or
`matrix`-like object including `SummarizedExperiment` objects, and the
number of samples can differ between the elements.



```r
#TODO: Find the right dataset to explain a non 1:1 sample relationship
```


For information have a look at the [intro vignette](https://bioconductor.org/packages/release/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html) of the `MultiAssayExperiment` package.  

 
   Option   Rows (features)    Cols (samples)               Recommended  
---------   --------------    ---------------  ------------------------
   assays  	     match              match       Data transformations  
   altExp             free              match    Alternative experiments  
MultiAssay            free      free (mapping)    Multi-omic experiments    

Table: (\#tab:options) **Recommended options for storing multiple data tables in microbiome studies** The _assays_ are best suited for data transformations (one-to-one match between samples and columns across the assays). The _alternative experiments_ are particularly suitable for alternative versions of the data that is of same type but may have a different number of features (e.g. taxonomic groups); this is for instance the case with taxonomic abundance tables agglomerated at different levels (e.g. genus vs. phyla) or alternative profiling technologies (e.g. amplicon sequencing vs. shallow shotgun metagenomics). For alternative experiments one-to-one match between samples (cols) is required but the alternative experiment tables can have different numbers of features (rows). Finally, elements of the _MultiAssayExperiment_ provide the most flexible way to incorporate multi-omic data tables with flexible numbers of samples and features. We recommend these conventions as the basis for methods development and application in microbiome studies.




## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R Under development (unstable) (2021-03-18 r80099)
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
 [1] mia_0.99.5                       TreeSummarizedExperiment_1.99.11
 [3] Biostrings_2.59.2                XVector_0.31.1                  
 [5] SingleCellExperiment_1.13.12     SummarizedExperiment_1.21.1     
 [7] Biobase_2.51.0                   GenomicRanges_1.43.3            
 [9] GenomeInfoDb_1.27.8              IRanges_2.25.6                  
[11] S4Vectors_0.29.9                 BiocGenerics_0.37.1             
[13] MatrixGenerics_1.3.1             matrixStats_0.58.0              
[15] BiocStyle_2.19.1                 rebook_1.1.16                   
[17] BiocManager_1.30.10             

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
[19] processx_3.5.0              bit_4.0.4                  
[21] compiler_4.1.0              graph_1.69.0               
[23] BiocNeighbors_1.9.4         DelayedArray_0.17.9        
[25] bookdown_0.21               sass_0.3.1                 
[27] scales_1.1.1                callr_3.5.1                
[29] stringr_1.4.0               digest_0.6.27              
[31] rmarkdown_2.7               scater_1.19.11             
[33] pkgconfig_2.0.3             htmltools_0.5.1.1          
[35] sparseMatrixStats_1.3.6     fastmap_1.1.0              
[37] rlang_0.4.10                RSQLite_2.2.4              
[39] DelayedMatrixStats_1.13.5   jquerylib_0.1.3            
[41] generics_0.1.0              jsonlite_1.7.2             
[43] BiocParallel_1.25.5         dplyr_1.0.5                
[45] RCurl_1.98-1.3              magrittr_2.0.1             
[47] BiocSingular_1.7.2          GenomeInfoDbData_1.2.4     
[49] scuttle_1.1.18              Matrix_1.3-2               
[51] Rcpp_1.0.6                  ggbeeswarm_0.6.0           
[53] munsell_0.5.0               fansi_0.4.2                
[55] DECIPHER_2.19.2             viridis_0.5.1              
[57] ape_5.4-1                   lifecycle_1.0.0            
[59] stringi_1.5.3               yaml_2.2.1                 
[61] MASS_7.3-53.1               debugme_1.1.0              
[63] zlibbioc_1.37.0             blob_1.2.1                 
[65] grid_4.1.0                  crayon_1.4.1               
[67] lattice_0.20-41             splines_4.1.0              
[69] beachmat_2.7.7              CodeDepends_0.6.5          
[71] knitr_1.31                  ps_1.6.0                   
[73] pillar_1.5.1                codetools_0.2-18           
[75] ScaledMatrix_0.99.2         XML_3.99-0.6               
[77] glue_1.4.2                  evaluate_0.14              
[79] vctrs_0.3.6                 tidyr_1.1.3                
[81] gtable_0.3.0                purrr_0.3.4                
[83] assertthat_0.2.1            cachem_1.0.4               
[85] ggplot2_3.3.3               xfun_0.22                  
[87] rsvd_1.0.3                  viridisLite_0.3.0          
[89] tibble_3.1.0                memoise_2.0.0              
[91] beeswarm_0.3.1              cluster_2.1.1              
[93] ellipsis_0.3.1             
```
</div>
