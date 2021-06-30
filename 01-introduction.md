# (PART) Introduction {-}

# Data Infrastructure {#data-introduction}

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

The
[`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)
(`SE`) is a widely used class for analyzing data obtained by common sequencing
techniques. `SE` is common structure for several Bioconductor packages that are
used for analyzing RNAseq, ChIp-Seq data. `SE` class is also used in R packages
for analyzing microarrays, flow cytometry, proteomics, single-cell sequencing
data and many more. The single-cell analysis is facilitated by
[SingelCellExperiment
class](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)
(`SCE`), which allows the user to store results of dimensionality reduction and
alternative experiments. Alternative experiments (`altExps`) can be differently
processed data within the analysis workflows.

Recently,
[`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)
(`TSE`)
were developed to extend the `SE` and `SCE` class for incorporating hierarchical
information (including phylogenetic tree) and reference sequences.

The `mia` package implements tools using these classes for analysis of
microbiome sequencing data.

## Installation

Install the development version from GitHub using `remotes` R package.  


```r
# install remotes 
#install.packages("remotes")
BiocManager::install("FelixErnst/mia")
```

### Packages    

1. `mia`    : Microbiome analysis tools   
2. `miaViz` : Microbiome analysis specific visualization

**See also:**    

[`microbiome`](https://bioconductor.org/packages/devel/bioc/html/microbiome.html)


## Background

The widely used `phyloseq` package and class were around before the [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)  
and the derived 
[`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html) 
class. Many methods for taxonomic profiling data are readily for the  `phyloseq` class structure. 

In order to facilitate the transition, we provide here a short description how `phyloseq` and `*Experiment` classes relate to 
each other.

`assays`     : This slot is similar to `otu_table` in `phyloseq`. In a `SummarizedExperiment`
               object multiple assays, raw counts, transformed counts can be stored. See also 
               [`MultiAssayExperiment`](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) for storing data from multiple `experiments` such as RNASeq, Proteomics, etc.       
`rowData`    : This slot is similar to `tax_table` in `phyloseq` to store taxonomic information.     
`colData`    : This slot is similar to `sample_data` in `phyloseq` to store information related to samples.    
`rowTree`    : This slot is similar to `phy_tree` in `phyloseq` to store phylogenetic tree.     

In this book, you will come across terms like `FeatureIDs` and `SampleIDs`.   
`FeatureIDs` : These are basically OTU/ASV ids which are row names in `assays` and `rowData`.    
`SampleIDs`  : As the name suggests, these are sample ids which are column names in `assays` and row names in `colData`.  

`FeatureIDs` and `SampleIDs` are used but the technical terms `rownames` and 
`colnames` are encouraged to be used, since they relate to actual objects we 
work with.



## Loading experimental microbiome data


### Importing data from external files

Microbiome (taxonomic) profiling data is commonly distributed in
various file formats. You can import such external data files as a
(Tree)SummarizedExperiment object but the details depend on the file
format. Here, we provide examples for common formats.

A specific import function is provided for QIIME2 files (see
`help(mia::loadFromQIIME2)`).

TODO: link to Mothur and other importers.

**CSV data tables** can be imported with the standard R functions,
  then converted to the desired format. For detailed examples, you can
  check the [Bioconductor course
  material](https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html)
  by Martin Morgan. The following example reads abundance tables,
  taxonomic mapping tables, and sample metadata, assuming that the
  input data files are properly prepared with appropriate row and
  column names.


```r
counts <- read.csv(count_file)   # Abundance table (e.g. ASV data; to assay data)
tax <- read.csv(tax_file)        # Taxonomy table (to rowData)
samples <- read.csv(sample_file) # Sample data (to colData)
se <- SummarizedExperiment(assays = list(counts = counts),
                           colData = samples,
                           rowData = tax)
```



### Conversions between data formats in R

If the data is has already been imported in R in another format, it
can be readily converted into `TreeSummarizedExperiment`, as shown in our next
example. Note that similar conversion functions to
`TreeSummarizedExperiment` are available for multiple data formats via
the `mia` package (see makeTreeSummarizedExperimentFrom* for phyloseq,
Biom, and DADA2).


```r
library(mia)

# phyloseq example data
data(GlobalPatterns, package="phyloseq") 
GlobalPatterns_phyloseq <- GlobalPatterns

# convert phyloseq to TSE
GlobalPatterns_TSE <- makeTreeSummarizedExperimentFromphyloseq(GlobalPatterns_phyloseq) 
```


We can also convert `TreeSummarizedExperiment` objects into `phyloseq`
with respect to the shared components that are supported by both
formats (i.e. taxonomic abundance table, sample metadata, taxonomic
table, phylogenetic tree, sequence information). This is useful for
instance when additional methods are available for `phyloseq`.

TODO: conversion function from TSE to phyloseq





## Metadata

## Microbiome and tree data specific aspects


```r
library(mia)
data("GlobalPatterns")
se <- GlobalPatterns 
se
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

### Assays  

The `assays` slot contains the experimental data as count matrices. Multiple 
matrices can be stored the result of `assays` is actually a list of matrices.


```r
assays(se)
```

```
## List of length 1
## names(1): counts
```

Individual assays can be accessed via `assay`


```r
assay(se, "counts")[1:5,1:7]
```

```
##        CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr M11Plmr
## 549322   0   0   0       0       0       0       0
## 522457   0   0   0       0       0       0       0
## 951      0   0   0       0       0       0       1
## 244423   0   0   0       0       0       0       0
## 586076   0   0   0       0       0       0       0
```

To illustrate the use of multiple assays, the relative abundance data can be 
calcualted and stored along the original count data using `relAbundanceCounts`.


```r
se <- relAbundanceCounts(se)
assays(se)
```

```
## List of length 2
## names(2): counts relabundance
```

Now there are two assays available in the `se` object, `counts` and 
`relabundance`.


```r
assay(se, "relabundance")[1:5,1:7]
```

```
##        CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr   M11Plmr
## 549322   0   0   0       0       0       0 0.000e+00
## 522457   0   0   0       0       0       0 0.000e+00
## 951      0   0   0       0       0       0 2.305e-06
## 244423   0   0   0       0       0       0 0.000e+00
## 586076   0   0   0       0       0       0 0.000e+00
```

### colData

`colData` contains data on the samples.


```r
colData(se)
```

```
## DataFrame with 26 rows and 7 columns
##         X.SampleID   Primer Final_Barcode Barcode_truncated_plus_T
##           <factor> <factor>      <factor>                 <factor>
## CL3        CL3      ILBC_01        AACGCA                   TGCGTT
## CC1        CC1      ILBC_02        AACTCG                   CGAGTT
## SV1        SV1      ILBC_03        AACTGT                   ACAGTT
## M31Fcsw    M31Fcsw  ILBC_04        AAGAGA                   TCTCTT
## M11Fcsw    M11Fcsw  ILBC_05        AAGCTG                   CAGCTT
## ...            ...      ...           ...                      ...
## TS28         TS28   ILBC_25        ACCAGA                   TCTGGT
## TS29         TS29   ILBC_26        ACCAGC                   GCTGGT
## Even1        Even1  ILBC_27        ACCGCA                   TGCGGT
## Even2        Even2  ILBC_28        ACCTCG                   CGAGGT
## Even3        Even3  ILBC_29        ACCTGT                   ACAGGT
##         Barcode_full_length SampleType
##                    <factor>   <factor>
## CL3             CTAGCGTGCGT      Soil 
## CC1             CATCGACGAGT      Soil 
## SV1             GTACGCACAGT      Soil 
## M31Fcsw         TCGACATCTCT      Feces
## M11Fcsw         CGACTGCAGCT      Feces
## ...                     ...        ...
## TS28            GCATCGTCTGG      Feces
## TS29            CTAGTCGCTGG      Feces
## Even1           TGACTCTGCGG      Mock 
## Even2           TCTGATCGAGG      Mock 
## Even3           AGAGAGACAGG      Mock 
##                                        Description
##                                           <factor>
## CL3     Calhoun South Carolina Pine soil, pH 4.9  
## CC1     Cedar Creek Minnesota, grassland, pH 6.1  
## SV1     Sevilleta new Mexico, desert scrub, pH 8.3
## M31Fcsw M3, Day 1, fecal swab, whole body study   
## M11Fcsw M1, Day 1, fecal swab, whole body study   
## ...                                            ...
## TS28                                       Twin #1
## TS29                                       Twin #2
## Even1                                      Even1  
## Even2                                      Even2  
## Even3                                      Even3
```

### rowData

`rowData` contains data on the features of the samples analyzed. Of particular
interest for the microbiome field this is used to store taxonomic information.


```r
rowData(se)
```

```
## DataFrame with 19216 rows and 7 columns
##            Kingdom        Phylum        Class        Order        Family
##        <character>   <character>  <character>  <character>   <character>
## 549322     Archaea Crenarchaeota Thermoprotei           NA            NA
## 522457     Archaea Crenarchaeota Thermoprotei           NA            NA
## 951        Archaea Crenarchaeota Thermoprotei Sulfolobales Sulfolobaceae
## 244423     Archaea Crenarchaeota        Sd-NA           NA            NA
## 586076     Archaea Crenarchaeota        Sd-NA           NA            NA
## ...            ...           ...          ...          ...           ...
## 278222    Bacteria           SR1           NA           NA            NA
## 463590    Bacteria           SR1           NA           NA            NA
## 535321    Bacteria           SR1           NA           NA            NA
## 200359    Bacteria           SR1           NA           NA            NA
## 271582    Bacteria           SR1           NA           NA            NA
##              Genus                Species
##        <character>            <character>
## 549322          NA                     NA
## 522457          NA                     NA
## 951     Sulfolobus Sulfolobusacidocalda..
## 244423          NA                     NA
## 586076          NA                     NA
## ...            ...                    ...
## 278222          NA                     NA
## 463590          NA                     NA
## 535321          NA                     NA
## 200359          NA                     NA
## 271582          NA                     NA
```

### rowTree  

Phylogenetic trees also play an important role for the microbiome field. The 
`TreeSummarizedExperiment` class is able to keep track of feature and node
relations via two functions, `rowTree` and `rowLinks`.

A tree can be accessed via `rowTree` as `phylo` object.       

```r
rowTree(se)
```

```
## 
## Phylogenetic tree with 19216 tips and 19215 internal nodes.
## 
## Tip labels:
##   549322, 522457, 951, 244423, 586076, 246140, ...
## Node labels:
##   , 0.858.4, 1.000.154, 0.764.3, 0.995.2, 1.000.2, ...
## 
## Rooted; includes branch lengths.
```

The links to the individual features are available through `rowLinks`.


```r
rowLinks(se)
```

```
## LinkDataFrame with 19216 rows and 5 columns
##           nodeLab   nodeNum nodeLab_alias    isLeaf   whichTree
##       <character> <integer>   <character> <logical> <character>
## 1          549322         1       alias_1      TRUE       phylo
## 2          522457         2       alias_2      TRUE       phylo
## 3             951         3       alias_3      TRUE       phylo
## 4          244423         4       alias_4      TRUE       phylo
## 5          586076         5       alias_5      TRUE       phylo
## ...           ...       ...           ...       ...         ...
## 19212      278222     19212   alias_19212      TRUE       phylo
## 19213      463590     19213   alias_19213      TRUE       phylo
## 19214      535321     19214   alias_19214      TRUE       phylo
## 19215      200359     19215   alias_19215      TRUE       phylo
## 19216      271582     19216   alias_19216      TRUE       phylo
```

Please note that there can be a 1:1 relationship between tree nodes and 
features, but this is not a must have. This means there can be features, which
are not linked to nodes, and nodes, which are not linked to features. To change
the links in an existing object, the `changeTree` function is available.

## Data conversion

Sometimes custom solutions are need for analyzing the data. `mia` contains a 
few functions to help in these situations.

### Tidy data

For several custom analysis and visualization packages such as those from the 
`tidyverse` the `SE` data can be converted to long data.frame format with 
`meltAssay`.    


```r
library(mia)
molten_se <- meltAssay(se,
                       add_row_data = TRUE,
                       add_col_data = TRUE,
                       abund_values = "relabundance")
molten_se
```

```
## # A tibble: 499,616 x 17
##    FeatureID SampleID counts Kingdom Phylum   Class   Order Family Genus Species
##    <fct>     <fct>     <dbl> <chr>   <chr>    <chr>   <chr> <chr>  <chr> <chr>  
##  1 549322    CL3           0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  2 549322    CC1           0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  3 549322    SV1           0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  4 549322    M31Fcsw       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  5 549322    M11Fcsw       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  6 549322    M31Plmr       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  7 549322    M11Plmr       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  8 549322    F21Plmr       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
##  9 549322    M31Tong       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
## 10 549322    M11Tong       0 Archaea Crenarc… Thermo… <NA>  <NA>   <NA>  <NA>   
## # … with 499,606 more rows, and 7 more variables: X.SampleID <fct>,
## #   Primer <fct>, Final_Barcode <fct>, Barcode_truncated_plus_T <fct>,
## #   Barcode_full_length <fct>, SampleType <fct>, Description <fct>
```

## Conclusion

Some wrapping up...

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
 [1] mia_1.1.5                      TreeSummarizedExperiment_2.1.3
 [3] Biostrings_2.61.1              XVector_0.33.0                
 [5] SingleCellExperiment_1.15.1    SummarizedExperiment_1.23.1   
 [7] Biobase_2.53.0                 GenomicRanges_1.45.0          
 [9] GenomeInfoDb_1.29.2            IRanges_2.27.0                
[11] S4Vectors_0.31.0               BiocGenerics_0.39.1           
[13] MatrixGenerics_1.5.0           matrixStats_0.59.0            
[15] BiocStyle_2.21.3               rebook_1.3.0                  

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0            colorspace_2.0-2           
  [3] ellipsis_0.3.2              scuttle_1.3.0              
  [5] BiocNeighbors_1.11.0        rstudioapi_0.13            
  [7] bit64_4.0.5                 fansi_0.5.0                
  [9] decontam_1.13.0             splines_4.1.0              
 [11] codetools_0.2-18            sparseMatrixStats_1.5.0    
 [13] cachem_1.0.5                knitr_1.33                 
 [15] scater_1.21.2               ade4_1.7-17                
 [17] phyloseq_1.37.0             jsonlite_1.7.2             
 [19] cluster_2.1.2               graph_1.71.2               
 [21] BiocManager_1.30.16         compiler_4.1.0             
 [23] assertthat_0.2.1            Matrix_1.3-4               
 [25] fastmap_1.1.0               lazyeval_0.2.2             
 [27] cli_2.5.0                   BiocSingular_1.9.1         
 [29] htmltools_0.5.1.1           tools_4.1.0                
 [31] igraph_1.2.6                rsvd_1.0.5                 
 [33] gtable_0.3.0                glue_1.4.2                 
 [35] GenomeInfoDbData_1.2.6      reshape2_1.4.4             
 [37] dplyr_1.0.7                 Rcpp_1.0.6                 
 [39] jquerylib_0.1.4             rhdf5filters_1.5.0         
 [41] vctrs_0.3.8                 multtest_2.49.0            
 [43] ape_5.5                     nlme_3.1-152               
 [45] DECIPHER_2.21.0             iterators_1.0.13           
 [47] DelayedMatrixStats_1.15.0   xfun_0.24                  
 [49] stringr_1.4.0               ps_1.6.0                   
 [51] beachmat_2.9.0              lifecycle_1.0.0            
 [53] irlba_2.3.3                 XML_3.99-0.6               
 [55] zlibbioc_1.39.0             MASS_7.3-54                
 [57] scales_1.1.1                biomformat_1.21.0          
 [59] parallel_4.1.0              rhdf5_2.37.0               
 [61] yaml_2.2.1                  memoise_2.0.0              
 [63] gridExtra_2.3               ggplot2_3.3.5              
 [65] sass_0.4.0                  stringi_1.6.2              
 [67] RSQLite_2.2.7               foreach_1.5.1              
 [69] ScaledMatrix_1.1.0          tidytree_0.3.4             
 [71] permute_0.9-5               filelock_1.0.2             
 [73] BiocParallel_1.27.0         rlang_0.4.11               
 [75] pkgconfig_2.0.3             bitops_1.0-7               
 [77] evaluate_0.14               lattice_0.20-44            
 [79] Rhdf5lib_1.15.1             purrr_0.3.4                
 [81] treeio_1.17.2               CodeDepends_0.6.5          
 [83] bit_4.0.4                   tidyselect_1.1.1           
 [85] plyr_1.8.6                  magrittr_2.0.1             
 [87] bookdown_0.22               R6_2.5.0                   
 [89] generics_0.1.0              DelayedArray_0.19.1        
 [91] DBI_1.1.1                   mgcv_1.8-36                
 [93] pillar_1.6.1                survival_3.2-11            
 [95] RCurl_1.98-1.3              tibble_3.1.2               
 [97] dir.expiry_1.1.0            crayon_1.4.1               
 [99] utf8_1.2.1                  rmarkdown_2.9              
[101] viridis_0.6.1               grid_4.1.0                 
[103] data.table_1.14.0           blob_1.2.1                 
[105] vegan_2.5-7                 digest_0.6.27              
[107] tidyr_1.1.3                 munsell_0.5.0              
[109] DirichletMultinomial_1.35.0 beeswarm_0.4.0             
[111] viridisLite_0.4.0           vipor_0.4.5                
[113] bslib_0.2.5.1              
```
</div>
