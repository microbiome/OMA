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
(`TSE`),
[`MicrobiomeExperiment`](https://github.com/FelixErnst/MicrobiomeExperiment)
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

The `phyloseq` package and class was around before the [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)  
and the derived 
[`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html) 
class. To make the transition easy 
here is short description how `phyloseq` and `*Experiment` classes relate to 
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
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: a phylo (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

### Assays  

The `assays` slots contains the experimental data as count matrices. Multiple 
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
## LinkDataFrame with 19216 rows and 4 columns
##            nodeLab nodeLab_alias   nodeNum    isLeaf
##        <character>   <character> <integer> <logical>
## 549322      549322       alias_1         1      TRUE
## 522457      522457       alias_2         2      TRUE
## 951            951       alias_3         3      TRUE
## 244423      244423       alias_4         4      TRUE
## 586076      586076       alias_5         5      TRUE
## ...            ...           ...       ...       ...
## 278222      278222   alias_19212     19212      TRUE
## 463590      463590   alias_19213     19213      TRUE
## 535321      535321   alias_19214     19214      TRUE
## 200359      200359   alias_19215     19215      TRUE
## 271582      271582   alias_19216     19216      TRUE
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
##    FeatureID SampleID relabundance Kingdom Phylum Class Order Family Genus
##    <fct>     <fct>           <dbl> <chr>   <chr>  <chr> <chr> <chr>  <chr>
##  1 549322    CL3                 0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  2 549322    CC1                 0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  3 549322    SV1                 0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  4 549322    M31Fcsw             0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  5 549322    M11Fcsw             0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  6 549322    M31Plmr             0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  7 549322    M11Plmr             0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  8 549322    F21Plmr             0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
##  9 549322    M31Tong             0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
## 10 549322    M11Tong             0 Archaea Crena… Ther… <NA>  <NA>   <NA> 
## # … with 499,606 more rows, and 8 more variables: Species <chr>,
## #   X.SampleID <fct>, Primer <fct>, Final_Barcode <fct>,
## #   Barcode_truncated_plus_T <fct>, Barcode_full_length <fct>,
## #   SampleType <fct>, Description <fct>
```

## Conclusion

Some wrapping up...

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
 [1] mia_0.0.0.9007                   MicrobiomeExperiment_0.99.0.9014
 [3] Biostrings_2.58.0                XVector_0.30.0                  
 [5] TreeSummarizedExperiment_1.6.0   SingleCellExperiment_1.12.0     
 [7] SummarizedExperiment_1.20.0      Biobase_2.50.0                  
 [9] GenomicRanges_1.42.0             GenomeInfoDb_1.26.1             
[11] IRanges_2.24.0                   S4Vectors_0.28.0                
[13] BiocGenerics_0.36.0              MatrixGenerics_1.2.0            
[15] matrixStats_0.57.0               BiocStyle_2.18.1                
[17] rebook_1.0.0                     BiocManager_1.30.10             

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5                ape_5.4-1                
 [3] lattice_0.20-41           tidyr_1.1.2              
 [5] ps_1.4.0                  utf8_1.1.4               
 [7] assertthat_0.2.1          digest_0.6.27            
 [9] R6_2.5.0                  evaluate_0.14            
[11] pillar_1.4.7              sparseMatrixStats_1.2.0  
[13] zlibbioc_1.36.0           rlang_0.4.9              
[15] callr_3.5.1               Matrix_1.2-18            
[17] rmarkdown_2.5             BiocParallel_1.24.1      
[19] stringr_1.4.0             RCurl_1.98-1.2           
[21] beachmat_2.6.2            DelayedArray_0.16.0      
[23] compiler_4.0.3            xfun_0.19                
[25] pkgconfig_2.0.3           CodeDepends_0.6.5        
[27] htmltools_0.5.0           tidyselect_1.1.0         
[29] tibble_3.0.4              GenomeInfoDbData_1.2.4   
[31] bookdown_0.21             codetools_0.2-18         
[33] XML_3.99-0.5              fansi_0.4.1              
[35] crayon_1.3.4              dplyr_1.0.2              
[37] bitops_1.0-6              grid_4.0.3               
[39] nlme_3.1-150              lifecycle_0.2.0          
[41] magrittr_2.0.1            graph_1.68.0             
[43] cli_2.2.0                 stringi_1.5.3            
[45] scuttle_1.0.3             DelayedMatrixStats_1.12.1
[47] ellipsis_0.3.1            generics_0.1.0           
[49] vctrs_0.3.5               tools_4.0.3              
[51] glue_1.4.2                purrr_0.3.4              
[53] processx_3.4.4            yaml_2.2.1               
[55] knitr_1.30               
```
</div>
