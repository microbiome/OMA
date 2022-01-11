# Taxonomic Information {#taxonomic-information}

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


```r
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns 
```

Taxonomic information is a key part of analyzing microbiome data and without
it, any type of data analysis probably will not make much sense. However,
the degree of detail of taxonomic information differs depending on the dataset
and annotation data used.

Therefore, the mia package expects a loose assembly of taxonomic information
and assumes certain key aspects:

* Taxonomic information is given as character vectors or factors in the 
`rowData` of a `SummarizedExperiment` object.
* The columns containing the taxonomic information must be named `domain`,
`kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` or with
a capital first letter.
* the columns must be given in the order shown above
* column can be omited, but the order must remain

## Assigning taxonomic information.

There are a number of methods to assign taxonomic information. We like to give
a short introduction about the methods available without ranking one over the 
other. This has to be your choice based on the result for the individual 
dataset.

### dada2

The dada2 package [@R-dada2] implements the `assignTaxonomy` function, which 
takes as input the ASV sequences associated with each row of data and a training
dataset. For more information visit the 
[dada2 homepage](https://benjjneb.github.io/dada2/assign.html).

### DECIPHER

The DECIPHER package [@R-DECIPHER] implements the `IDTAXA` algorithm to assign
either taxonomic information or function information. For `mia`
only the first option is of interest for now and more information can be
found on the [DECIPHER website](http://www2.decipher.codes/Classification.html).

## Functions to access taxonomic information

`checkTaxonomy` checks whether the taxonomic information is usable for `mia`


```r
checkTaxonomy(tse)
```

```
## [1] TRUE
```

Since the `rowData` can contain other data, `taxonomyRanks` will return the 
columns `mia` assumes to contain the taxonomic information.


```r
taxonomyRanks(tse)
```

```
## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
```

This can then be used to subset the `rowData` to columns needed.


```r
rowData(tse)[,taxonomyRanks(tse)]
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

`taxonomyRankEmpty` checks for empty values in the given `rank` and returns a 
logical vector of `length(x)`.


```r
all(!taxonomyRankEmpty(tse, rank = "Kingdom"))
```

```
## [1] TRUE
```

```r
table(taxonomyRankEmpty(tse, rank = "Genus"))
```

```
## 
## FALSE  TRUE 
##  8008 11208
```

```r
table(taxonomyRankEmpty(tse, rank = "Species"))
```

```
## 
## FALSE  TRUE 
##  1413 17803
```

`getTaxonomyLabels` is a multi-purpose function, which turns taxonomic
information into a character vector of `length(x)`


```r
head(getTaxonomyLabels(tse))
```

```
## [1] "Class:Thermoprotei"               "Class:Thermoprotei_1"            
## [3] "Species:Sulfolobusacidocaldarius" "Class:Sd-NA"                     
## [5] "Class:Sd-NA_1"                    "Class:Sd-NA_2"
```

By default, this will use the lowest non-empty information to construct a
string with the following scheme `level:value`. If all levels are the same,
this part is omitted, but can be added by setting `with_rank = TRUE`.


```r
phylum <- !is.na(rowData(tse)$Phylum) & 
    vapply(data.frame(apply(rowData(tse)[,taxonomyRanks(tse)[3:7]],1L,is.na)),all,logical(1))
head(getTaxonomyLabels(tse[phylum,]))
```

```
## [1] "Crenarchaeota"    "Crenarchaeota_1"  "Crenarchaeota_2"  "Actinobacteria"  
## [5] "Actinobacteria_1" "Spirochaetes"
```

```r
head(getTaxonomyLabels(tse[phylum,], with_rank = TRUE))
```

```
## [1] "Phylum:Crenarchaeota"    "Phylum:Crenarchaeota_1" 
## [3] "Phylum:Crenarchaeota_2"  "Phylum:Actinobacteria"  
## [5] "Phylum:Actinobacteria_1" "Phylum:Spirochaetes"
```

By default the return value of `getTaxonomyLabels` contains only unique elements
by passing it through `make.unique`. This step can be omitted by setting 
`make_unique = FALSE`.


```r
head(getTaxonomyLabels(tse[phylum,], with_rank = TRUE, make_unique = FALSE))
```

```
## [1] "Phylum:Crenarchaeota"  "Phylum:Crenarchaeota"  "Phylum:Crenarchaeota" 
## [4] "Phylum:Actinobacteria" "Phylum:Actinobacteria" "Phylum:Spirochaetes"
```

To apply the loop resolving function `resolveLoop` from the
`TreeSummarizedExperiment` package [@R-TreeSummarizedExperiment] within
`getTaxonomyLabels`, set `resolve_loops = TRUE`.

### Generate a taxonomic tree on the fly

To create a taxonomic tree, `taxonomyTree` used the information and returns a
`phylo` object. Duplicate information from the `rowData` is removed.


```r
taxonomyTree(tse)
```

```
## 
## Phylogenetic tree with 1645 tips and 1089 internal nodes.
## 
## Tip labels:
##   Species:Cenarchaeumsymbiosum, Species:pIVWA5, Species:CandidatusNitrososphaeragargensis, Species:SCA1145, Species:SCA1170, Species:Sulfolobusacidocaldarius, ...
## Node labels:
##   root:ALL, Kingdom:Archaea, Phylum:Crenarchaeota, Class:C2, Class:Sd-NA, Class:Thaumarchaeota, ...
## 
## Rooted; includes branch lengths.
```


```r
tse <- addTaxonomyTree(tse)
tse
```

```
## class: TreeSummarizedExperiment 
## dim: 19216 26 
## metadata(0):
## assays(1): counts
## rownames(19216): Class:Thermoprotei Class:Thermoprotei ... Phylum:SR1
##   Phylum:SR1
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: 1 phylo tree(s) (1645 leaves)
## colLinks: NULL
## colTree: NULL
```

The implementation is based on the `toTree` function from the
`TreeSummarizedExperiment` package [@R-TreeSummarizedExperiment].

## Data agglomeration {#data-agglomeration}

One of the main applications of taxonomic information in regards to count data
is to agglomerate count data on taxonomic levels and track the influence of 
changing conditions through these levels. For this `mia` contains the
`agglomerateByRank` function. The ideal location to store the agglomerated data
is as an alternative experiment.


```r
tse <- relAbundanceCounts(tse)
altExp(tse, "Family") <- agglomerateByRank(tse, rank = "Family",
                                           agglomerateTree = TRUE)
altExp(tse, "Family")
```

```
## class: TreeSummarizedExperiment 
## dim: 603 26 
## metadata(0):
## assays(2): counts relabundance
## rownames(603): Class:Thermoprotei Family:Sulfolobaceae ...
##   Family:Thermodesulfobiaceae Phylum:SR1
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (603 rows)
## rowTree: 1 phylo tree(s) (496 leaves)
## colLinks: NULL
## colTree: NULL
```

If multiple assays (counts and relabundance) exist, both will be agglomerated.


```r
assayNames(tse)
```

```
## [1] "counts"       "relabundance"
```

```r
assayNames(altExp(tse, "Family"))
```

```
## [1] "counts"       "relabundance"
```


```r
assay(altExp(tse, "Family"), "relabundance")[1:5,1:7]
```

```
##                            CL3       CC1 SV1 M31Fcsw M11Fcsw M31Plmr   M11Plmr
## Class:Thermoprotei   0.0000000 0.000e+00   0       0       0       0 0.000e+00
## Family:Sulfolobaceae 0.0000000 0.000e+00   0       0       0       0 2.305e-06
## Class:Sd-NA          0.0000000 0.000e+00   0       0       0       0 0.000e+00
## Order:NRP-J          0.0001991 2.070e-04   0       0       0       0 6.914e-06
## Family:SAGMA-X       0.0000000 6.165e-06   0       0       0       0 0.000e+00
```
  

```r
assay(altExp(tse, "Family"), "counts")[1:5,1:7]
```

```
##                      CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr M11Plmr
## Class:Thermoprotei     0   0   0       0       0       0       0
## Family:Sulfolobaceae   0   0   0       0       0       0       1
## Class:Sd-NA            0   0   0       0       0       0       0
## Order:NRP-J          172 235   0       0       0       0       3
## Family:SAGMA-X         0   7   0       0       0       0       0
```

`altExpNames` now consists of `Family` level data. This can be extended to use 
any level present in Kingdom, Phylum, Class, Order, Family, Genus, Species.   

## Data transformation

Data transformation is a very common procedure in microbiome analysis. 
In transformation, each data point is replaced with transformed value that is 
calculated by applying transformation formula to the data point. Transformation 
can be used, for example, to normalize skewed data, or to reduce weight of bigger 
values compared to smaller values. 

In mia package, transformations are applied to abundance data. The transformed 
abundance table is stored back to 'assays'. mia includes transformation 
functions for sample-wise or column-wise transformation ('transformSamples()'), 
and for feature-wise or row-wise transformation ('transformFeatures()'). 

For complete list of available transformations and parameters, see function 
[help](https://microbiome.github.io/mia/reference/transformCounts.html).


```r
tse <- transformSamples(x = tse, abund_values = "counts", method = "clr", 
                        pseudocount = 1, name = "clr_transformation")

head(assay(tse, "clr_transformation"))
```

```
##                                      CL3    CC1     SV1 M31Fcsw M11Fcsw M31Plmr
## Class:Thermoprotei               -0.9552 -1.124 -0.7435 -0.2916 -0.2652  -0.356
## Class:Thermoprotei               -0.9552 -1.124 -0.7435 -0.2916 -0.2652  -0.356
## Species:Sulfolobusacidocaldarius -0.9552 -1.124 -0.7435 -0.2916 -0.2652  -0.356
## Class:Sd-NA                      -0.9552 -1.124 -0.7435 -0.2916 -0.2652  -0.356
## Class:Sd-NA                      -0.9552 -1.124 -0.7435 -0.2916 -0.2652  -0.356
## Class:Sd-NA                      -0.9552 -1.124 -0.7435 -0.2916 -0.2652  -0.356
##                                  M11Plmr F21Plmr M31Tong M11Tong LMEpi24M
## Class:Thermoprotei               -0.4713 -0.2645 -0.2547 -0.1572   -0.359
## Class:Thermoprotei               -0.4713 -0.2645 -0.2547 -0.1572   -0.359
## Species:Sulfolobusacidocaldarius  0.2219 -0.2645 -0.2547 -0.1572   -0.359
## Class:Sd-NA                      -0.4713 -0.2645 -0.2547 -0.1572   -0.359
## Class:Sd-NA                      -0.4713 -0.2645 -0.2547 -0.1572   -0.359
## Class:Sd-NA                      -0.4713 -0.2645 -0.2547 -0.1572   -0.359
##                                  SLEpi20M  AQC1cm  AQC4cm  AQC7cm     NP2
## Class:Thermoprotei                 0.3704  2.6250  3.7862  4.0751  0.4502
## Class:Thermoprotei                -0.3228 -0.7072  0.2697  1.1459 -0.2429
## Species:Sulfolobusacidocaldarius  -0.3228 -0.7072 -0.8289 -0.8001 -0.2429
## Class:Sd-NA                       -0.3228 -0.7072  2.3066  2.6011 -0.2429
## Class:Sd-NA                       -0.3228 -0.7072  0.2697 -0.1069 -0.2429
## Class:Sd-NA                       -0.3228 -0.7072 -0.1357  0.5862 -0.2429
##                                     NP3     NP5 TRRsed1 TRRsed2 TRRsed3    TS28
## Class:Thermoprotei               -0.433 -0.3606 -0.2677 -0.4828 -0.4384 -0.2691
## Class:Thermoprotei               -0.433 -0.3606 -0.2677 -0.4828 -0.4384 -0.2691
## Species:Sulfolobusacidocaldarius -0.433 -0.3606 -0.2677 -0.4828 -0.4384 -0.2691
## Class:Sd-NA                      -0.433 -0.3606 -0.2677 -0.4828 -0.4384 -0.2691
## Class:Sd-NA                      -0.433 -0.3606 -0.2677 -0.4828 -0.4384 -0.2691
## Class:Sd-NA                      -0.433 -0.3606 -0.2677 -0.4828 -0.4384 -0.2691
##                                     TS29   Even1   Even2   Even3
## Class:Thermoprotei               -0.2569 -0.3481 -0.2534 -0.2382
## Class:Thermoprotei               -0.2569 -0.3481 -0.2534 -0.2382
## Species:Sulfolobusacidocaldarius -0.2569 -0.3481 -0.2534 -0.2382
## Class:Sd-NA                      -0.2569 -0.3481 -0.2534 -0.2382
## Class:Sd-NA                      -0.2569 -0.3481 -0.2534 -0.2382
## Class:Sd-NA                      -0.2569 -0.3481 -0.2534 -0.2382
```

-   In 'pa' transformation, 'threshold' specifies the value that divides observations to
be absent or present. By default, it is 0.


```r
tse <- transformFeatures(tse, method = "pa", threshold = 10)

head(assay(tse, "pa"))
```

```
##                                  CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr M11Plmr
## Class:Thermoprotei                 0   0   0       0       0       0       0
## Class:Thermoprotei                 0   0   0       0       0       0       0
## Species:Sulfolobusacidocaldarius   0   0   0       0       0       0       0
## Class:Sd-NA                        0   0   0       0       0       0       0
## Class:Sd-NA                        0   0   0       0       0       0       0
## Class:Sd-NA                        0   0   0       0       0       0       0
##                                  F21Plmr M31Tong M11Tong LMEpi24M SLEpi20M
## Class:Thermoprotei                     0       0       0        0        0
## Class:Thermoprotei                     0       0       0        0        0
## Species:Sulfolobusacidocaldarius       0       0       0        0        0
## Class:Sd-NA                            0       0       0        0        0
## Class:Sd-NA                            0       0       0        0        0
## Class:Sd-NA                            0       0       0        0        0
##                                  AQC1cm AQC4cm AQC7cm NP2 NP3 NP5 TRRsed1
## Class:Thermoprotei                    1      1      1   0   0   0       0
## Class:Thermoprotei                    0      0      0   0   0   0       0
## Species:Sulfolobusacidocaldarius      0      0      0   0   0   0       0
## Class:Sd-NA                           0      1      1   0   0   0       0
## Class:Sd-NA                           0      0      0   0   0   0       0
## Class:Sd-NA                           0      0      0   0   0   0       0
##                                  TRRsed2 TRRsed3 TS28 TS29 Even1 Even2 Even3
## Class:Thermoprotei                     0       0    0    0     0     0     0
## Class:Thermoprotei                     0       0    0    0     0     0     0
## Species:Sulfolobusacidocaldarius       0       0    0    0     0     0     0
## Class:Sd-NA                            0       0    0    0     0     0     0
## Class:Sd-NA                            0       0    0    0     0     0     0
## Class:Sd-NA                            0       0    0    0     0     0     0
```


```r
# list of abundance tables that assays slot contains
assays(tse)
```

```
## List of length 4
## names(4): counts relabundance clr_transformation pa
```

## Pick specific  

Retrieving of specific elements that are required for specific analysis. For
instance, extracting abundances for a specific taxa in all samples or all taxa 
in one sample.  

### Abundances of all taxa in specific sample 

```r
taxa.abund.cc1 <- getAbundanceSample(tse, 
                                     sample_id = "CC1",
                                     abund_values = "counts")
taxa.abund.cc1[1:10]
```

```
##               Class:Thermoprotei               Class:Thermoprotei 
##                                0                                0 
## Species:Sulfolobusacidocaldarius                      Class:Sd-NA 
##                                0                                0 
##                      Class:Sd-NA                      Class:Sd-NA 
##                                0                                0 
##                      Order:NRP-J                      Order:NRP-J 
##                                1                                0 
##                      Order:NRP-J                      Order:NRP-J 
##                              194                                5
```

### Abundances of specific taxa in all samples   


```r
taxa.abundances <- getAbundanceFeature(tse, 
                                      feature_id = "Phylum:Bacteroidetes",
                                      abund_values = "counts")
taxa.abundances[1:10]
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr M11Plmr F21Plmr M31Tong M11Tong 
##       2      18       2       0       0       0       0       1       0       0
```


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.2 (2021-11-01)
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
 [1] mia_1.3.13                     MultiAssayExperiment_1.20.0   
 [3] TreeSummarizedExperiment_2.1.4 Biostrings_2.62.0             
 [5] XVector_0.34.0                 SingleCellExperiment_1.16.0   
 [7] SummarizedExperiment_1.24.0    Biobase_2.54.0                
 [9] GenomicRanges_1.46.1           GenomeInfoDb_1.30.0           
[11] IRanges_2.28.0                 S4Vectors_0.32.3              
[13] BiocGenerics_0.40.0            MatrixGenerics_1.6.0          
[15] matrixStats_0.61.0-9001        BiocStyle_2.22.0              
[17] rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0            colorspace_2.0-2           
  [3] ellipsis_0.3.2              scuttle_1.4.0              
  [5] BiocNeighbors_1.12.0        ggrepel_0.9.1              
  [7] bit64_4.0.5                 fansi_1.0.0                
  [9] decontam_1.14.0             splines_4.1.2              
 [11] codetools_0.2-18            sparseMatrixStats_1.6.0    
 [13] cachem_1.0.6                knitr_1.37                 
 [15] scater_1.22.0               jsonlite_1.7.2             
 [17] cluster_2.1.2               graph_1.72.0               
 [19] BiocManager_1.30.16         compiler_4.1.2             
 [21] assertthat_0.2.1            Matrix_1.4-0               
 [23] fastmap_1.1.0               lazyeval_0.2.2             
 [25] BiocSingular_1.10.0         htmltools_0.5.2            
 [27] tools_4.1.2                 rsvd_1.0.5                 
 [29] gtable_0.3.0                glue_1.6.0                 
 [31] GenomeInfoDbData_1.2.7      reshape2_1.4.4             
 [33] dplyr_1.0.7                 Rcpp_1.0.7                 
 [35] jquerylib_0.1.4             vctrs_0.3.8                
 [37] ape_5.6-1                   nlme_3.1-153               
 [39] DECIPHER_2.22.0             DelayedMatrixStats_1.16.0  
 [41] xfun_0.29                   stringr_1.4.0              
 [43] beachmat_2.10.0             lifecycle_1.0.1            
 [45] irlba_2.3.5                 XML_3.99-0.8               
 [47] zlibbioc_1.40.0             MASS_7.3-54                
 [49] scales_1.1.1                parallel_4.1.2             
 [51] yaml_2.2.1                  memoise_2.0.1              
 [53] gridExtra_2.3               ggplot2_3.3.5              
 [55] yulab.utils_0.0.4           sass_0.4.0                 
 [57] stringi_1.7.6               RSQLite_2.2.9              
 [59] ScaledMatrix_1.2.0          permute_0.9-5              
 [61] tidytree_0.3.6              filelock_1.0.2             
 [63] BiocParallel_1.28.3         rlang_0.4.12               
 [65] pkgconfig_2.0.3             bitops_1.0-7               
 [67] evaluate_0.14               lattice_0.20-45            
 [69] purrr_0.3.4                 treeio_1.18.1              
 [71] CodeDepends_0.6.5           bit_4.0.4                  
 [73] tidyselect_1.1.1            plyr_1.8.6                 
 [75] magrittr_2.0.1              bookdown_0.24              
 [77] R6_2.5.1                    generics_0.1.1             
 [79] DelayedArray_0.20.0         DBI_1.1.2                  
 [81] mgcv_1.8-38                 pillar_1.6.4               
 [83] RCurl_1.98-1.5              tibble_3.1.6               
 [85] dir.expiry_1.2.0            crayon_1.4.2               
 [87] utf8_1.2.2                  rmarkdown_2.11             
 [89] viridis_0.6.2               grid_4.1.2                 
 [91] blob_1.2.2                  vegan_2.5-7                
 [93] digest_0.6.29               tidyr_1.1.4                
 [95] munsell_0.5.0               DirichletMultinomial_1.36.0
 [97] beeswarm_0.4.0              viridisLite_0.4.0          
 [99] vipor_0.4.5                 bslib_0.3.1                
```
</div>
