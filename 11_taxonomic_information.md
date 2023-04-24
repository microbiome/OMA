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

The dada2 package [@Callahan2016dada2] implements the `assignTaxonomy`
function, which takes as input the ASV sequences associated with each
row of data and a training dataset. For more information visit the
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

By default the return value of `getTaxonomyLabels` contains only
unique elements by passing it through `make.unique`. This step can be
omitted by setting `make_unique = FALSE`.


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

The function `getUniqueTaxa` gives a list of unique taxa for the
specified taxonomic rank.


```r
head(getUniqueTaxa(tse, rank = "Phylum"))
```

```
## [1] "Crenarchaeota"  "Euryarchaeota"  "Actinobacteria" "Spirochaetes"  
## [5] "MVP-15"         "Proteobacteria"
```


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
## metadata(1): agglomerated_by_rank
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
any taxonomic level listed in `mia::taxonomyRanks(tse)`.   


## Data transformation

Data transformations are common in microbiome analysis. Examples
include the logarithmic transformation, calculation of relative
abundances (percentages), and compositionality-aware transformations
such as the centered log-ratio transformation (clr).

In mia package, transformations are applied to abundance data. The transformed 
abundance table is stored back to 'assays'. mia includes transformation 
function ('transformCounts()') which applies sample-wise or column-wise transformation when MARGIN = 'samples', feature-wise or row-wise transformation when MARGIN = 'features'.

For a complete list of available transformations and parameters, see function 
[help](https://microbiome.github.io/mia/reference/transformCounts.html).


```r
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformCounts(tse, assay_name = "pseudo", method = "relabundance")
tse <- transformCounts(x = tse, assay_name = "relabundance", method = "clr", 
                        pseudocount = 1, name = "clr_transformation")

head(assay(tse, "clr_transformation"))
```

```
##                                         CL3        CC1        SV1    M31Fcsw
## Class:Thermoprotei               -5.078e-05 -5.105e-05 -5.055e-05 -4.975e-05
## Class:Thermoprotei               -5.078e-05 -5.105e-05 -5.055e-05 -4.975e-05
## Species:Sulfolobusacidocaldarius -5.078e-05 -5.105e-05 -5.055e-05 -4.975e-05
## Class:Sd-NA                      -5.078e-05 -5.105e-05 -5.055e-05 -4.975e-05
## Class:Sd-NA                      -5.078e-05 -5.105e-05 -5.055e-05 -4.975e-05
## Class:Sd-NA                      -5.078e-05 -5.105e-05 -5.055e-05 -4.975e-05
##                                     M11Fcsw    M31Plmr    M11Plmr    F21Plmr
## Class:Thermoprotei               -4.947e-05 -4.931e-05 -4.879e-05 -4.671e-05
## Class:Thermoprotei               -4.947e-05 -4.931e-05 -4.879e-05 -4.671e-05
## Species:Sulfolobusacidocaldarius -4.947e-05 -4.931e-05 -4.658e-05 -4.671e-05
## Class:Sd-NA                      -4.947e-05 -4.931e-05 -4.879e-05 -4.671e-05
## Class:Sd-NA                      -4.947e-05 -4.931e-05 -4.879e-05 -4.671e-05
## Class:Sd-NA                      -4.947e-05 -4.931e-05 -4.879e-05 -4.671e-05
##                                     M31Tong    M11Tong   LMEpi24M   SLEpi20M
## Class:Thermoprotei               -4.846e-05 -4.257e-05 -4.756e-05 -4.837e-05
## Class:Thermoprotei               -4.846e-05 -4.257e-05 -4.756e-05 -4.918e-05
## Species:Sulfolobusacidocaldarius -4.846e-05 -4.257e-05 -4.756e-05 -4.918e-05
## Class:Sd-NA                      -4.846e-05 -4.257e-05 -4.756e-05 -4.918e-05
## Class:Sd-NA                      -4.846e-05 -4.257e-05 -4.756e-05 -4.918e-05
## Class:Sd-NA                      -4.846e-05 -4.257e-05 -4.756e-05 -4.918e-05
##                                      AQC1cm     AQC4cm     AQC7cm        NP2
## Class:Thermoprotei               -2.385e-05 -4.438e-06  2.787e-05 -4.731e-05
## Class:Thermoprotei               -4.660e-05 -4.568e-05 -4.428e-05 -4.915e-05
## Species:Sulfolobusacidocaldarius -4.660e-05 -4.652e-05 -4.777e-05 -4.915e-05
## Class:Sd-NA                      -4.660e-05 -3.726e-05 -3.090e-05 -4.915e-05
## Class:Sd-NA                      -4.660e-05 -4.568e-05 -4.719e-05 -4.915e-05
## Class:Sd-NA                      -4.660e-05 -4.610e-05 -4.603e-05 -4.915e-05
##                                         NP3        NP5    TRRsed1    TRRsed2
## Class:Thermoprotei               -5.068e-05 -5.083e-05 -3.909e-05 -4.927e-05
## Class:Thermoprotei               -5.068e-05 -5.083e-05 -3.909e-05 -4.927e-05
## Species:Sulfolobusacidocaldarius -5.068e-05 -5.083e-05 -3.909e-05 -4.927e-05
## Class:Sd-NA                      -5.068e-05 -5.083e-05 -3.909e-05 -4.927e-05
## Class:Sd-NA                      -5.068e-05 -5.083e-05 -3.909e-05 -4.927e-05
## Class:Sd-NA                      -5.068e-05 -5.083e-05 -3.909e-05 -4.927e-05
##                                     TRRsed3       TS28       TS29      Even1
## Class:Thermoprotei               -4.829e-05 -5.016e-05 -4.934e-05 -5.046e-05
## Class:Thermoprotei               -4.829e-05 -5.016e-05 -4.934e-05 -5.046e-05
## Species:Sulfolobusacidocaldarius -4.829e-05 -5.016e-05 -4.934e-05 -5.046e-05
## Class:Sd-NA                      -4.829e-05 -5.016e-05 -4.934e-05 -5.046e-05
## Class:Sd-NA                      -4.829e-05 -5.016e-05 -4.934e-05 -5.046e-05
## Class:Sd-NA                      -4.829e-05 -5.016e-05 -4.934e-05 -5.046e-05
##                                       Even2      Even3
## Class:Thermoprotei               -5.017e-05 -5.034e-05
## Class:Thermoprotei               -5.017e-05 -5.034e-05
## Species:Sulfolobusacidocaldarius -5.017e-05 -5.034e-05
## Class:Sd-NA                      -5.017e-05 -5.034e-05
## Class:Sd-NA                      -5.017e-05 -5.034e-05
## Class:Sd-NA                      -5.017e-05 -5.034e-05
```

-   In 'pa' transformation, abundance table is converted to present/absent table.


```r
tse <- transformCounts(tse, method = "pa")

head(assay(tse, "pa"))
```

```
##                                  CL3 CC1 SV1 M31Fcsw M11Fcsw M31Plmr M11Plmr
## Class:Thermoprotei                 0   0   0       0       0       0       0
## Class:Thermoprotei                 0   0   0       0       0       0       0
## Species:Sulfolobusacidocaldarius   0   0   0       0       0       0       1
## Class:Sd-NA                        0   0   0       0       0       0       0
## Class:Sd-NA                        0   0   0       0       0       0       0
## Class:Sd-NA                        0   0   0       0       0       0       0
##                                  F21Plmr M31Tong M11Tong LMEpi24M SLEpi20M
## Class:Thermoprotei                     0       0       0        0        1
## Class:Thermoprotei                     0       0       0        0        0
## Species:Sulfolobusacidocaldarius       0       0       0        0        0
## Class:Sd-NA                            0       0       0        0        0
## Class:Sd-NA                            0       0       0        0        0
## Class:Sd-NA                            0       0       0        0        0
##                                  AQC1cm AQC4cm AQC7cm NP2 NP3 NP5 TRRsed1
## Class:Thermoprotei                    1      1      1   1   0   0       0
## Class:Thermoprotei                    0      1      1   0   0   0       0
## Species:Sulfolobusacidocaldarius      0      0      0   0   0   0       0
## Class:Sd-NA                           0      1      1   0   0   0       0
## Class:Sd-NA                           0      1      1   0   0   0       0
## Class:Sd-NA                           0      1      1   0   0   0       0
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
## List of length 5
## names(5): counts relabundance pseudo clr_transformation pa
```

## Pick specific  

Retrieving of specific elements that are required for specific analysis. For
instance, extracting abundances for a specific taxa in all samples or all taxa 
in one sample.  

### Abundances of all taxa in specific sample 

```r
taxa.abund.cc1 <- getAbundanceSample(tse, 
                                     sample_id = "CC1",
                                     assay_name = "counts")
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
                                      assay_name = "counts")
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
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] mia_1.7.11                     MultiAssayExperiment_1.24.0   
 [3] TreeSummarizedExperiment_2.1.4 Biostrings_2.66.0             
 [5] XVector_0.38.0                 SingleCellExperiment_1.20.1   
 [7] SummarizedExperiment_1.28.0    Biobase_2.58.0                
 [9] GenomicRanges_1.50.2           GenomeInfoDb_1.34.9           
[11] IRanges_2.32.0                 S4Vectors_0.36.2              
[13] BiocGenerics_0.44.0            MatrixGenerics_1.10.0         
[15] matrixStats_0.63.0-9003        BiocStyle_2.24.0              
[17] rebook_1.6.0                  

loaded via a namespace (and not attached):
 [1] nlme_3.1-162                bitops_1.0-7               
 [3] DirichletMultinomial_1.40.0 bit64_4.0.5                
 [5] filelock_1.0.2              tools_4.2.1                
 [7] vegan_2.6-4                 utf8_1.2.3                 
 [9] R6_2.5.1                    irlba_2.3.5.1              
[11] vipor_0.4.5                 mgcv_1.8-42                
[13] DBI_1.1.3                   lazyeval_0.2.2             
[15] colorspace_2.1-0            permute_0.9-7              
[17] withr_2.5.0                 tidyselect_1.2.0           
[19] gridExtra_2.3               bit_4.0.5                  
[21] compiler_4.2.1              graph_1.74.0               
[23] cli_3.6.1                   BiocNeighbors_1.16.0       
[25] DelayedArray_0.24.0         bookdown_0.33              
[27] scales_1.2.1                stringr_1.5.0              
[29] digest_0.6.31               yulab.utils_0.0.6          
[31] rmarkdown_2.21              scater_1.26.1              
[33] pkgconfig_2.0.3             htmltools_0.5.5            
[35] decontam_1.18.0             sparseMatrixStats_1.10.0   
[37] fastmap_1.1.1               rlang_1.1.0                
[39] RSQLite_2.3.1               DelayedMatrixStats_1.20.0  
[41] generics_0.1.3              jsonlite_1.8.4             
[43] BiocParallel_1.32.6         dplyr_1.1.2                
[45] RCurl_1.98-1.12             magrittr_2.0.3             
[47] BiocSingular_1.14.0         GenomeInfoDbData_1.2.9     
[49] scuttle_1.8.4               Matrix_1.5-4               
[51] Rcpp_1.0.10                 ggbeeswarm_0.7.1           
[53] munsell_0.5.0               fansi_1.0.4                
[55] DECIPHER_2.26.0             ape_5.7-1                  
[57] viridis_0.6.2               lifecycle_1.0.3            
[59] stringi_1.7.12              yaml_2.3.7                 
[61] MASS_7.3-58.3               zlibbioc_1.44.0            
[63] plyr_1.8.8                  blob_1.2.4                 
[65] grid_4.2.1                  parallel_4.2.1             
[67] ggrepel_0.9.3               crayon_1.5.2               
[69] dir.expiry_1.4.0            lattice_0.21-8             
[71] splines_4.2.1               beachmat_2.14.2            
[73] CodeDepends_0.6.5           knitr_1.42                 
[75] pillar_1.9.0                reshape2_1.4.4             
[77] codetools_0.2-19            ScaledMatrix_1.6.0         
[79] XML_3.99-0.14               glue_1.6.2                 
[81] evaluate_0.20               BiocManager_1.30.20        
[83] vctrs_0.6.2                 treeio_1.22.0              
[85] gtable_0.3.3                purrr_1.0.1                
[87] tidyr_1.3.0                 cachem_1.0.7               
[89] ggplot2_3.4.2               xfun_0.39                  
[91] rsvd_1.0.5                  tidytree_0.4.2             
[93] viridisLite_0.4.1           tibble_3.2.1               
[95] memoise_2.0.1               beeswarm_0.4.0             
[97] cluster_2.1.4              
```
</div>
