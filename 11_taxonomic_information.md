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

The function `getUniqueTaxa` gives a list of unique taxa for the specified taxonomic rank.


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
tse <- transformSamples(tse, method = "relabundance")
tse <- transformSamples(x = tse, assay_name = "relabundance", method = "clr", 
                        pseudocount = 1, name = "clr_transformation")

head(assay(tse, "clr_transformation"))
```

```
##                                        CL3        CC1        SV1    M31Fcsw
## Class:Thermoprotei               -5.19e-05 -5.192e-05 -5.194e-05 -5.035e-05
## Class:Thermoprotei               -5.19e-05 -5.192e-05 -5.194e-05 -5.035e-05
## Species:Sulfolobusacidocaldarius -5.19e-05 -5.192e-05 -5.194e-05 -5.035e-05
## Class:Sd-NA                      -5.19e-05 -5.192e-05 -5.194e-05 -5.035e-05
## Class:Sd-NA                      -5.19e-05 -5.192e-05 -5.194e-05 -5.035e-05
## Class:Sd-NA                      -5.19e-05 -5.192e-05 -5.194e-05 -5.035e-05
##                                     M11Fcsw   M31Plmr    M11Plmr    F21Plmr
## Class:Thermoprotei               -4.991e-05 -5.06e-05 -5.091e-05 -5.148e-05
## Class:Thermoprotei               -4.991e-05 -5.06e-05 -5.091e-05 -5.148e-05
## Species:Sulfolobusacidocaldarius -4.991e-05 -5.06e-05 -4.860e-05 -5.148e-05
## Class:Sd-NA                      -4.991e-05 -5.06e-05 -5.091e-05 -5.148e-05
## Class:Sd-NA                      -4.991e-05 -5.06e-05 -5.091e-05 -5.148e-05
## Class:Sd-NA                      -4.991e-05 -5.06e-05 -5.091e-05 -5.148e-05
##                                    M31Tong   M11Tong   LMEpi24M   SLEpi20M
## Class:Thermoprotei               -4.89e-05 -5.05e-05 -4.796e-05 -4.910e-05
## Class:Thermoprotei               -4.89e-05 -5.05e-05 -4.796e-05 -4.993e-05
## Species:Sulfolobusacidocaldarius -4.89e-05 -5.05e-05 -4.796e-05 -4.993e-05
## Class:Sd-NA                      -4.89e-05 -5.05e-05 -4.796e-05 -4.993e-05
## Class:Sd-NA                      -4.89e-05 -5.05e-05 -4.796e-05 -4.993e-05
## Class:Sd-NA                      -4.89e-05 -5.05e-05 -4.796e-05 -4.993e-05
##                                      AQC1cm     AQC4cm     AQC7cm        NP2
## Class:Thermoprotei               -2.419e-05 -4.442e-06  2.822e-05 -4.900e-05
## Class:Thermoprotei               -4.731e-05 -4.602e-05 -4.475e-05 -5.091e-05
## Species:Sulfolobusacidocaldarius -4.731e-05 -4.686e-05 -4.828e-05 -5.091e-05
## Class:Sd-NA                      -4.731e-05 -3.753e-05 -3.121e-05 -5.091e-05
## Class:Sd-NA                      -4.731e-05 -4.602e-05 -4.769e-05 -5.091e-05
## Class:Sd-NA                      -4.731e-05 -4.644e-05 -4.651e-05 -5.091e-05
##                                         NP3        NP5    TRRsed1    TRRsed2
## Class:Thermoprotei               -5.133e-05 -5.141e-05 -5.185e-05 -5.116e-05
## Class:Thermoprotei               -5.133e-05 -5.141e-05 -5.185e-05 -5.116e-05
## Species:Sulfolobusacidocaldarius -5.133e-05 -5.141e-05 -5.185e-05 -5.116e-05
## Class:Sd-NA                      -5.133e-05 -5.141e-05 -5.185e-05 -5.116e-05
## Class:Sd-NA                      -5.133e-05 -5.141e-05 -5.185e-05 -5.116e-05
## Class:Sd-NA                      -5.133e-05 -5.141e-05 -5.185e-05 -5.116e-05
##                                     TRRsed3       TS28      TS29      Even1
## Class:Thermoprotei               -5.158e-05 -5.117e-05 -5.01e-05 -5.125e-05
## Class:Thermoprotei               -5.158e-05 -5.117e-05 -5.01e-05 -5.125e-05
## Species:Sulfolobusacidocaldarius -5.158e-05 -5.117e-05 -5.01e-05 -5.125e-05
## Class:Sd-NA                      -5.158e-05 -5.117e-05 -5.01e-05 -5.125e-05
## Class:Sd-NA                      -5.158e-05 -5.117e-05 -5.01e-05 -5.125e-05
## Class:Sd-NA                      -5.158e-05 -5.117e-05 -5.01e-05 -5.125e-05
##                                       Even2      Even3
## Class:Thermoprotei               -5.114e-05 -5.122e-05
## Class:Thermoprotei               -5.114e-05 -5.122e-05
## Species:Sulfolobusacidocaldarius -5.114e-05 -5.122e-05
## Class:Sd-NA                      -5.114e-05 -5.122e-05
## Class:Sd-NA                      -5.114e-05 -5.122e-05
## Class:Sd-NA                      -5.114e-05 -5.122e-05
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
R version 4.2.0 (2022-04-22)
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
 [1] mia_1.3.34                     MultiAssayExperiment_1.22.0   
 [3] TreeSummarizedExperiment_2.1.4 Biostrings_2.64.0             
 [5] XVector_0.36.0                 SingleCellExperiment_1.18.0   
 [7] SummarizedExperiment_1.26.1    Biobase_2.56.0                
 [9] GenomicRanges_1.48.0           GenomeInfoDb_1.32.2           
[11] IRanges_2.30.0                 S4Vectors_0.34.0              
[13] BiocGenerics_0.42.0            MatrixGenerics_1.8.1          
[15] matrixStats_0.62.0-9000        BiocStyle_2.24.0              
[17] rebook_1.6.0                  

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0            colorspace_2.0-3           
 [3] ellipsis_0.3.2              scuttle_1.6.2              
 [5] BiocNeighbors_1.14.0        ggrepel_0.9.1              
 [7] bit64_4.0.5                 fansi_1.0.3                
 [9] decontam_1.16.0             splines_4.2.0              
[11] codetools_0.2-18            sparseMatrixStats_1.8.0    
[13] cachem_1.0.6                knitr_1.39                 
[15] scater_1.24.0               jsonlite_1.8.0             
[17] cluster_2.1.3               graph_1.74.0               
[19] BiocManager_1.30.18         compiler_4.2.0             
[21] assertthat_0.2.1            Matrix_1.4-1               
[23] fastmap_1.1.0               lazyeval_0.2.2             
[25] cli_3.3.0                   BiocSingular_1.12.0        
[27] htmltools_0.5.3             tools_4.2.0                
[29] rsvd_1.0.5                  gtable_0.3.0               
[31] glue_1.6.2                  GenomeInfoDbData_1.2.8     
[33] reshape2_1.4.4              dplyr_1.0.9                
[35] Rcpp_1.0.9                  vctrs_0.4.1                
[37] ape_5.6-2                   nlme_3.1-158               
[39] DECIPHER_2.24.0             DelayedMatrixStats_1.18.0  
[41] xfun_0.31                   stringr_1.4.0              
[43] beachmat_2.12.0             lifecycle_1.0.1            
[45] irlba_2.3.5                 XML_3.99-0.10              
[47] zlibbioc_1.42.0             MASS_7.3-58                
[49] scales_1.2.0                parallel_4.2.0             
[51] yaml_2.3.5                  memoise_2.0.1              
[53] gridExtra_2.3               ggplot2_3.3.6              
[55] yulab.utils_0.0.5           stringi_1.7.8              
[57] RSQLite_2.2.15              ScaledMatrix_1.4.0         
[59] tidytree_0.3.9              permute_0.9-7              
[61] filelock_1.0.2              BiocParallel_1.30.3        
[63] rlang_1.0.4                 pkgconfig_2.0.3            
[65] bitops_1.0-7                evaluate_0.15              
[67] lattice_0.20-45             purrr_0.3.4                
[69] treeio_1.20.1               CodeDepends_0.6.5          
[71] bit_4.0.4                   tidyselect_1.1.2           
[73] plyr_1.8.7                  magrittr_2.0.3             
[75] bookdown_0.27               R6_2.5.1                   
[77] generics_0.1.3              DelayedArray_0.22.0        
[79] DBI_1.1.3                   mgcv_1.8-40                
[81] pillar_1.8.0                RCurl_1.98-1.7             
[83] tibble_3.1.7                dir.expiry_1.4.0           
[85] crayon_1.5.1                utf8_1.2.2                 
[87] rmarkdown_2.14              viridis_0.6.2              
[89] grid_4.2.0                  blob_1.2.3                 
[91] vegan_2.6-2                 digest_0.6.29              
[93] tidyr_1.2.0                 munsell_0.5.0              
[95] DirichletMultinomial_1.38.0 beeswarm_0.4.0             
[97] viridisLite_0.4.0           vipor_0.4.5                
```
</div>
