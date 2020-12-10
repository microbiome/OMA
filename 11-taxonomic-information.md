# Taxonomic information {#taxonomic-information}

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
data("GlobalPatterns")
se <- GlobalPatterns 
```

Taxonomic information are a key part of analyzing microbiome data and without
it, any type of data analysis probably will not make much sense. However,
the degree of detail of taxonomic information differs depending on the dataset
and annotation data used.

Therefore, the mia package expects a loose assembly of taxonomic information
and assumes certain key aspects:

* Taxonomic information is given as character vectors or factors in the 
`rowData` of an `SummarizedExperiment` object.
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
[dada2 website](https://benjjneb.github.io/dada2/assign.html).

### DECIPHER

The DECIPHER package [@R-DECIPHER] implements the `IDTAXA` algorithm to assign
either taxonomic information or function information. For `mia`
only the first option is of interest for now and more information can be
found on the [DECIPHER website](http://www2.decipher.codes/Classification.html)

## Functions to access taxonomic information

`checkTaxonomy` checks whether the taxonomic information is usable for `mia`


```r
checkTaxonomy(se)
```

```
## [1] TRUE
```

Since the `rowData` can contain other data, `taxonomyRanks` will return the 
columns `mia` assumes to contain the taxonomic information.


```r
taxonomyRanks(se)
```

```
## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
```

This can then be used to subset the `rowData` to columns needed.


```r
rowData(se)[,taxonomyRanks(se)]
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
all(!taxonomyRankEmpty(se, rank = "Kingdom"))
```

```
## [1] TRUE
```

```r
table(taxonomyRankEmpty(se, rank = "Genus"))
```

```
## 
## FALSE  TRUE 
##  8008 11208
```

```r
table(taxonomyRankEmpty(se, rank = "Species"))
```

```
## 
## FALSE  TRUE 
##  1413 17803
```

`getTaxonomyLabels` is a multi-purpose function, which turns taxonomic
information into a character vector of `length(x)`


```r
head(getTaxonomyLabels(se))
```

```
## [1] "Class:Thermoprotei"               "Class:Thermoprotei_1"            
## [3] "Species:Sulfolobusacidocaldarius" "Class:Sd-NA"                     
## [5] "Class:Sd-NA_1"                    "Class:Sd-NA_2"
```

By default this will used the lowest non-empty information to construct a
string with the following scheme `level:value`. If all levels are the same
this part is omited, but can be added by setting `with_rank = TRUE`


```r
phylum <- !is.na(rowData(se)$Phylum) & 
    vapply(data.frame(apply(rowData(se)[,taxonomyRanks(se)[3:7]],1L,is.na)),all,logical(1))
head(getTaxonomyLabels(se[phylum,]))
```

```
## [1] "Crenarchaeota"    "Crenarchaeota_1"  "Crenarchaeota_2"  "Actinobacteria"  
## [5] "Actinobacteria_1" "Spirochaetes"
```

```r
head(getTaxonomyLabels(se[phylum,], with_rank = TRUE))
```

```
## [1] "Phylum:Crenarchaeota"    "Phylum:Crenarchaeota_1" 
## [3] "Phylum:Crenarchaeota_2"  "Phylum:Actinobacteria"  
## [5] "Phylum:Actinobacteria_1" "Phylum:Spirochaetes"
```

By default the return value of `getTaxonomyLabels` contains only unique elements
by passing it through `make.unique`. This step can be omited by setting 
`make_unique = FALSE`


```r
head(getTaxonomyLabels(se[phylum,], with_rank = TRUE, make_unique = FALSE))
```

```
## [1] "Phylum:Crenarchaeota"  "Phylum:Crenarchaeota"  "Phylum:Crenarchaeota" 
## [4] "Phylum:Actinobacteria" "Phylum:Actinobacteria" "Phylum:Spirochaetes"
```

To apply the loop resolving function `resolveLoop` from the
`TreeSummarizedExperiment` package [@R-TreeSummarizedExperiment] within
`getTaxonomyLabels`, set `resolve_loops = TRUE`.

### Generate a taxonomic tree on the fly

To create a taxonomic tree `taxonomyTree` used the information and returns a
`phylo` object. Duplicate information from the `rowData` are removed.


```r
taxonomyTree(se)
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
se <- addTaxonomyTree(se)
se
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
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: a phylo (1645 leaves)
## colLinks: NULL
## colTree: NULL
```

The implementation is based on the the `toTree` function from the
`TreeSummarizedExperiment` package [@R-TreeSummarizedExperiment].

## Data agglomeration {#data-agglomeration}

One of the main applications of taxonomic information in regards to count data
is to agglomerate count data on taxonomic levels and track the influence of 
changing conditions through these levels. For this `mia` contains the
`agglomerateByRank` function. The ideal location to store the agglomerated data
is as an alternative experiment.


```r
se <- relAbundanceCounts(se)
altExp(se, "Family") <- agglomerateByRank(se, rank = "Family",
                                          agglomerateTree = TRUE)
altExp(se, "Family")
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
## altExpNames(0):
## rowLinks: a LinkDataFrame (603 rows)
## rowTree: a phylo (876 leaves)
## colLinks: NULL
## colTree: NULL
```

If multiple assays (counts and relabundance) exists, both will be agglomerated.


```r
assayNames(se)
```

```
## [1] "counts"       "relabundance"
```

```r
assayNames(altExp(se, "Family"))
```

```
## [1] "counts"       "relabundance"
```


```r
assay(altExp(se, "Family"), "relabundance")[1:5,1:7]
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
assay(altExp(se, "Family"), "counts")[1:5,1:7]
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

## Pick specific  

Retrieving of specific elements are required for specific analysis. For
instance, extracting abundances for a specific taxa in all samples or all taxa 
in one sample.  

### Abundances of all taxa in specific sample 

```r
taxa.abund.cc1 <- getAbundanceSample(se, 
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
taxa.abundances <- getAbundanceFeature(se, 
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
 [1] mia_0.98.15                      MicrobiomeExperiment_0.99.0.9014
 [3] Biostrings_2.58.0                XVector_0.30.0                  
 [5] TreeSummarizedExperiment_1.6.2   SingleCellExperiment_1.12.0     
 [7] SummarizedExperiment_1.20.0      Biobase_2.50.0                  
 [9] GenomicRanges_1.42.0             GenomeInfoDb_1.26.2             
[11] IRanges_2.24.0                   S4Vectors_0.28.1                
[13] BiocGenerics_0.36.0              MatrixGenerics_1.2.0            
[15] matrixStats_0.57.0               BiocStyle_2.18.1                
[17] rebook_1.0.0                     BiocManager_1.30.10             

loaded via a namespace (and not attached):
 [1] viridis_0.5.1             tidyr_1.1.2              
 [3] BiocSingular_1.6.0        viridisLite_0.3.0        
 [5] DelayedMatrixStats_1.12.1 scuttle_1.0.3            
 [7] vipor_0.4.5               GenomeInfoDbData_1.2.4   
 [9] yaml_2.2.1                pillar_1.4.7             
[11] lattice_0.20-41           glue_1.4.2               
[13] beachmat_2.6.2            digest_0.6.27            
[15] colorspace_2.0-0          htmltools_0.5.0          
[17] Matrix_1.2-18             XML_3.99-0.5             
[19] pkgconfig_2.0.3           bookdown_0.21            
[21] zlibbioc_1.36.0           purrr_0.3.4              
[23] scales_1.1.1              processx_3.4.5           
[25] BiocParallel_1.24.1       tibble_3.0.4             
[27] generics_0.1.0            ggplot2_3.3.2            
[29] ellipsis_0.3.1            magrittr_2.0.1           
[31] crayon_1.3.4              CodeDepends_0.6.5        
[33] evaluate_0.14             ps_1.5.0                 
[35] nlme_3.1-150              beeswarm_0.2.3           
[37] graph_1.68.0              tools_4.0.3              
[39] scater_1.18.3             lifecycle_0.2.0          
[41] stringr_1.4.0             munsell_0.5.0            
[43] DelayedArray_0.16.0       irlba_2.3.3              
[45] callr_3.5.1               compiler_4.0.3           
[47] rsvd_1.0.3                rlang_0.4.9              
[49] grid_4.0.3                RCurl_1.98-1.2           
[51] BiocNeighbors_1.8.2       bitops_1.0-6             
[53] rmarkdown_2.5             gtable_0.3.0             
[55] codetools_0.2-18          R6_2.5.0                 
[57] gridExtra_2.3             knitr_1.30               
[59] dplyr_1.0.2               ape_5.4-1                
[61] stringi_1.5.3             ggbeeswarm_0.6.0         
[63] Rcpp_1.0.5                vctrs_0.3.5              
[65] tidyselect_1.1.0          xfun_0.19                
[67] sparseMatrixStats_1.2.0  
```
</div>
