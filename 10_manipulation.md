# (PART) Focus Topics {-}


# Data Manipulation {#datamanipulation}

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


## Tidying and subsetting

### Tidy data

For several custom analysis and visualization packages, such as those from
`tidyverse`, the `SE` data can be converted to a long data.frame format with 
`meltAssay`.    




```r
library(mia)
data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
tse <- transformAssay(tse, MARGIN = "samples", method="relabundance")
molten_tse <- mia::meltAssay(tse,
                        add_row_data = TRUE,
                        add_col_data = TRUE,
                        assay.type = "relabundance")
molten_tse
```

```
## # A tibble: 499,616 x 17
##    FeatureID SampleID relabundance Kingdom Phylum       Class Order Family Genus
##    <fct>     <fct>           <dbl> <chr>   <chr>        <chr> <chr> <chr>  <chr>
##  1 549322    CL3                 0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  2 549322    CC1                 0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  3 549322    SV1                 0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  4 549322    M31Fcsw             0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  5 549322    M11Fcsw             0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  6 549322    M31Plmr             0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  7 549322    M11Plmr             0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  8 549322    F21Plmr             0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
##  9 549322    M31Tong             0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
## 10 549322    M11Tong             0 Archaea Crenarchaeo~ Ther~ <NA>  <NA>   <NA> 
## # i 499,606 more rows
## # i 8 more variables: Species <chr>, X.SampleID <fct>, Primer <fct>,
## #   Final_Barcode <fct>, Barcode_truncated_plus_T <fct>,
## #   Barcode_full_length <fct>, SampleType <fct>, Description <fct>
```

### Subsetting

**Subsetting** data helps to draw the focus of analysis on particular
  sets of samples and / or features. When dealing with large datasets, 
  the subset of interest can be extracted and investigated
  separately. This might improve performance and reduce the
  computational load.

Load:

* mia
* dplyr
* knitr
* data `GlobalPatterns`



Let us store `GlobalPatterns` into `tse` and check its original number of features (rows) and samples (columns). **Note**: when subsetting by sample, expect the number of columns to decrease; when subsetting by feature, expect the number of rows to decrease.


```r
# Store data into se and check dimensions
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
# Show dimensions (features x samples)
dim(tse) 
```

```
## [1] 19216    26
```

#### Subset by sample (column-wise)

For the sake of demonstration, here we will extract a subset containing only the samples of human origin (feces, skin or tongue), stored as `SampleType` within `colData(tse)` and also in `tse`.

First, we would like to see all the possible values that `SampleType` can take on and how frequent those are: 


```r
# Inspect possible values for SampleType
unique(tse$SampleType)
```

```
## [1] Soil               Feces              Skin               Tongue            
## [5] Freshwater         Freshwater (creek) Ocean              Sediment (estuary)
## [9] Mock              
## 9 Levels: Feces Freshwater Freshwater (creek) Mock ... Tongue
```

```r
# Show the frequency of each value
tse$SampleType %>% table()
```
\begin{table}
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{l|r}
\hline
. & Freq\\
\hline
Feces & 4\\
\hline
Freshwater & 2\\
\hline
Freshwater (creek) & 3\\
\hline
Mock & 3\\
\hline
Ocean & 3\\
\hline
Sediment (estuary) & 3\\
\hline
Skin & 3\\
\hline
Soil & 3\\
\hline
Tongue & 2\\
\hline
\end{tabular}}
\end{table}

**Note**: after subsetting, expect the number of columns to equal the
  sum of the frequencies of the samples that you are interested
  in. For instance, `ncols = Feces + Skin + Tongue = 4 + 3 + 2 = 9`.

Next, we _logical index_ across the columns of `tse` (make sure to
leave the first index empty to select all rows) and filter for the
samples of human origin. For this, we use the information on the
samples from the meta data `colData(tse)`.


```r
# Subset by sample
tse_subset_by_sample <- tse[ , tse$SampleType %in% c("Feces", "Skin", "Tongue")]

# Show dimensions
dim(tse_subset_by_sample)
```

```
## [1] 19216     9
```

As a sanity check, the new object `tse_subset_by_sample` should have
the original number of features (rows) and a number of samples
(columns) equal to the sum of the samples of interest (in this case
9).

Several characteristics can be used to subset by sample:

* origin
* sampling time
* sequencing method
* DNA / RNA barcode
* cohort

#### Subset by feature (row-wise)

Similarly, here we will extract a subset containing only the features
that belong to the phyla Actinobacteria and Chlamydiae, stored as
`Phylum` within `rowData(tse)`. However, subsetting by feature implies
a few more obstacles, such as the presence of `NA` elements and the
possible need for agglomeration.

As previously, we would first like to see all the possible values that
`Phylum` can take on and how frequent those are:
  

```r
# Inspect possible values for phylum
unique(rowData(tse)$Phylum)
```

```
##  [1] "Crenarchaeota"    "Euryarchaeota"    "Actinobacteria"   "Spirochaetes"    
##  [5] "MVP-15"           "Proteobacteria"   "SBR1093"          "Fusobacteria"    
##  [9] "Tenericutes"      "ZB3"              "Cyanobacteria"    "GOUTA4"          
## [13] "TG3"              "Chlorobi"         "Bacteroidetes"    "Caldithrix"      
## [17] "KSB1"             "SAR406"           "LCP-89"           "Thermi"          
## [21] "Gemmatimonadetes" "Fibrobacteres"    "GN06"             "AC1"             
## [25] "TM6"              "OP8"              "Elusimicrobia"    "NC10"            
## [29] "SPAM"             NA                 "Acidobacteria"    "CCM11b"          
## [33] "Nitrospirae"      "NKB19"            "BRC1"             "Hyd24-12"        
## [37] "WS3"              "PAUC34f"          "GN04"             "GN12"            
## [41] "Verrucomicrobia"  "Lentisphaerae"    "LD1"              "Chlamydiae"      
## [45] "OP3"              "Planctomycetes"   "Firmicutes"       "OP9"             
## [49] "WPS-2"            "Armatimonadetes"  "SC3"              "TM7"             
## [53] "GN02"             "SM2F11"           "ABY1_OD1"         "ZB2"             
## [57] "OP11"             "Chloroflexi"      "SC4"              "WS1"             
## [61] "GAL15"            "AD3"              "WS2"              "Caldiserica"     
## [65] "Thermotogae"      "Synergistetes"    "SR1"
```

```r
# Show the frequency of each value
rowData(tse)$Phylum %>% table()
```
\begin{table}
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{l|r}
\hline
. & Freq\\
\hline
ABY1\_OD1 & 7\\
\hline
AC1 & 1\\
\hline
Acidobacteria & 1021\\
\hline
Actinobacteria & 1631\\
\hline
AD3 & 9\\
\hline
Armatimonadetes & 61\\
\hline
Bacteroidetes & 2382\\
\hline
BRC1 & 13\\
\hline
Caldiserica & 3\\
\hline
Caldithrix & 10\\
\hline
CCM11b & 2\\
\hline
Chlamydiae & 21\\
\hline
Chlorobi & 64\\
\hline
Chloroflexi & 437\\
\hline
Crenarchaeota & 106\\
\hline
Cyanobacteria & 393\\
\hline
Elusimicrobia & 31\\
\hline
Euryarchaeota & 102\\
\hline
Fibrobacteres & 7\\
\hline
Firmicutes & 4356\\
\hline
Fusobacteria & 37\\
\hline
GAL15 & 2\\
\hline
Gemmatimonadetes & 191\\
\hline
GN02 & 8\\
\hline
GN04 & 7\\
\hline
GN06 & 2\\
\hline
GN12 & 1\\
\hline
GOUTA4 & 11\\
\hline
Hyd24-12 & 4\\
\hline
KSB1 & 6\\
\hline
LCP-89 & 2\\
\hline
LD1 & 2\\
\hline
Lentisphaerae & 21\\
\hline
MVP-15 & 5\\
\hline
NC10 & 9\\
\hline
Nitrospirae & 74\\
\hline
NKB19 & 16\\
\hline
OP11 & 6\\
\hline
OP3 & 30\\
\hline
OP8 & 27\\
\hline
OP9 & 4\\
\hline
PAUC34f & 3\\
\hline
Planctomycetes & 638\\
\hline
Proteobacteria & 6416\\
\hline
SAR406 & 21\\
\hline
SBR1093 & 9\\
\hline
SC3 & 8\\
\hline
SC4 & 8\\
\hline
SM2F11 & 5\\
\hline
SPAM & 22\\
\hline
Spirochaetes & 124\\
\hline
SR1 & 5\\
\hline
Synergistetes & 7\\
\hline
Tenericutes & 143\\
\hline
TG3 & 5\\
\hline
Thermi & 46\\
\hline
Thermotogae & 1\\
\hline
TM6 & 27\\
\hline
TM7 & 32\\
\hline
Verrucomicrobia & 470\\
\hline
WPS-2 & 20\\
\hline
WS1 & 5\\
\hline
WS2 & 2\\
\hline
WS3 & 70\\
\hline
ZB2 & 2\\
\hline
ZB3 & 2\\
\hline
\end{tabular}}
\end{table}

**Note**: after subsetting, expect the number of columns to equal the
  sum of the frequencies of the feature(s) that you are interested
  in. For instance, `nrows = Actinobacteria + Chlamydiae = 1631 + 21 =
  1652`.

Depending on your research question, you might or might not need to
agglomerate the data in the first place: if you want to find the
abundance of each and every feature that belongs to Actinobacteria and
Chlamydiae, agglomeration is not needed; if you want to find the total
abundance of all features that belong to Actinobacteria or
Chlamydiae, agglomeration is recommended.

##### Non-agglomerated data

Next, we _logical index_ across the rows of `tse` (make sure to leave
the second index empty to select all columns) and filter for the
features that fall in either Actinobacteria or Chlamydiae group. For this,
we use the information on the samples from the metadata
`rowData(tse)`.

The first term with the `%in%` operator includes all the features
of interest, whereas the second term after the AND operator `&`
filters out all features that have an `NA` in place of the phylum variable.


```r
# Subset by feature
tse_subset_by_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), ]

# Show dimensions
dim(tse_subset_by_feature)
```

```
## [1] 1652   26
```

As a sanity check, the new object, `tse_subset_by_feature`, should have the original number of samples (columns) and a number of features (rows) equal to the sum of the features of interest (in this case, 1652).

##### Agglomerated data

When total abundances of certain phyla are of relevance, the data is initially agglomerated by Phylum. Then, similar steps as in the case of non-agglomerated data are followed.


```r
# Agglomerate by phylum
tse_phylum <- tse %>% mergeFeaturesByRank(rank = "Phylum")

# Subset by feature and remove NAs
tse_phylum_subset_by_feature <- tse_phylum[rowData(tse_phylum)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse_phylum)$Phylum), ]

# Show dimensions
dim(tse_phylum_subset_by_feature)
```

```
## [1]  2 26
```

**Note**: as data was agglomerated, the number of rows should equal the
  number of phyla used to index (in this case, just 2).

Alternatively:


```r
# Store features of interest into phyla
phyla <- c("Phylum:Actinobacteria", "Phylum:Chlamydiae")
# subset by feature
tse_phylum_subset_by_feature <- tse_phylum[phyla, ]
# Show dimensions
dim(tse_subset_by_feature)
```

```
## [1] 1652   26
```

The code above returns the non-agglomerated version of the data.

Fewer characteristics can be used to subset by feature:

* Taxonomic rank
* Meta-taxonomic group

For subsetting by kingdom, agglomeration does not apply, whereas for
the other ranks it can be applied if necessary.

#### Subset by sample and feature

Finally, we can subset data by sample and feature at once. The
resulting subset contains all the samples of human origin and all the
features of phyla Actinobacteria or Chlamydiae.


```r
# Subset by sample and feature and remove NAs
tse_subset_by_sample_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), tse$SampleType %in% c("Feces", "Skin", "Tongue")]

# Show dimensions
dim(tse_subset_by_sample_feature)
```

```
## [1] 1652    9
```

**Note**: the dimensions of `tse_subset_by_sample_feature` agree with
  those of the previous subsets (9 columns filtered by sample and 1652
  rows filtered by feature).

If a study was to consider and quantify the presence of Actinobacteria
as well as Chlamydiae in different sites of the human body,
`tse_subset_by_sample_feature` might be a suitable subset to start
with.

#### Remove empty columns and rows

Sometimes data might contain, e.g., features that are not present in any of the  samples.
This can occur, for example, after the data subsetting. In certain analyses, we might want to
remove those instances.


```r
# Agglomerate data at Genus level 
tse_genus <- mergeFeaturesByRank(tse, rank = "Genus")
# List bacteria that we want to include
genera <- c("Class:Thermoprotei", "Genus:Sulfolobus", "Genus:Sediminicola")
# Subset data
tse_genus_sub <- tse_genus[genera, ]

tse_genus_sub
```

```
## class: TreeSummarizedExperiment 
## dim: 3 26 
## metadata(1): agglomerated_by_rank
## assays(1): counts
## rownames(3): Class:Thermoprotei Genus:Sulfolobus Genus:Sediminicola
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (3 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```


```r
# List total counts of each sample
colSums(assay(tse_genus_sub, "counts"))
```

```
##      CL3      CC1      SV1  M31Fcsw  M11Fcsw  M31Plmr  M11Plmr  F21Plmr 
##        1        0        0        1        1        0        4        1 
##  M31Tong  M11Tong LMEpi24M SLEpi20M   AQC1cm   AQC4cm   AQC7cm      NP2 
##        7        3        0        2       64      105      136      222 
##      NP3      NP5  TRRsed1  TRRsed2  TRRsed3     TS28     TS29    Even1 
##     6433     1154        2        2        2        0        0        0 
##    Even2    Even3 
##        2        0
```

Now we can see that certain samples do not include any bacteria. We can remove those.


```r
# Remove samples that do not contain any bacteria
tse_genus_sub <- tse_genus_sub[ , colSums(assay(tse_genus_sub, "counts")) != 0 ]
tse_genus_sub
```

```
## class: TreeSummarizedExperiment 
## dim: 3 18 
## metadata(1): agglomerated_by_rank
## assays(1): counts
## rownames(3): Class:Thermoprotei Genus:Sulfolobus Genus:Sediminicola
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(18): CL3 M31Fcsw ... TRRsed3 Even2
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (3 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

The same action can also be applied to the features.


```r
# Take only those samples that are collected from feces, skin, or tongue
tse_genus_sub <- tse_genus[ , tse_genus$SampleType %in% c("Feces", "Skin", "Tongue")]

tse_genus_sub
```

```
## class: TreeSummarizedExperiment 
## dim: 1516 9 
## metadata(1): agglomerated_by_rank
## assays(1): counts
## rownames(1516): Class:Thermoprotei Genus:Sulfolobus ...
##   Genus:Coprothermobacter Phylum:SR1
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(9): M31Fcsw M11Fcsw ... TS28 TS29
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (1516 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```


```r
# What is the number of bacteria that are not present?
sum(rowSums(assay(tse_genus_sub, "counts")) == 0)
```

```
## [1] 435
```

We can see that there are bacteria that are not present in these samples we chose.
We can remove those bacteria from the data. 


```r
# Take only those bacteria that are present
tse_genus_sub <- tse_genus_sub[rowSums(assay(tse_genus_sub, "counts")) > 0, ]

tse_genus_sub
```

```
## class: TreeSummarizedExperiment 
## dim: 1081 9 
## metadata(1): agglomerated_by_rank
## assays(1): counts
## rownames(1081): Genus:Sulfolobus Order:NRP-J ...
##   Genus:Coprothermobacter Phylum:SR1
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(9): M31Fcsw M11Fcsw ... TS28 TS29
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (1081 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

### Splitting

You can split the data based on variables by using the functions `splitByRanks` 
and `splitOn`.

`splitByRanks` splits the data based on taxonomic ranks. Since the elements of the output list
share columns, they can be stored into `altExp`. 


```r
altExps(tse) <- splitByRanks(tse)
altExps(tse)
```

```
## List of length 7
## names(7): Kingdom Phylum Class Order Family Genus Species
```

If you want to split the data based on another variable than taxonomic rank, use 
`splitOn`. It works for row-wise and column-wise splitting.


```r
splitOn(tse, "SampleType")
```

```
## List of length 9
## names(9): Soil Feces Skin Tongue ... Ocean Sediment (estuary) Mock
```

## Add or modify data

The information contained by the `colData` of a `TreeSE` can be modified by
accessing the desired variables.


```r
# modify the Description entries
colData(tse)$Description <- paste(colData(tse)$Description, "modified description")

# view modified variable
head(tse$Description)
```

```
## [1] "Calhoun South Carolina Pine soil, pH 4.9 modified description"  
## [2] "Cedar Creek Minnesota, grassland, pH 6.1 modified description"  
## [3] "Sevilleta new Mexico, desert scrub, pH 8.3 modified description"
## [4] "M3, Day 1, fecal swab, whole body study modified description"   
## [5] "M1, Day 1, fecal swab, whole body study  modified description"  
## [6] "M3, Day 1, right palm, whole body study modified description"
```
New information can also be added to the experiment by creating a new variable.


```r
# simulate new data
new_data <- runif(ncol(tse))

# store new data as new variable in colData
colData(tse)$NewVariable <- new_data

# view new variable
head(tse$NewVariable)
```

```
## [1] 0.8627 0.2058 0.9080 0.1576 0.7374 0.2652
```

## Merge data

`mia` package has `mergeSEs` function that merges multiple `SummarizedExperiment`
objects. For example, it is possible to combine multiple `TreeSE` objects which each
includes one sample. 

`mergeSEs` works like `dplyr` joining functions. In fact, there are available
`dplyr-like` aliases of `mergeSEs`, such as `full_join`.


```r
# Take subsets for demonstration purposes
tse1 <- tse[, 1]
tse2 <- tse[, 2]
tse3 <- tse[, 3]
tse4 <- tse[1:100, 4]
```


```r
# With inner join, we want to include all shared rows. When using mergeSEs function
# all samples are always preserved.
tse <- mergeSEs(list(tse1, tse2, tse3, tse4), join = "inner")
tse
```

```
## class: TreeSummarizedExperiment 
## dim: 100 4 
## metadata(0):
## assays(1): counts
## rownames(100): 239672 243675 ... 549322 951
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(4): CC1 CL3 M31Fcsw SV1
## colData names(8): X.SampleID Primer ... Description NewVariable
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (100 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```


```r
# Left join preserves all rows of the 1st object
tse <- mia::left_join(tse1, tse4, missing_values = 0)
tse
```

```
## class: TreeSummarizedExperiment 
## dim: 19216 2 
## metadata(0):
## assays(1): counts
## rownames(19216): 239672 243675 ... 239967 254851
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(2): CL3 M31Fcsw
## colData names(8): X.SampleID Primer ... Description NewVariable
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: 1 phylo tree(s) (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

### Additional functions
* [mapTaxonomy](https://microbiome.github.io/mia/reference/taxonomy-methods.html)
* [mergeFeatures/mergeSamples](https://microbiome.github.io/mia/reference/merge-methods.html)
