# Introduction {#intro}


## Installation  

Install the development version from GitHub using `remotes` R package.  


```r
# install remotes 
#install.packages("remotes")
remotes::install_github("FelixErnst/mia")
```

```
## Microbiom... (3012ee51f... -> cd40180ab...) [GitHub]
## generics     (0.0.2        -> 0.1.0       ) [CRAN]
## labeling     (0.3          -> 0.4.2       ) [CRAN]
## cpp11        (0.2.2        -> 0.2.4       ) [CRAN]
## 
##      checking for file ‘/tmp/Rtmpltr73W/remotes34b7ad6b18e/FelixErnst-MicrobiomeExperiment-cd40180/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/Rtmpltr73W/remotes34b7ad6b18e/FelixErnst-MicrobiomeExperiment-cd40180/DESCRIPTION’
##   ─  preparing ‘MicrobiomeExperiment’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
## ─  looking to see if a ‘data/datalist’ file should be added
##   ─  building ‘MicrobiomeExperiment_0.99.0.9014.tar.gz’ (480ms)
##      
##      checking for file ‘/tmp/Rtmpltr73W/remotes34b79f1cc91/FelixErnst-mia-d5ba1c2/DESCRIPTION’ ...  ✔  checking for file ‘/tmp/Rtmpltr73W/remotes34b79f1cc91/FelixErnst-mia-d5ba1c2/DESCRIPTION’
##   ─  preparing ‘mia’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
## ─  looking to see if a ‘data/datalist’ file should be added
##   ─  building ‘mia_0.0.0.9002.tar.gz’ (451ms)
##      
## 
```

## Package data  

The datasets in `mia` are conversions of the `phyloseq` datasets `GlobalPatterns` `enterotype`, `esophagus` and `soilrep`. 
Taxa contains example taxonomic data derived from the metagenomeFeatures package.


```r
library(mia)
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```
## Loading required package: SingleCellExperiment
```

```
## Loading required package: TreeSummarizedExperiment
```

```
## Loading required package: MicrobiomeExperiment
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```r
# Example how to load data
data("GlobalPatterns")

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
## altExpNames(0):
## rowLinks: a LinkDataFrame (19216 rows)
## rowTree: a phylo (19216 leaves)
## colLinks: NULL
## colTree: NULL
```

## Tools    

1. `mia`    : Microbiome analysis tools   
2. `miaViz` : Visualization support for `mia`   


**See also:**    

[`TreeSummarizedExperiment`](https://www.bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html)  
[`MicrobiomeExperiment`](https://github.com/FelixErnst/MicrobiomeExperiment)    [`SummarizedExperiment`](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)  
[`microbiome`](https://bioconductor.org/packages/devel/bioc/html/microbiome.html)    




