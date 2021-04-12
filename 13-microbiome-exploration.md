# Microbiome Exploration {#microbiome-exploration}

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

This chapter focuses on the exploration of microbiome data and establish 
commonly used descriptors of a microbiome. The main difference to quality
control is that the exploration assumes the technical aspects of the dataset
have been investigated to your satisfaction. Generally speaking at this point 
you should be quite certain, that the dataset doesn't suffer from severe 
technical biases or you should at least be aware of potential problems.

In reality you might need to go back and forth between QC and exploration, 
since you discover through exploration of your dataset technical aspects you 
need to check.


```r
library(mia)
data("GlobalPatterns")
se <- GlobalPatterns 
```

## Prevalence

Prevalence is a measurements, which describes in how many samples certain
microbes were detected.

Investigating the prevalence of microbes allows you either to focus on changes,
which pertain to most of the samples, or to focus on less often found microbes,
which are nonetheless abundantly found in a small number of samples.

On the raw data, the population prevalence (frequency) at a 1% relative
abundance threshold (`detection = 1/100` and `as_relative = TRUE`), can look
like this. The low prevalence in this examples can be explained by rather
different sample types as well as the in depth nature of the data.


```r
head(getPrevalence(se, detection = 1/100, sort = TRUE, as_relative = TRUE))
```

```
## 331820 158660  98605 326977 145149 114821 
## 0.2308 0.2308 0.1923 0.1923 0.1538 0.1538
```

The `detection` and `as_relative` can also be used to access, how many samples
do pass a threshold for raw counts. Here the population prevalence (frequency) 
at the absolute abundance threshold (`as_relative = FALSE`) at read count 1
(`detection = TRUE`) is accessed.


```r
head(getPrevalence(se, detection = 1, sort = TRUE, abund_values = "counts",
                   as_relative = FALSE))
```

```
## 145149 114821 108747 526804  98605 180658 
##      1      1      1      1      1      1
```

Note that, if the output should used for subsetting or storing the data in 
the `rowData`, set `sort = FALSE`.

### Prevalent microbiota analysis

To investigate the microbiome data at a selected taxonomic levels, two 
approaches are available.

First the data can be agglomerated to the taxonomic level and `getPrevalence` 
be used on the result.


```r
altExp(se,"Phylum") <- agglomerateByRank(se, "Phylum")
head(getPrevalence(altExp(se,"Phylum"), detection = 1/100, sort = TRUE,
                   abund_values = "counts", as_relative = TRUE))
```

```
##   Phylum:Bacteroidetes  Phylum:Proteobacteria  Phylum:Actinobacteria 
##                 1.0000                 0.9231                 0.8462 
##   Phylum:Cyanobacteria      Phylum:Firmicutes Phylum:Verrucomicrobia 
##                 0.6154                 0.5769                 0.4615
```

Alternatively the `rank` argument can be set, to perform the agglomeration on
the fly.


```r
altExp(se,"Phylum") <- agglomerateByRank(se, "Phylum")
head(getPrevalence(se, rank = "Phylum", detection = 1/100, sort = TRUE,
                   abund_values = "counts", as_relative = TRUE))
```

```
##   Bacteroidetes  Proteobacteria  Actinobacteria   Cyanobacteria      Firmicutes 
##          1.0000          0.9231          0.8462          0.6154          0.5769 
## Verrucomicrobia 
##          0.4615
```

The difference in the naming scheme, is that by default `na.rm = TRUE` is used
for agglomeration in `getPrevalence`, whereas the default for 
`agglomerateByRank` is `FALSE` to prevent accidental data loss.

If you only need the names of the prevalent taxa, `getPrevalentTaxa` is
available. This returns the taxa that exceed the given prevalence and detection
thresholds.


```r
getPrevalentTaxa(se, detection = 0, prevalence = 50/100)
prev <- getPrevalentTaxa(se, detection = 0, prevalence = 50/100,
                         rank = "Phylum", sort = TRUE)
prev
```

Note, that the `detection` and `prevalence` thresholds are not the same, since
`detection` can be applied to relative counts or absolute couts depending on 
whether `as_relative` is set `TRUE` or `FALSE`

TODO
See also related functions for the analysis of rare and variable taxa
(rareMembers; rareAbundance; lowAbundance). 

### Plotting prevalence

To plot the prevalence, the data is first added to the `rowData`.


```r
rowData(altExp(se,"Phylum"))$prevalence <- 
    getPrevalence(altExp(se,"Phylum"), detection = 1/100, sort = FALSE,
                  abund_values = "counts", as_relative = TRUE)
```

Then it can be plotted via the plotting functions from the `scater` package.
 

```r
library(scater)
plotRowData(altExp(se,"Phylum"), "prevalence", colour_by = "Phylum")
```

<img src="13-microbiome-exploration_files/figure-html/unnamed-chunk-5-1.png" width="672" />

Additionally, the prevalence can be plotted on the taxonomic tree using the
`miaViz` package.


```r
altExps(se) <- splitByRanks(se)
altExps(se) <-
   lapply(altExps(se),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 1/100, sort = FALSE,
                                abund_values = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopTaxa(altExp(se,"Phylum"),
                        method="prevalence",
                        top=10L,
                        abund_values="counts")
top_phyla_mean <- getTopTaxa(altExp(se,"Phylum"),
                             method="mean",
                             top=10L,
                             abund_values="counts")
x <- unsplitByRanks(se, ranks = taxonomyRanks(se)[1:6])
x <- addTaxonomyTree(x)
```
 
After some preparation the data is assembled and can be plotted via 
`plotRowTree`.


```r
library(miaViz)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
```

<div class="figure">
<img src="13-microbiome-exploration_files/figure-html/plot-prev-prev-1.png" alt="Prevalence of top phyla as judged by prevalence" width="672" />
<p class="caption">(\#fig:plot-prev-prev)Prevalence of top phyla as judged by prevalence</p>
</div>

```r
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
```

<div class="figure">
<img src="13-microbiome-exploration_files/figure-html/plot-prev-mean-1.png" alt="Prevalence of top phyla as judged by mean abundance" width="672" />
<p class="caption">(\#fig:plot-prev-mean)Prevalence of top phyla as judged by mean abundance</p>
</div>

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R Under development (unstable) (2021-04-08 r80148)
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
 [1] miaViz_0.99.6                    ggraph_2.0.5                    
 [3] scater_1.19.11                   ggplot2_3.3.3                   
 [5] mia_0.99.10                      TreeSummarizedExperiment_1.99.11
 [7] Biostrings_2.59.2                XVector_0.31.1                  
 [9] SingleCellExperiment_1.13.14     SummarizedExperiment_1.21.3     
[11] Biobase_2.51.0                   GenomicRanges_1.43.4            
[13] GenomeInfoDb_1.27.11             IRanges_2.25.7                  
[15] S4Vectors_0.29.15                BiocGenerics_0.37.1             
[17] MatrixGenerics_1.3.1             matrixStats_0.58.0              
[19] BiocStyle_2.19.2                 rebook_1.1.20                   

loaded via a namespace (and not attached):
  [1] ggtree_2.5.2                ggnewscale_0.4.5           
  [3] ggbeeswarm_0.6.0            colorspace_2.0-0           
  [5] ellipsis_0.3.1              scuttle_1.1.18             
  [7] BiocNeighbors_1.9.4         aplot_0.0.6                
  [9] farver_2.1.0                graphlayouts_0.7.1         
 [11] ggrepel_0.9.1               bit64_4.0.5                
 [13] fansi_0.4.2                 codetools_0.2-18           
 [15] splines_4.1.0               sparseMatrixStats_1.3.7    
 [17] cachem_1.0.4                knitr_1.31                 
 [19] polyclip_1.10-0             jsonlite_1.7.2             
 [21] cluster_2.1.1               graph_1.69.0               
 [23] ggforce_0.3.3               BiocManager_1.30.12        
 [25] compiler_4.1.0              rvcheck_0.1.8              
 [27] lazyeval_0.2.2              assertthat_0.2.1           
 [29] Matrix_1.3-2                fastmap_1.1.0              
 [31] tweenr_1.0.2                BiocSingular_1.7.2         
 [33] htmltools_0.5.1.1           tools_4.1.0                
 [35] igraph_1.2.6                rsvd_1.0.3                 
 [37] gtable_0.3.0                glue_1.4.2                 
 [39] GenomeInfoDbData_1.2.4      dplyr_1.0.5                
 [41] Rcpp_1.0.6                  jquerylib_0.1.3            
 [43] vctrs_0.3.7                 debugme_1.1.0              
 [45] ape_5.4-1                   nlme_3.1-152               
 [47] DECIPHER_2.19.2             DelayedMatrixStats_1.13.5  
 [49] xfun_0.22                   stringr_1.4.0              
 [51] beachmat_2.7.7              lifecycle_1.0.0            
 [53] irlba_2.3.3                 XML_3.99-0.6               
 [55] zlibbioc_1.37.0             MASS_7.3-53.1              
 [57] scales_1.1.1                tidygraph_1.2.0            
 [59] yaml_2.2.1                  memoise_2.0.0              
 [61] gridExtra_2.3               sass_0.3.1                 
 [63] stringi_1.5.3               RSQLite_2.2.6              
 [65] highr_0.8                   ScaledMatrix_0.99.2        
 [67] tidytree_0.3.3              permute_0.9-5              
 [69] filelock_1.0.2              BiocParallel_1.25.5        
 [71] rlang_0.4.10                pkgconfig_2.0.3            
 [73] bitops_1.0-6                evaluate_0.14              
 [75] lattice_0.20-41             purrr_0.3.4                
 [77] treeio_1.15.6               patchwork_1.1.1            
 [79] CodeDepends_0.6.5           labeling_0.4.2             
 [81] cowplot_1.1.1               bit_4.0.4                  
 [83] tidyselect_1.1.0            magrittr_2.0.1             
 [85] bookdown_0.21               R6_2.5.0                   
 [87] generics_0.1.0              DelayedArray_0.17.10       
 [89] DBI_1.1.1                   pillar_1.5.1               
 [91] withr_2.4.1                 mgcv_1.8-34                
 [93] RCurl_1.98-1.3              tibble_3.1.0               
 [95] dir.expiry_0.99.4           crayon_1.4.1               
 [97] utf8_1.2.1                  rmarkdown_2.7              
 [99] viridis_0.5.1               grid_4.1.0                 
[101] blob_1.2.1                  vegan_2.5-7                
[103] digest_0.6.27               tidyr_1.1.3                
[105] munsell_0.5.0               DirichletMultinomial_1.33.2
[107] beeswarm_0.3.1              viridisLite_0.3.0          
[109] vipor_0.4.5                 bslib_0.2.4                
```
</div>
