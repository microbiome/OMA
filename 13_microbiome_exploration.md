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

This chapter focuses on the exploration of microbiome data and establishes 
commonly used descriptors of a microbiome. The main difference to quality
control is that the exploration assumes that technical aspects of the dataset
have been investigated to your satisfaction. Generally speaking, at this point 
you should be quite certain that the dataset doesn't suffer from severe 
technical biases, or you should at least be aware of potential problems.

In reality you might need to go back and forth between QC and exploration, 
since you discover through exploration of your dataset technical aspects you 
need to check.


```r
library(mia)
data("GlobalPatterns", package="mia")
se <- GlobalPatterns 
```

## Prevalence

Prevalence is a measurement, which describes in how many samples certain
microbes were detected.

Investigating the prevalence of microbes allows you either to focus on changes,
which pertain to most of the samples, or to focus on less often found microbes,
which are nonetheless abundantly found in a small number of samples.

On the raw data, the population prevalence (frequency) at a 1% relative
abundance threshold (`detection = 1/100` and `as_relative = TRUE`), can look
like this. The low prevalence in this example can be explained by rather
different sample types as well as the in-depth nature of the data.


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
(`detection = 1`) is accessed.


```r
head(getPrevalence(se, detection = 1, sort = TRUE, abund_values = "counts",
                   as_relative = FALSE))
```

```
## 145149 114821 108747 526804  98605 180658 
##      1      1      1      1      1      1
```

Note that, if the output should be used for subsetting or storing the data in 
the `rowData`, set `sort = FALSE`.

### Prevalent microbiota analysis

To investigate the microbiome data at a selected taxonomic level, two 
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

Alternatively, the `rank` argument can be set to perform the agglomeration on
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

The difference in the naming scheme is that, by default, `na.rm = TRUE` is used
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

Note that the `detection` and `prevalence` thresholds are not the same, since
`detection` can be applied to relative counts or absolute counts depending on 
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

<img src="13_microbiome_exploration_files/figure-html/unnamed-chunk-5-1.png" width="672" />

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
<img src="13_microbiome_exploration_files/figure-html/plot-prev-prev-1.png" alt="Prevalence of top phyla as judged by prevalence" width="672" />
<p class="caption">(\#fig:plot-prev-prev)Prevalence of top phyla as judged by prevalence</p>
</div>

```r
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")
```

<div class="figure">
<img src="13_microbiome_exploration_files/figure-html/plot-prev-mean-1.png" alt="Prevalence of top phyla as judged by mean abundance" width="672" />
<p class="caption">(\#fig:plot-prev-mean)Prevalence of top phyla as judged by mean abundance</p>
</div>

## Abundance

Another approach for exploring data is by visualizing abundance. `miaViz` offers
the function `plotAbundanceDensity` where most abundant taxa can be plotted 
including several options.

In the following few demonstrations are shown, using the 
[@Lahti2014] dataset.

A Jitter plot based on relative abundance data, similar to the one presented at 
[@Salosensaari2021] supplementary figure 1, could be visualized as follows: 


```r
# Loading data
library(microbiomeDataSets)
tse <- atlas1006()

# Counts relative abundances
tse <- transformSamples(tse, method = "relabundance")

plotAbundanceDensity(tse, layout = "jitter", abund_values = "relabundance",
                     n = 40, point_size=1, point_shape=19, point_alpha=0.1) + 
    scale_x_log10(label=scales::percent)
```

<img src="13_microbiome_exploration_files/figure-html/unnamed-chunk-7-1.png" width="672" />

For instance, relative abundance values for the top 5 taxa can be visualized as a
density plot over a log scaled axis, using "nationality" as an overlaying information:


```r
plotAbundanceDensity(tse, layout = "density", abund_values = "relabundance",
                     n = 5, colour_by="nationality", point_alpha=1/10 ) +
    scale_x_log10()
```

<img src="13_microbiome_exploration_files/figure-html/unnamed-chunk-8-1.png" width="672" />
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
 [1] microbiomeDataSets_1.1.1       MultiAssayExperiment_1.19.5   
 [3] miaViz_1.1.4                   ggraph_2.0.5                  
 [5] scater_1.21.3                  ggplot2_3.3.5                 
 [7] scuttle_1.3.1                  mia_1.1.10                    
 [9] TreeSummarizedExperiment_2.1.4 Biostrings_2.61.2             
[11] XVector_0.33.0                 SingleCellExperiment_1.15.1   
[13] SummarizedExperiment_1.23.1    Biobase_2.53.0                
[15] GenomicRanges_1.45.0           GenomeInfoDb_1.29.3           
[17] IRanges_2.27.0                 S4Vectors_0.31.0              
[19] BiocGenerics_0.39.1            MatrixGenerics_1.5.3          
[21] matrixStats_0.60.0             BiocStyle_2.21.3              
[23] rebook_1.3.0                  

loaded via a namespace (and not attached):
  [1] AnnotationHub_3.1.5           BiocFileCache_2.1.1          
  [3] plyr_1.8.6                    igraph_1.2.6                 
  [5] lazyeval_0.2.2                splines_4.1.0                
  [7] BiocParallel_1.27.3           digest_0.6.27                
  [9] htmltools_0.5.1.1             viridis_0.6.1                
 [11] fansi_0.5.0                   magrittr_2.0.1               
 [13] memoise_2.0.0                 ScaledMatrix_1.1.0           
 [15] cluster_2.1.2                 DECIPHER_2.21.0              
 [17] graphlayouts_0.7.1            colorspace_2.0-2             
 [19] rappdirs_0.3.3                blob_1.2.2                   
 [21] ggrepel_0.9.1                 xfun_0.25                    
 [23] dplyr_1.0.7                   crayon_1.4.1                 
 [25] RCurl_1.98-1.3                jsonlite_1.7.2               
 [27] graph_1.71.2                  ape_5.5                      
 [29] glue_1.4.2                    polyclip_1.10-0              
 [31] gtable_0.3.0                  zlibbioc_1.39.0              
 [33] DelayedArray_0.19.1           BiocSingular_1.9.1           
 [35] scales_1.1.1                  DBI_1.1.1                    
 [37] Rcpp_1.0.7                    xtable_1.8-4                 
 [39] viridisLite_0.4.0             decontam_1.13.0              
 [41] tidytree_0.3.4                bit_4.0.4                    
 [43] rsvd_1.0.5                    httr_1.4.2                   
 [45] dir.expiry_1.1.0              ellipsis_0.3.2               
 [47] pkgconfig_2.0.3               XML_3.99-0.6                 
 [49] farver_2.1.0                  dbplyr_2.1.1                 
 [51] CodeDepends_0.6.5             sass_0.4.0                   
 [53] utf8_1.2.2                    AnnotationDbi_1.55.1         
 [55] later_1.2.0                   tidyselect_1.1.1             
 [57] labeling_0.4.2                rlang_0.4.11                 
 [59] reshape2_1.4.4                BiocVersion_3.14.0           
 [61] munsell_0.5.0                 tools_4.1.0                  
 [63] cachem_1.0.5                  DirichletMultinomial_1.35.0  
 [65] generics_0.1.0                RSQLite_2.2.7                
 [67] ExperimentHub_2.1.4           evaluate_0.14                
 [69] stringr_1.4.0                 fastmap_1.1.0                
 [71] yaml_2.2.1                    ggtree_3.1.3                 
 [73] knitr_1.33                    bit64_4.0.5                  
 [75] tidygraph_1.2.0               purrr_0.3.4                  
 [77] KEGGREST_1.33.0               nlme_3.1-152                 
 [79] sparseMatrixStats_1.5.2       mime_0.11                    
 [81] aplot_0.0.6                   compiler_4.1.0               
 [83] png_0.1-7                     interactiveDisplayBase_1.31.2
 [85] curl_4.3.2                    beeswarm_0.4.0               
 [87] filelock_1.0.2                treeio_1.17.2                
 [89] tibble_3.1.3                  tweenr_1.0.2                 
 [91] bslib_0.2.5.1                 stringi_1.7.3                
 [93] highr_0.9                     lattice_0.20-44              
 [95] Matrix_1.3-4                  vegan_2.5-7                  
 [97] permute_0.9-5                 vctrs_0.3.8                  
 [99] pillar_1.6.2                  lifecycle_1.0.0              
[101] BiocManager_1.30.16           jquerylib_0.1.4              
[103] BiocNeighbors_1.11.0          cowplot_1.1.1                
[105] bitops_1.0-7                  irlba_2.3.3                  
[107] httpuv_1.6.1                  patchwork_1.1.1              
[109] R6_2.5.0                      promises_1.2.0.1             
[111] bookdown_0.22                 gridExtra_2.3                
[113] vipor_0.4.5                   codetools_0.2-18             
[115] MASS_7.3-54                   assertthat_0.2.1             
[117] withr_2.4.2                   GenomeInfoDbData_1.2.6       
[119] mgcv_1.8-36                   parallel_4.1.0               
[121] grid_4.1.0                    beachmat_2.9.1               
[123] tidyr_1.1.3                   rmarkdown_2.10               
[125] DelayedMatrixStats_1.15.2     rvcheck_0.1.8                
[127] ggnewscale_0.4.5              ggforce_0.3.3                
[129] shiny_1.6.0                   ggbeeswarm_0.6.0             
```
</div>
