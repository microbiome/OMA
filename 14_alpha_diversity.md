# Community diversity {#community-diversity}

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


Diversity estimates are a central topic in microbiome data analysis. 

There are three commonly employed levels of diversity measurements,
which are trying to put a number on different aspects of the questions
associated with diversity [@Whittaker1960].

Many different ways for estimating such diversity measurements have been 
described in the literature. Which measurement is best or applicable for your 
samples, is not the aim of the following sections.


```r
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
```


**_Alpha diversity_**, also sometimes interchangeably used with the
term **_species diversity_**, summarizes the distribution of species
abundances in a given sample into a single number that depends on
species richness and evenness. Diversity indices measure the overall
community heterogeneity. A number of ecological diversity measures are
available. The Hill coefficient combines many standard indices into a
single equation that provides observed richness, inverse Simpson, and
Shannon diversity, and generalized diversity as special cases. In
general, diversity increases together with increasing richness and
evenness. Sometimes richness, phylogenetic diversity, evenness, dominance, 
and rarity are considered to be variants of alpha diversity.

**Richness** refers to the total number of species in a community
(sample). The simplest richness index is the number of observed
species (observed richness). Assuming limited sampling from the
community, however, this may underestimate the true species
richness. Several estimators are available, including for instance ACE
[@Chao1992] and Chao1 [@Chao1984]. Richness estimates are unaffected
by species abundances.
  
**Phylogenetic diversity** was first proposed by [@Faith1992]. Unlike the 
  diversity measures mentioned above, Phylogenetic diversity (PD) 
  measure incorporates information from phylogenetic relationships 
  stored in `phylo` tree between species in a community (sample). The 
  Faith's PD is calculated as the sum of branch length of all species in 
  a community (sample).

**Evenness** focuses on species abundances, and can thus complement
  the number of species. A typical evenness index is the Pielou's
  evenness, which is Shannon diversity normalized by the observed
  richness.

**Dominance** indices are in general negatively correlated with
  diversity, and sometimes used in ecological literature. High
  dominance is obtained when one or few species have a high share of
  the total species abundance in the community.  
  
**Rarity** indices characterize the concentration of taxa at low abundance. 
  Prevalence and detection thresholds determine rare taxa whose total concentration
  is represented as a rarity index.
  
## Estimation 

Alpha diversity can be estimated with wrapper functions that interact
with other packages implementing the calculation, such as _`vegan`_
[@R-vegan].


### Richness {#richness}

Richness gives the number of features present within a community and can be calculated with `estimateRichness`. Each of the estimate diversity/richness/evenness/dominance functions adds the calculated measure(s) to the `colData` of the `SummarizedExperiment` under the given column `name`. Here, we calculate `observed` features as a measure of richness.     


```r
tse <- mia::estimateRichness(tse, 
                             assay_name = "counts", 
                             index = "observed", 
                             name="observed")

head(colData(tse)$observed)
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr 
##    6964    7679    5729    2667    2574    3214
```
This allows access to the values to be analyzed directly from the `colData`, for example
by plotting them using `plotColData` from the _`scater`_ package [@R-scater].


```r
library(scater)
plotColData(tse, 
            "observed", 
            "SampleType", 
            colour_by = "Final_Barcode") +
    theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  ylab(expression(Richness[Observed]))
```

![(\#fig:plot-div-shannon)Shannon diversity estimates plotted grouped by sample type with colour-labeled barcode.](14_alpha_diversity_files/figure-latex/plot-div-shannon-1.pdf) 

### Diversity {#estimate-diversity}  

The main function, `estimateDiversity`, calculates the selected
diversity index based on the selected assay data.  


```r
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "shannon", 
                              name = "shannon")
head(colData(tse)$shannon)
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr 
##   6.577   6.777   6.498   3.828   3.288   4.289
```

Alpha diversities can be visualized with boxplot. Here, Shannon index is compared 
between different sample type groups. Individual data points are visualized by 
plotting them as points with `geom_jitter`.

`geom_signif` is used to test whether these differences are statistically significant.
It adds p-values to plot.


```r
if( !require(ggsignif) ){
  install.packages(ggsignif)
}
library(ggplot2)
library(patchwork)
library(ggsignif)

# Subsets the data. Takes only those samples that are from feces, skin, or tongue,
# and creates data frame from the collected data
df <- as.data.frame(colData(tse)[colData(tse)$SampleType %in% 
                                  c("Feces", "Skin", "Tongue"), ])

# Changes old levels with new levels
df$SampleType <- factor(df$SampleType)

# For significance testing, all different combinations are determined
comb <- split(t(combn(levels(df$SampleType), 2)), 
           seq(nrow(t(combn(levels(df$SampleType), 2)))))

ggplot(df, aes(x = SampleType, y = shannon)) +
  # Outliers are removed, because otherwise each data point would be plotted twice; 
  # as an outlier of boxplot and as a point of dotplot.
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb, map_signif_level = FALSE) +
  theme(text = element_text(size = 10))
```

![](14_alpha_diversity_files/figure-latex/visualize-shannon-1.pdf)<!-- --> 

### Faith phylogenetic diversity {#faith-diversity}

The Faith index is returned by the function `estimateFaith`.


```r
tse <- mia::estimateFaith(tse,
                          assay_name = "counts")
head(colData(tse)$faith)
```

```
## [1] 250.5 262.3 208.5 117.9 119.8 135.8
```

**Note**: because `tse` is a `TreeSummarizedExperiment` object, its phylogenetic tree is used by default. However, the optional argument `tree` must be provided if `tse` does not contain one.

Below a visual comparison between shannon and faith indices is shown with a violin plot.


```r
plots <- lapply(c("shannon", "faith"),
                plotColData,
                object = tse, colour_by = "SampleType")
plots[[1]] + plots[[2]] +
  plot_layout(guides = "collect")
```

![](14_alpha_diversity_files/figure-latex/phylo-div-2-1.pdf)<!-- --> 
 
Alternatively, the phylogenetic diversity can be calculated by `mia::estimateDiversity`. This is a faster re-implementation of   
the widely used function in _`picante`_ [@R-picante, @Kembel2010].  

Load `picante` R package and get the `phylo` stored in `rowTree`.


```r
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "faith", 
                              name = "faith")
```

### Evenness  

Evenness can be calculated with `estimateEvenness`.  


```r
tse <- estimateEvenness(tse, 
                        assay_name = "counts", 
                        index="simpson")
head(colData(tse)$simpson)
```

```
## [1] 0.026871 0.027197 0.047049 0.005179 0.004304 0.005011
```


### Dominance  

Dominance can be calculated with `estimateDominance`. Here, the `Relative index` is calculated which is the relative abundance of the most dominant species in the sample.   


```r
tse <- estimateDominance(tse, 
                         assay_name = "counts", 
                         index="relative")

head(colData(tse)$relative)
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr 
## 0.03910 0.03226 0.01690 0.22981 0.21778 0.22329
```

### Rarity  

`mia` package provides one rarity index called log-modulo skewness. It can be 
calculated with `estimateDiversity`.


```r
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "log_modulo_skewness")

head(colData(tse)$log_modulo_skewness)
```

```
## [1] 2.061 2.061 2.061 2.061 2.061 2.061
```

### Divergence

Divergence can be evaluated with `estimateDivergence`. Reference and algorithm for the calculation of divergence can be specified as `reference` and `FUN`, respectively. 


```r
tse <- mia::estimateDivergence(tse,
                               assay_name = "counts",
                               reference = "median",
                               FUN = vegan::vegdist)
```


## Visualization

A plot comparing all the diversity measures calculated above and stored in `colData` can then be constructed directly.


```r
plots <- lapply(c("observed", "shannon", "simpson", "relative", "faith", "log_modulo_skewness"),
                plotColData,
                object = tse,
                x = "SampleType",
                colour_by = "SampleType")

plots <- lapply(plots, "+", 
                theme(axis.text.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.ticks.x = element_blank()))

((plots[[1]] | plots[[2]] | plots[[3]]) / 
(plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect")
```

![](14_alpha_diversity_files/figure-latex/plot-all-diversities-1.pdf)<!-- --> 

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
 [1] patchwork_1.1.2                ggsignif_0.6.4                
 [3] scater_1.26.1                  ggplot2_3.4.0                 
 [5] scuttle_1.8.0                  mia_1.5.17                    
 [7] MultiAssayExperiment_1.24.0    TreeSummarizedExperiment_2.1.4
 [9] Biostrings_2.66.0              XVector_0.38.0                
[11] SingleCellExperiment_1.20.0    SummarizedExperiment_1.28.0   
[13] Biobase_2.58.0                 GenomicRanges_1.50.1          
[15] GenomeInfoDb_1.34.3            IRanges_2.32.0                
[17] S4Vectors_0.36.0               BiocGenerics_0.44.0           
[19] MatrixGenerics_1.10.0          matrixStats_0.63.0-9003       
[21] BiocStyle_2.24.0               rebook_1.6.0                  

loaded via a namespace (and not attached):
 [1] ggbeeswarm_0.6.0            colorspace_2.0-3           
 [3] BiocNeighbors_1.16.0        farver_2.1.1               
 [5] ggrepel_0.9.2               bit64_4.0.5                
 [7] fansi_1.0.3                 decontam_1.18.0            
 [9] codetools_0.2-18            splines_4.2.1              
[11] sparseMatrixStats_1.10.0    cachem_1.0.6               
[13] knitr_1.41                  jsonlite_1.8.3             
[15] cluster_2.1.4               graph_1.74.0               
[17] BiocManager_1.30.19         compiler_4.2.1             
[19] assertthat_0.2.1            Matrix_1.5-3               
[21] fastmap_1.1.0               lazyeval_0.2.2             
[23] cli_3.4.1                   BiocSingular_1.14.0        
[25] htmltools_0.5.3             tools_4.2.1                
[27] rsvd_1.0.5                  gtable_0.3.1               
[29] glue_1.6.2                  GenomeInfoDbData_1.2.9     
[31] reshape2_1.4.4              dplyr_1.0.10               
[33] Rcpp_1.0.9                  vctrs_0.5.1                
[35] ape_5.6-2                   nlme_3.1-160               
[37] DECIPHER_2.26.0             DelayedMatrixStats_1.20.0  
[39] xfun_0.35                   stringr_1.4.1              
[41] beachmat_2.14.0             lifecycle_1.0.3            
[43] irlba_2.3.5.1               XML_3.99-0.12              
[45] zlibbioc_1.44.0             MASS_7.3-58.1              
[47] scales_1.2.1                parallel_4.2.1             
[49] yaml_2.3.6                  memoise_2.0.1              
[51] gridExtra_2.3               yulab.utils_0.0.5          
[53] stringi_1.7.8               RSQLite_2.2.19             
[55] highr_0.9                   ScaledMatrix_1.6.0         
[57] tidytree_0.4.1              permute_0.9-7              
[59] filelock_1.0.2              BiocParallel_1.32.1        
[61] rlang_1.0.6                 pkgconfig_2.0.3            
[63] bitops_1.0-7                evaluate_0.18              
[65] lattice_0.20-45             purrr_0.3.5                
[67] labeling_0.4.2              treeio_1.22.0              
[69] CodeDepends_0.6.5           cowplot_1.1.1              
[71] bit_4.0.5                   tidyselect_1.2.0           
[73] plyr_1.8.8                  magrittr_2.0.3             
[75] bookdown_0.30               R6_2.5.1                   
[77] generics_0.1.3              DelayedArray_0.24.0        
[79] DBI_1.1.3                   withr_2.5.0                
[81] pillar_1.8.1                mgcv_1.8-41                
[83] RCurl_1.98-1.9              tibble_3.1.8               
[85] dir.expiry_1.4.0            crayon_1.5.2               
[87] utf8_1.2.2                  rmarkdown_2.18             
[89] viridis_0.6.2               grid_4.2.1                 
[91] blob_1.2.3                  vegan_2.6-4                
[93] digest_0.6.30               tidyr_1.2.1                
[95] munsell_0.5.0               DirichletMultinomial_1.40.0
[97] beeswarm_0.4.0              viridisLite_0.4.1          
[99] vipor_0.4.5                
```
</div>
