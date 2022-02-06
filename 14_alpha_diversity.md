# Community diversity

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


### Richness  

Richness gives the number of features present within a community and can be calculated with `estimateRichness`. Each of the estimate diversity/richness/evenness/dominance functions adds the calculated measure(s) to the `colData` of the `SummarizedExperiment` under the given column `name`. Here, we calculate `observed` features as a measure of richness.     


```r
tse <- mia::estimateRichness(tse, 
                             abund_values = "counts", 
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
            colour_by = "SampleType") +
    theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  ylab(expression(Richness[Observed]))
```

![(\#fig:plot-div-shannon)Shannon diversity estimates plotted grouped by sample type.](14_alpha_diversity_files/figure-latex/plot-div-shannon-1.pdf) 

### Diversity  

**Non-Phylogenetic measures**  
The main function, `estimateDiversity`, calculates the selected
diversity index based on the selected assay data.  


```r
tse <- mia::estimateDiversity(tse, 
                              abund_values = "counts",
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

### Faith phylogenetic diversity

The Faith index is returned by the function `estimateFaith`.


```r
tse <- mia::estimateFaith(tse,
                          abund_values = "counts")
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
                object = tse)
plots[[1]] + plots[[2]]
```

![](14_alpha_diversity_files/figure-latex/phylo-div-2-1.pdf)<!-- --> 
 
Alternatively, the phylogenetic diversity can be calculated by `mia::estimateDiversity`. This is a faster re-implementation of   
the widely used function in _`picante`_ [@R-picante, @Kembel2010].  

Load `picante` R package and get the `phylo` stored in `rowTree`.


```r
tse <- mia::estimateDiversity(tse, 
                              abund_values = "counts",
                              index = "faith", 
                              name = "faith")
```

### Evenness  

Evenness can be calculated with `estimateEvenness`.  


```r
tse <- estimateEvenness(tse, 
                        abund_values = "counts", 
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
                         abund_values = "counts", 
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
                              abund_values = "counts",
                              index = "log_modulo_skewness")

head(colData(tse)$log_modulo_skewness)
```

```
## [1] 2.061 2.061 2.061 2.061 2.061 2.061
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
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

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
 [1] patchwork_1.1.1                ggsignif_0.6.3                
 [3] scater_1.22.0                  ggplot2_3.3.5                 
 [5] scuttle_1.4.0                  mia_1.3.15                    
 [7] MultiAssayExperiment_1.20.0    TreeSummarizedExperiment_2.1.4
 [9] Biostrings_2.62.0              XVector_0.34.0                
[11] SingleCellExperiment_1.16.0    SummarizedExperiment_1.24.0   
[13] Biobase_2.54.0                 GenomicRanges_1.46.1          
[15] GenomeInfoDb_1.30.1            IRanges_2.28.0                
[17] S4Vectors_0.32.3               BiocGenerics_0.40.0           
[19] MatrixGenerics_1.6.0           matrixStats_0.61.0-9001       
[21] BiocStyle_2.22.0               rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0            colorspace_2.0-2           
  [3] ellipsis_0.3.2              BiocNeighbors_1.12.0       
  [5] farver_2.1.0                ggrepel_0.9.1              
  [7] bit64_4.0.5                 fansi_1.0.2                
  [9] decontam_1.14.0             splines_4.1.2              
 [11] codetools_0.2-18            sparseMatrixStats_1.6.0    
 [13] cachem_1.0.6                knitr_1.37                 
 [15] jsonlite_1.7.3              cluster_2.1.2              
 [17] graph_1.72.0                BiocManager_1.30.16        
 [19] compiler_4.1.2              assertthat_0.2.1           
 [21] Matrix_1.4-0                fastmap_1.1.0              
 [23] lazyeval_0.2.2              cli_3.1.1                  
 [25] BiocSingular_1.10.0         htmltools_0.5.2            
 [27] tools_4.1.2                 rsvd_1.0.5                 
 [29] gtable_0.3.0                glue_1.6.1                 
 [31] GenomeInfoDbData_1.2.7      reshape2_1.4.4             
 [33] dplyr_1.0.7                 Rcpp_1.0.8                 
 [35] vctrs_0.3.8                 ape_5.6-1                  
 [37] nlme_3.1-155                DECIPHER_2.22.0            
 [39] DelayedMatrixStats_1.16.0   xfun_0.29                  
 [41] stringr_1.4.0               beachmat_2.10.0            
 [43] lifecycle_1.0.1             irlba_2.3.5                
 [45] XML_3.99-0.8                zlibbioc_1.40.0            
 [47] MASS_7.3-55                 scales_1.1.1               
 [49] parallel_4.1.2              yaml_2.2.2                 
 [51] memoise_2.0.1               gridExtra_2.3              
 [53] yulab.utils_0.0.4           stringi_1.7.6              
 [55] RSQLite_2.2.9               highr_0.9                  
 [57] ScaledMatrix_1.2.0          tidytree_0.3.7             
 [59] permute_0.9-7               filelock_1.0.2             
 [61] BiocParallel_1.28.3         rlang_1.0.0                
 [63] pkgconfig_2.0.3             bitops_1.0-7               
 [65] evaluate_0.14               lattice_0.20-45            
 [67] purrr_0.3.4                 labeling_0.4.2             
 [69] treeio_1.18.1               CodeDepends_0.6.5          
 [71] cowplot_1.1.1               bit_4.0.4                  
 [73] tidyselect_1.1.1            plyr_1.8.6                 
 [75] magrittr_2.0.2              bookdown_0.24              
 [77] R6_2.5.1                    generics_0.1.2             
 [79] DelayedArray_0.20.0         DBI_1.1.2                  
 [81] withr_2.4.3                 mgcv_1.8-38                
 [83] pillar_1.7.0                RCurl_1.98-1.5             
 [85] tibble_3.1.6                dir.expiry_1.2.0           
 [87] crayon_1.4.2                utf8_1.2.2                 
 [89] rmarkdown_2.11              viridis_0.6.2              
 [91] grid_4.1.2                  blob_1.2.2                 
 [93] vegan_2.5-7                 digest_0.6.29              
 [95] tidyr_1.2.0                 munsell_0.5.0              
 [97] DirichletMultinomial_1.36.0 beeswarm_0.4.0             
 [99] viridisLite_0.4.0           vipor_0.4.5                
```
</div>
