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

<div class="figure">
<img src="14_alpha_diversity_files/figure-html/plot-div-shannon-1.png" alt="Shannon diversity estimates plotted grouped by sample type." width="672" />
<p class="caption">(\#fig:plot-div-shannon)Shannon diversity estimates plotted grouped by sample type.</p>
</div>

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

<img src="14_alpha_diversity_files/figure-html/visualize-shannon-1.png" width="672" />

**Phylogenetic diversity**  

The phylogenetic diversity is calculated by `mia::estimateDiversity`. This is a faster re-implementation of   
the widely used function in _`picante`_ [@R-picante, @Kembel2010].  

Load `picante` R package and get the `phylo` stored in `rowTree`. 

```r
tse <- mia::estimateDiversity(tse, 
                              abund_values = "counts",
                              index = "faith", 
                              name = "faith")
head(colData(tse)$faith)
```

```
## [1] 250.5 262.3 208.5 117.9 119.8 135.8
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
plots <- lapply(c("observed", "shannon","simpson", "relative", "faith","log_modulo_skewness"),
                plotColData,
                object = tse,
                x = "SampleType",
                colour_by = "SampleType")

plots <- lapply(plots,"+", theme(axis.text.x = element_text(angle=45,hjust=1)))
ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol = 3, common.legend = TRUE, legend = "right")
```

<img src="14_alpha_diversity_files/figure-html/plot-all-diversities-1.png" width="672" />

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
 [1] ggsignif_0.6.3                 scater_1.22.0                 
 [3] ggplot2_3.3.5                  scuttle_1.4.0                 
 [5] mia_1.3.8                      MultiAssayExperiment_1.20.0   
 [7] TreeSummarizedExperiment_2.1.4 Biostrings_2.62.0             
 [9] XVector_0.34.0                 SingleCellExperiment_1.16.0   
[11] SummarizedExperiment_1.24.0    Biobase_2.54.0                
[13] GenomicRanges_1.46.0           GenomeInfoDb_1.30.0           
[15] IRanges_2.28.0                 S4Vectors_0.32.2              
[17] BiocGenerics_0.40.0            MatrixGenerics_1.6.0          
[19] matrixStats_0.61.0-9001        BiocStyle_2.22.0              
[21] rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0            colorspace_2.0-2           
  [3] ellipsis_0.3.2              BiocNeighbors_1.12.0       
  [5] ggpubr_0.4.0                farver_2.1.0               
  [7] ggrepel_0.9.1               bit64_4.0.5                
  [9] fansi_0.5.0                 decontam_1.14.0            
 [11] splines_4.1.2               codetools_0.2-18           
 [13] sparseMatrixStats_1.6.0     cachem_1.0.6               
 [15] knitr_1.36                  jsonlite_1.7.2             
 [17] broom_0.7.10                cluster_2.1.2              
 [19] graph_1.72.0                BiocManager_1.30.16        
 [21] compiler_4.1.2              backports_1.3.0            
 [23] assertthat_0.2.1            Matrix_1.3-4               
 [25] fastmap_1.1.0               lazyeval_0.2.2             
 [27] BiocSingular_1.10.0         htmltools_0.5.2            
 [29] tools_4.1.2                 rsvd_1.0.5                 
 [31] gtable_0.3.0                glue_1.5.0                 
 [33] GenomeInfoDbData_1.2.7      reshape2_1.4.4             
 [35] dplyr_1.0.7                 Rcpp_1.0.7                 
 [37] carData_3.0-4               jquerylib_0.1.4            
 [39] vctrs_0.3.8                 ape_5.5                    
 [41] nlme_3.1-153                DECIPHER_2.22.0            
 [43] DelayedMatrixStats_1.16.0   xfun_0.28                  
 [45] stringr_1.4.0               beachmat_2.10.0            
 [47] lifecycle_1.0.1             irlba_2.3.3                
 [49] rstatix_0.7.0               XML_3.99-0.8               
 [51] zlibbioc_1.40.0             MASS_7.3-54                
 [53] scales_1.1.1                parallel_4.1.2             
 [55] yaml_2.2.1                  memoise_2.0.0              
 [57] gridExtra_2.3               yulab.utils_0.0.4          
 [59] sass_0.4.0                  stringi_1.7.5              
 [61] RSQLite_2.2.8               highr_0.9                  
 [63] ScaledMatrix_1.2.0          permute_0.9-5              
 [65] tidytree_0.3.6              filelock_1.0.2             
 [67] BiocParallel_1.28.0         rlang_0.4.12               
 [69] pkgconfig_2.0.3             bitops_1.0-7               
 [71] evaluate_0.14               lattice_0.20-45            
 [73] purrr_0.3.4                 labeling_0.4.2             
 [75] treeio_1.18.1               CodeDepends_0.6.5          
 [77] cowplot_1.1.1               bit_4.0.4                  
 [79] tidyselect_1.1.1            plyr_1.8.6                 
 [81] magrittr_2.0.1              bookdown_0.24              
 [83] R6_2.5.1                    generics_0.1.1             
 [85] DelayedArray_0.20.0         DBI_1.1.1                  
 [87] withr_2.4.2                 mgcv_1.8-38                
 [89] pillar_1.6.4                abind_1.4-5                
 [91] RCurl_1.98-1.5              tibble_3.1.6               
 [93] dir.expiry_1.2.0            car_3.0-12                 
 [95] crayon_1.4.2                utf8_1.2.2                 
 [97] rmarkdown_2.11              viridis_0.6.2              
 [99] grid_4.1.2                  blob_1.2.2                 
[101] vegan_2.5-7                 digest_0.6.28              
[103] tidyr_1.1.4                 munsell_0.5.0              
[105] DirichletMultinomial_1.36.0 beeswarm_0.4.0             
[107] viridisLite_0.4.0           vipor_0.4.5                
[109] bslib_0.3.1                
```
</div>
