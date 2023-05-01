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
  install.packages("ggsignif")
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
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] patchwork_1.1.2                ggsignif_0.6.4                
 [3] scater_1.28.0                  ggplot2_3.4.2                 
 [5] scuttle_1.10.0                 mia_1.9.2                     
 [7] MultiAssayExperiment_1.26.0    TreeSummarizedExperiment_2.1.4
 [9] Biostrings_2.68.0              XVector_0.40.0                
[11] SingleCellExperiment_1.22.0    SummarizedExperiment_1.30.0   
[13] Biobase_2.60.0                 GenomicRanges_1.52.0          
[15] GenomeInfoDb_1.36.0            IRanges_2.34.0                
[17] S4Vectors_0.38.0               BiocGenerics_0.46.0           
[19] MatrixGenerics_1.12.0          matrixStats_0.63.0-9003       
[21] BiocStyle_2.28.0               rebook_1.9.0                  

loaded via a namespace (and not attached):
 [1] DBI_1.1.3                   bitops_1.0-7               
 [3] gridExtra_2.3               permute_0.9-7              
 [5] CodeDepends_0.6.5           rlang_1.1.1                
 [7] magrittr_2.0.3              compiler_4.3.0             
 [9] RSQLite_2.3.1               mgcv_1.8-42                
[11] dir.expiry_1.8.0            DelayedMatrixStats_1.22.0  
[13] vctrs_0.6.2                 reshape2_1.4.4             
[15] stringr_1.5.0               pkgconfig_2.0.3            
[17] crayon_1.5.2                fastmap_1.1.1              
[19] labeling_0.4.2              utf8_1.2.3                 
[21] rmarkdown_2.21              graph_1.78.0               
[23] ggbeeswarm_0.7.2            DirichletMultinomial_1.42.0
[25] purrr_1.0.1                 bit_4.0.5                  
[27] xfun_0.39                   zlibbioc_1.46.0            
[29] cachem_1.0.7                beachmat_2.16.0            
[31] jsonlite_1.8.4              blob_1.2.4                 
[33] highr_0.10                  DelayedArray_0.25.0        
[35] BiocParallel_1.34.0         cluster_2.1.4              
[37] irlba_2.3.5.1               parallel_4.3.0             
[39] R6_2.5.1                    stringi_1.7.12             
[41] Rcpp_1.0.10                 bookdown_0.33              
[43] knitr_1.42                  DECIPHER_2.28.0            
[45] splines_4.3.0               Matrix_1.5-4               
[47] tidyselect_1.2.0            yaml_2.3.7                 
[49] viridis_0.6.2               vegan_2.6-4                
[51] codetools_0.2-19            lattice_0.21-8             
[53] tibble_3.2.1                plyr_1.8.8                 
[55] withr_2.5.0                 treeio_1.24.0              
[57] evaluate_0.20               pillar_1.9.0               
[59] BiocManager_1.30.20         filelock_1.0.2             
[61] generics_0.1.3              RCurl_1.98-1.12            
[63] sparseMatrixStats_1.12.0    munsell_0.5.0              
[65] scales_1.2.1                tidytree_0.4.2             
[67] glue_1.6.2                  lazyeval_0.2.2             
[69] tools_4.3.0                 BiocNeighbors_1.18.0       
[71] ScaledMatrix_1.7.1          XML_3.99-0.14              
[73] cowplot_1.1.1               grid_4.3.0                 
[75] tidyr_1.3.0                 ape_5.7-1                  
[77] colorspace_2.1-0            nlme_3.1-162               
[79] GenomeInfoDbData_1.2.10     beeswarm_0.4.0             
[81] BiocSingular_1.16.0         vipor_0.4.5                
[83] cli_3.6.1                   rsvd_1.0.5                 
[85] fansi_1.0.4                 viridisLite_0.4.1          
[87] dplyr_1.1.2                 gtable_0.3.3               
[89] yulab.utils_0.0.6           digest_0.6.31              
[91] ggrepel_0.9.3               farver_2.1.1               
[93] decontam_1.20.0             memoise_2.0.1              
[95] htmltools_0.5.5             lifecycle_1.0.3            
[97] bit64_4.0.5                 MASS_7.3-59                
```
</div>
