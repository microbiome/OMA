# Alpha Diversity

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
data("GlobalPatterns")
se <- GlobalPatterns
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
  
**Phylogenetic diversity** was first proposed by [@Faith1992], unlike the 
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
  
## Estimating alpha diversity

Alpha diversity can be estimated with wrapper functions that interact
with other packages implementing the calculation, such as _`vegan`_
[@R-vegan].

### Richness  

Richness gives the number of features present within a community and can be calculated with `estimateRichness`. Each of the estimate diversity/richness/evenness/dominance functions adds the calculated measure(s) to the `colData` of the `SummarizedExperiment` under the given column `name`. Here, we calculate `observed` features as a measure of richness.     


```r
se <- mia::estimateRichness(se, 
                       abund_values = "counts", 
                       index = "observed", 
                       name="observed")

head(colData(se)$observed)
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr 
##    6964    7679    5729    2667    2574    3214
```
This allows access to the values to be analyzed directly from the `colData`, for example
by plotting them using `plotColData` from the _`scater`_ package [@R-scater].


```r
library(scater)
plotColData(se, 
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
se <- mia::estimateDiversity(se, 
                             abund_values = "counts",
                             index = "shannon", 
                             name = "shannon")
head(colData(se)$shannon)
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr 
##   6.577   6.777   6.498   3.828   3.288   4.289
```

Alpha diversities can be visualized with boxplot. Here, Shannon index is compared 
between different sample type groups. Individual data points are visualized by 
plotting them as points with `geom_jitter`.

`geom_signif` is used to test, if these differences are statistically significant.
It adds p-values to plot.


```r
if( !require(ggsignif) ){
  install.packages(ggsignif)
}
library(ggplot2)
library(ggsignif)

# Subsets the data. Takes only those samples that are from feces, skin, or tongue,
# and creates data frame from the collected data
df <- as.data.frame(colData(se)[colData(se)$SampleType %in% 
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
the widely function in _`picante`_ [@R-picante, @Kembel2010].  

Load `picante` R package and get the `phylo` stored in `rowTree`. 

```r
se <- mia::estimateDiversity(se, 
                             abund_values = "counts",
                             index = "faith", 
                             name = "faith")
head(colData(se)$faith)
```

```
## [1] 250.5 262.3 208.5 117.9 119.8 135.8
```

### Evenness  

Evenness can be calculated with `estimateEvenness`.  


```r
se <- estimateEvenness(se, 
                       abund_values = "counts", 
                       index="simpson")
head(colData(se)$simpson)
```

```
## [1] 0.026871 0.027197 0.047049 0.005179 0.004304 0.005011
```


### Dominance  

Dominance can be calculated with `estimateDominance`. Here, the `Relative index` is calculated which is the relative abundance of the most dominant species in the sample.   


```r
se <- estimateDominance(se, 
                       abund_values = "counts", 
                       index="relative")

head(colData(se)$relative)
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr 
## 0.03910 0.03226 0.01690 0.22981 0.21778 0.22329
```

### Rarity  

`mia` package provides one rarity index called log-modulo skewness. It can be 
calculated with `estimateDiversity`.


```r
se <- mia::estimateDiversity(se, 
                             abund_values = "counts",
                             index = "log_modulo_skewness")

head(colData(se)$log_modulo_skewness)
```

```
## [1] 2.061 2.061 2.061 2.061 2.061 2.061
```


## Visualize alpha diversities  

A plot comparing all the diversity measures calculated above and stored in `colData` can then be constructed directly.  

```r
plots <- lapply(c("observed", "shannon","simpson", "relative", "faith"),
                plotColData,
                object = se,
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
 [1] ggsignif_0.6.2                 scater_1.21.2                 
 [3] ggplot2_3.3.5                  scuttle_1.3.0                 
 [5] mia_1.1.7                      TreeSummarizedExperiment_2.1.3
 [7] Biostrings_2.61.1              XVector_0.33.0                
 [9] SingleCellExperiment_1.15.1    SummarizedExperiment_1.23.1   
[11] Biobase_2.53.0                 GenomicRanges_1.45.0          
[13] GenomeInfoDb_1.29.3            IRanges_2.27.0                
[15] S4Vectors_0.31.0               BiocGenerics_0.39.1           
[17] MatrixGenerics_1.5.1           matrixStats_0.60.0            
[19] BiocStyle_2.21.3               rebook_1.3.0                  

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                backports_1.2.1            
  [3] plyr_1.8.6                  lazyeval_0.2.2             
  [5] splines_4.1.0               BiocParallel_1.27.2        
  [7] digest_0.6.27               htmltools_0.5.1.1          
  [9] viridis_0.6.1               fansi_0.5.0                
 [11] magrittr_2.0.1              memoise_2.0.0              
 [13] ScaledMatrix_1.1.0          cluster_2.1.2              
 [15] DECIPHER_2.21.0             openxlsx_4.2.4             
 [17] colorspace_2.0-2            blob_1.2.2                 
 [19] haven_2.4.1                 xfun_0.24                  
 [21] dplyr_1.0.7                 crayon_1.4.1               
 [23] RCurl_1.98-1.3              jsonlite_1.7.2             
 [25] graph_1.71.2                ape_5.5                    
 [27] glue_1.4.2                  gtable_0.3.0               
 [29] zlibbioc_1.39.0             DelayedArray_0.19.1        
 [31] car_3.0-11                  BiocSingular_1.9.1         
 [33] abind_1.4-5                 scales_1.1.1               
 [35] DBI_1.1.1                   rstatix_0.7.0              
 [37] Rcpp_1.0.7                  viridisLite_0.4.0          
 [39] decontam_1.13.0             tidytree_0.3.4             
 [41] foreign_0.8-81              bit_4.0.4                  
 [43] rsvd_1.0.5                  dir.expiry_1.1.0           
 [45] ellipsis_0.3.2              pkgconfig_2.0.3            
 [47] XML_3.99-0.6                farver_2.1.0               
 [49] CodeDepends_0.6.5           sass_0.4.0                 
 [51] utf8_1.2.2                  tidyselect_1.1.1           
 [53] labeling_0.4.2              rlang_0.4.11               
 [55] reshape2_1.4.4              cellranger_1.1.0           
 [57] munsell_0.5.0               tools_4.1.0                
 [59] cachem_1.0.5                DirichletMultinomial_1.35.0
 [61] generics_0.1.0              RSQLite_2.2.7              
 [63] broom_0.7.8                 evaluate_0.14              
 [65] stringr_1.4.0               fastmap_1.1.0              
 [67] yaml_2.2.1                  knitr_1.33                 
 [69] bit64_4.0.5                 zip_2.2.0                  
 [71] purrr_0.3.4                 nlme_3.1-152               
 [73] sparseMatrixStats_1.5.0     compiler_4.1.0             
 [75] beeswarm_0.4.0              filelock_1.0.2             
 [77] curl_4.3.2                  treeio_1.17.2              
 [79] tibble_3.1.3                bslib_0.2.5.1              
 [81] stringi_1.7.3               highr_0.9                  
 [83] forcats_0.5.1               lattice_0.20-44            
 [85] Matrix_1.3-4                vegan_2.5-7                
 [87] permute_0.9-5               vctrs_0.3.8                
 [89] pillar_1.6.1                lifecycle_1.0.0            
 [91] BiocManager_1.30.16         jquerylib_0.1.4            
 [93] BiocNeighbors_1.11.0        data.table_1.14.0          
 [95] cowplot_1.1.1               bitops_1.0-7               
 [97] irlba_2.3.3                 R6_2.5.0                   
 [99] bookdown_0.22               gridExtra_2.3              
[101] rio_0.5.27                  vipor_0.4.5                
[103] codetools_0.2-18            MASS_7.3-54                
[105] assertthat_0.2.1            withr_2.4.2                
[107] GenomeInfoDbData_1.2.6      mgcv_1.8-36                
[109] parallel_4.1.0              hms_1.1.0                  
[111] grid_4.1.0                  beachmat_2.9.0             
[113] tidyr_1.1.3                 rmarkdown_2.9              
[115] DelayedMatrixStats_1.15.0   carData_3.0-4              
[117] ggpubr_0.4.0                ggbeeswarm_0.6.0           
```
</div>
