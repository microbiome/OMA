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
data("GlobalPatterns", package="mia")
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
the widely used function in _`picante`_ [@R-picante, @Kembel2010].  

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
plots <- lapply(c("observed", "shannon","simpson", "relative", "faith","log_modulo_skewness"),
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
R version 4.1.1 (2021-08-10)
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
 [1] ggsignif_0.6.3                 scater_1.21.8                 
 [3] ggplot2_3.3.5                  scuttle_1.3.1                 
 [5] mia_1.1.17                     TreeSummarizedExperiment_2.1.4
 [7] Biostrings_2.61.2              XVector_0.33.0                
 [9] SingleCellExperiment_1.15.2    SummarizedExperiment_1.23.5   
[11] Biobase_2.53.0                 GenomicRanges_1.45.0          
[13] GenomeInfoDb_1.29.8            IRanges_2.27.2                
[15] S4Vectors_0.31.5               BiocGenerics_0.39.2           
[17] MatrixGenerics_1.5.4           matrixStats_0.61.0-9001       
[19] BiocStyle_2.21.4               rebook_1.3.1                  

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                backports_1.2.1            
  [3] plyr_1.8.6                  lazyeval_0.2.2             
  [5] splines_4.1.1               BiocParallel_1.27.17       
  [7] digest_0.6.28               htmltools_0.5.2            
  [9] viridis_0.6.2               fansi_0.5.0                
 [11] magrittr_2.0.1              memoise_2.0.0              
 [13] ScaledMatrix_1.1.0          cluster_2.1.2              
 [15] DECIPHER_2.21.0             openxlsx_4.2.4             
 [17] colorspace_2.0-2            blob_1.2.2                 
 [19] ggrepel_0.9.1               haven_2.4.3                
 [21] xfun_0.26                   dplyr_1.0.7                
 [23] crayon_1.4.1                RCurl_1.98-1.5             
 [25] jsonlite_1.7.2              graph_1.71.2               
 [27] ape_5.5                     glue_1.4.2                 
 [29] gtable_0.3.0                zlibbioc_1.39.0            
 [31] DelayedArray_0.19.4         car_3.0-11                 
 [33] BiocSingular_1.9.1          abind_1.4-5                
 [35] scales_1.1.1                DBI_1.1.1                  
 [37] rstatix_0.7.0               Rcpp_1.0.7                 
 [39] viridisLite_0.4.0           decontam_1.13.0            
 [41] tidytree_0.3.5              foreign_0.8-81             
 [43] bit_4.0.4                   rsvd_1.0.5                 
 [45] dir.expiry_1.1.0            ellipsis_0.3.2             
 [47] pkgconfig_2.0.3             XML_3.99-0.8               
 [49] farver_2.1.0                CodeDepends_0.6.5          
 [51] sass_0.4.0                  utf8_1.2.2                 
 [53] tidyselect_1.1.1            labeling_0.4.2             
 [55] rlang_0.4.11                reshape2_1.4.4             
 [57] cellranger_1.1.0            munsell_0.5.0              
 [59] tools_4.1.1                 cachem_1.0.6               
 [61] DirichletMultinomial_1.35.0 generics_0.1.0             
 [63] RSQLite_2.2.8               broom_0.7.9                
 [65] evaluate_0.14               stringr_1.4.0              
 [67] fastmap_1.1.0               yaml_2.2.1                 
 [69] knitr_1.36                  bit64_4.0.5                
 [71] zip_2.2.0                   purrr_0.3.4                
 [73] nlme_3.1-153                sparseMatrixStats_1.5.3    
 [75] compiler_4.1.1              beeswarm_0.4.0             
 [77] filelock_1.0.2              curl_4.3.2                 
 [79] treeio_1.17.2               tibble_3.1.5               
 [81] bslib_0.3.1                 stringi_1.7.5              
 [83] highr_0.9                   forcats_0.5.1              
 [85] lattice_0.20-45             Matrix_1.3-4               
 [87] vegan_2.5-7                 permute_0.9-5              
 [89] vctrs_0.3.8                 pillar_1.6.3               
 [91] lifecycle_1.0.1             BiocManager_1.30.16        
 [93] jquerylib_0.1.4             BiocNeighbors_1.11.0       
 [95] data.table_1.14.2           cowplot_1.1.1              
 [97] bitops_1.0-7                irlba_2.3.3                
 [99] R6_2.5.1                    bookdown_0.24              
[101] gridExtra_2.3               rio_0.5.27                 
[103] vipor_0.4.5                 codetools_0.2-18           
[105] MASS_7.3-54                 assertthat_0.2.1           
[107] withr_2.4.2                 GenomeInfoDbData_1.2.7     
[109] mgcv_1.8-38                 parallel_4.1.1             
[111] hms_1.1.1                   grid_4.1.1                 
[113] beachmat_2.9.1              tidyr_1.1.4                
[115] rmarkdown_2.11              DelayedMatrixStats_1.15.4  
[117] carData_3.0-4               ggpubr_0.4.0               
[119] ggbeeswarm_0.6.0           
```
</div>
