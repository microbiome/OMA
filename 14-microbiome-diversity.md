# Microbiome Diversity {#microbiome-diversity}

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

## Alpha diversity

**_Alpha diversity_**, also sometimes interchangeably used with the
term **_species diversity_**, summarizes the distribution of species
abundances in a given sample into a single number that depends on
species richness and evenness. Diversity indices measure the overall
community heterogeneity. A number of ecological diversity measures are
available. The Hill coefficient combines many standard indices into a
single equation that provides observed richness, inverse Simpson, and
Shannon diversity, and generalized diversity as special cases. In
general, diversity increases together with increasing richness and
evenness. Sometimes richness, evenness, and dominance are considered
to be variants of alpha diversity.

**Richness** refers to the total number of species in a community
(sample). The simplest richness index is the number of observed
species (observed richness). Assuming limited sampling from the
community, however, this may underestimate the true species
richness. Several estimators are available, including for instance ACE
[@reference] and Chao1 [@reference]. Richness estimates are unaffected
by species abundances.

**Evenness** focuses on species abundances, and can thus complement
  the number of species. A typical evenness index is the Pielou's
  evenness, which is Shannon diversity normalized by the observed
  richness.

**Dominance** indices are in general negatively correlated with
  diversity, and sometimes used in ecological literature. High
  dominance is obtained when one or few species have a high share of
  the total species abundance in the community.


### Estimating alpha diversity

Alpha diversity can be estimated with wrapper functions that interact
with other packages implementing the calculation, such as _`vegan`_
[@R-vegan].

The main function, `estimateAlphaDiversity`, calculates the selected
diversity index based on the selected assay data and adds it to the
`colData` of the `SummarizedExperiment` under the given column `name`.


```r
se <- mia::estimateAlphaDiversity(se, abund_values = "counts",
                             index = "shannon", name = "shannon")
head(colData(se)$shannon)
```

```
##     CL3     CC1     SV1 M31Fcsw M11Fcsw M31Plmr 
##   6.577   6.777   6.498   3.828   3.288   4.289
```

This allows the values to analyzed directly from the `colData`, for example
by plotting them using `plotColData` from the _`scater`_ package [@R-scater].


```r
library(scater)
plotColData(se, "shannon", "SampleType", colour_by = "SampleType") +
    theme(axis.text.x = element_text(angle=45,hjust=1))
```

<div class="figure">
<img src="14-microbiome-diversity_files/figure-html/plot-div-shannon-1.png" alt="Shannon diversity estimates plotted grouped by sample type." width="672" />
<p class="caption">(\#fig:plot-div-shannon)Shannon diversity estimates plotted grouped by sample type.</p>
</div>

All available indices will be calculated by default...


```r
se <- estimateAlphaDiversity(se)
```

.. and a plot comparing them can then be constructed directly.


```r
plots <- lapply(c("shannon","simpson","inv_simpson", "richness", "chao1", "ACE"),
                plotColData,
                object = se,
                x = "SampleType",
                colour_by = "SampleType")
plots <- lapply(plots,"+", theme(axis.text.x = element_text(angle=45,hjust=1)))
ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol = 3, common.legend = TRUE, legend = "right")
```

<img src="14-microbiome-diversity_files/figure-html/unnamed-chunk-4-1.png" width="672" />

## Beta diversity

Where alpha diversity focuses on community variation within a
community (sample), beta diversity quantifies (dis-)similarites
between communities (samples).

Technically, beta diversities are usually represented as `dist`
objects, which contain triangular data describing the distance between
each pair of samples. These distances can be further subjected to
ordination with methods such as multi-dimensional scaling (also known
as PCoA) to retrieve reduced dimensions for further evaluation or
visualization. [TODO we could here instead link to ordination chapter?]


### Estimating beta diversity

In the following examples distances are calculated by variable
functions supplied to the `FUN` argument. The function can defined by
the user. It must return a `dist` function, which can then be used to
calculate reduced dimension either via ordination methods (such as MDS
or NMDS), and the results can be stored in the `reducedDim`.

This entire process is wrapped in the `runMDS` and `runNMDS`
functions.


```r
se <- runMDS2(se, FUN = vegan::vegdist, name = "MDS_BC")
```

Sample similarities can be visualized on a lower-dimensional display
(typically 2D) using the `plotReducedDim` function in the
`scater`package. This provides also further tools to incorporate
additional information using variations in colour, shape or size.


```r
plotReducedDim(se, "MDS_BC", colour_by = "SampleType")
```

<div class="figure">
<img src="14-microbiome-diversity_files/figure-html/plot-mds-bray-curtis-1.png" alt="MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset." width="672" />
<p class="caption">(\#fig:plot-mds-bray-curtis)MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset.</p>
</div>

With additional tools from the `ggplot2` universe, comparisons can be 
performed informing on the applicability to visualize sample similarities in a 
meaningful way.


```r
se <- runMDS2(se, FUN = vegan::vegdist, name = "MDS_euclidean",
             method = "euclidean")
se <- runNMDS(se, FUN = vegan::vegdist, name = "NMDS_BC")
```

```
## initial  value 47.733208 
## iter   5 value 33.853364
## iter  10 value 32.891200
## final  value 32.823570 
## converged
```

```r
se <- runNMDS(se, FUN = vegan::vegdist, name = "NMDS_euclidean",
              method = "euclidean")
```

```
## initial  value 31.882673 
## final  value 31.882673 
## converged
```

```r
plots <- lapply(c("MDS_BC","MDS_euclidean","NMDS_BC","NMDS_euclidean"),
                plotReducedDim, object = se, colour_by = "SampleType")
ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol = 2, common.legend = TRUE,
                  legend = "right")
```

<div class="figure">
<img src="14-microbiome-diversity_files/figure-html/plot-mds-nmds-comparison-1.png" alt="Comparison of MDS and NMDS plots based on the Bray-Curtis or euclidean distances on the GlobalPattern dataset." width="672" />
<p class="caption">(\#fig:plot-mds-nmds-comparison)Comparison of MDS and NMDS plots based on the Bray-Curtis or euclidean distances on the GlobalPattern dataset.</p>
</div>

The _UniFrac_ method is a special case, as it requires data on the
relationship of features in form on a `phylo` tree. `calculateUniFrac`
performs the calculation to return a `dist` object, which can again be
used within `runMDS`.



```r
se <- runMDS2(se, FUN = calculateUniFrac, name = "UniFrac",
              tree = rowTree(se),
              ntop = nrow(se))
```


```r
plotReducedDim(se, "UniFrac", colour_by = "SampleType")
```

<div class="figure">
<img src="14-microbiome-diversity_files/figure-html/plot-unifrac-1.png" alt="UniFrac distances scaled by MDS of the GlobalPattern dataset." width="672" />
<p class="caption">(\#fig:plot-unifrac)UniFrac distances scaled by MDS of the GlobalPattern dataset.</p>
</div>

### Other ordination methods

Other dimension reduction methods, such as `PCA`, `t-SNE` and `UMAP` are 
inherited directly from the `scater` package.


```r
se <- runPCA(se, name = "PCA", exprs_values = "counts", ncomponents = 10)
```


```r
plotReducedDim(se, "PCA", colour_by = "SampleType")
```

<div class="figure">
<img src="14-microbiome-diversity_files/figure-html/plot-pca-1.png" alt="PCA plot on the GlobalPatterns data set containing sample from different sources." width="672" />
<p class="caption">(\#fig:plot-pca)PCA plot on the GlobalPatterns data set containing sample from different sources.</p>
</div>

As mentioned before, applicability of the different methods depends on your
sample set.

FIXME: let us switch to UMAP for the examples?


```r
se <- runTSNE(se, name = "TSNE", exprs_values = "counts", ncomponents = 3)
```


```r
plotReducedDim(se, "TSNE", colour_by = "SampleType", ncomponents = c(1:3))
```

<div class="figure">
<img src="14-microbiome-diversity_files/figure-html/plot-tsne-1.png" alt="t-SNE plot on the GlobalPatterns data set containing sample from different sources." width="672" />
<p class="caption">(\#fig:plot-tsne)t-SNE plot on the GlobalPatterns data set containing sample from different sources.</p>
</div>



## Community comparisons [TODO combine with the material above for simplicity?]


A typical comparison of community composition starts with a visual
comparison of the groups on a 2D ordination.

Let us load an example data set:


```r
library(microbiomeDataSets)
se.lahti <- LahtiMData()
```


Then we estimate relative abundances and MDS ordination based on
Bray-Curtis (BC) dissimilarity between the groups, and visualize the
results.



```r
se.lahti <- relAbundanceCounts(se.lahti)
se.lahti <- runNMDS(se.lahti, FUN = vegan::vegdist, name = "BC", nmdsFUN = "monoMDS",
                    exprs_values = "relabundance",
                    keep_dist = TRUE)
```


```r
plotReducedDim(se.lahti, "BC", colour_by = "group")
```

<img src="14-microbiome-diversity_files/figure-html/unnamed-chunk-11-1.png" width="672" />

No clear difference between the groups can be visually observed.


### Testing differences in community composition between sample groups

The permutational analysis of variance (PERMANOVA) [@Anderson2001] is
a widely used non-parametric multivariate method that can be used to
estimate the actual statistical significance of differences in the
observed community composition between two groups of
samples.

PERMANOVA evaluates the hypothesis that the centroids and dispersion
of the community are equivalent between the compared groups. A small
p-value indicates that the compared groups have, on average, a
different community composition.

This method is implemented in the `vegan` package (function `adonis`).


```r
library(vegan)
permanova <- vegan::adonis(t(assay(se.lahti,"relabundance")) ~ group,
                           data = colData(se.lahti),
                           permutations = 9999)

# P-value
print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])
```

```
## [1] 0.2737
```

In this case, the community composition is not significantly different
between the groups.

Let us visualize the model coefficients for species that exhibit the
largest differences between the groups. This gives some insights into
how the groups tend to differ from each other in terms of community
composition.



```r
coef <- coefficients(permanova)["group1",]
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))
```





```r
ggplot(data.frame(x = top.coef,
                  y = factor(names(top.coef),
                                      unique(names(top.coef)))),
        aes(x = x, y = y)) +
    geom_bar(stat="identity") +
    labs(x="",y="",title="Top Taxa") +
    theme_bw()
```

<img src="14-microbiome-diversity_files/figure-html/plot-top-coef-anova-1.png" width="672" />

In the above example, the largest differences between the two groups
can be attributed to Bacteroides intestinalis (elevated in the first
group) and Faecalibacterium prausnitzii (elevated in the second
group), and many other co-varying species.



### Checking the homogeneity condition 

It is important to note that the application of PERMANOVA assumes
homogeneous group dispersions (variances). This can be tested with the
PERMDISP2 method [@Anderson2006].



```r
anova(vegan::betadisper(attr(reducedDim(se.lahti,"BC"),"dist"),
                        colData(se.lahti)$group))
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df Sum Sq Mean Sq F value Pr(>F)
## Groups     1  0.000 0.00002       0   0.95
## Residuals 42  0.158 0.00376
```

In our example, the groups have similar dispersions, and PERMANOVA is
an appropriate choice for comparing community compositions.


## Further reading

In certain settings beta diversities might be used to group samples without
prior knowledge. For this we want to point to excellent resources on 
[how to extract information from the clusters](http://bioconductor.org/books/release/OSCA/clustering.html).

See also [community typing](15-microbiome-community.md).

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R Under development (unstable) (2021-02-25 r80035)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.1 LTS

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
 [1] vegan_2.5-7                      lattice_0.20-41                 
 [3] permute_0.9-5                    microbiomeDataSets_0.99.5       
 [5] MultiAssayExperiment_1.17.13     scater_1.19.11                  
 [7] ggplot2_3.3.3                    mia_0.98.40                     
 [9] TreeSummarizedExperiment_1.99.11 Biostrings_2.59.2               
[11] XVector_0.31.1                   SingleCellExperiment_1.13.12    
[13] SummarizedExperiment_1.21.1      Biobase_2.51.0                  
[15] GenomicRanges_1.43.3             GenomeInfoDb_1.27.6             
[17] IRanges_2.25.6                   S4Vectors_0.29.7                
[19] BiocGenerics_0.37.1              MatrixGenerics_1.3.1            
[21] matrixStats_0.58.0               BiocStyle_2.19.1                
[23] rebook_1.1.15                    BiocManager_1.30.10             

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                  backports_1.2.1              
  [3] AnnotationHub_2.23.2          BiocFileCache_1.15.1         
  [5] splines_4.1.0                 BiocParallel_1.25.4          
  [7] digest_0.6.27                 htmltools_0.5.1.1            
  [9] viridis_0.5.1                 fansi_0.4.2                  
 [11] magrittr_2.0.1                memoise_2.0.0                
 [13] ScaledMatrix_0.99.2           cluster_2.1.1                
 [15] DECIPHER_2.19.2               openxlsx_4.2.3               
 [17] colorspace_2.0-0              rappdirs_0.3.3               
 [19] blob_1.2.1                    haven_2.3.1                  
 [21] xfun_0.21                     dplyr_1.0.4                  
 [23] callr_3.5.1                   crayon_1.4.1                 
 [25] RCurl_1.98-1.2                jsonlite_1.7.2               
 [27] graph_1.69.0                  ape_5.4-1                    
 [29] glue_1.4.2                    gtable_0.3.0                 
 [31] zlibbioc_1.37.0               DelayedArray_0.17.9          
 [33] car_3.0-10                    BiocSingular_1.7.2           
 [35] Rhdf5lib_1.13.4               HDF5Array_1.19.8             
 [37] abind_1.4-5                   scales_1.1.1                 
 [39] DBI_1.1.1                     rstatix_0.7.0                
 [41] Rcpp_1.0.6                    xtable_1.8-4                 
 [43] viridisLite_0.3.0             foreign_0.8-81               
 [45] bit_4.0.4                     rsvd_1.0.3                   
 [47] httr_1.4.2                    ellipsis_0.3.1               
 [49] pkgconfig_2.0.3               XML_3.99-0.5                 
 [51] farver_2.1.0                  scuttle_1.1.18               
 [53] dbplyr_2.1.0                  CodeDepends_0.6.5            
 [55] sass_0.3.1                    utf8_1.1.4                   
 [57] AnnotationDbi_1.53.1          later_1.1.0.1                
 [59] tidyselect_1.1.0              labeling_0.4.2               
 [61] rlang_0.4.10                  BiocVersion_3.13.1           
 [63] munsell_0.5.0                 cellranger_1.1.0             
 [65] tools_4.1.0                   cachem_1.0.4                 
 [67] ExperimentHub_1.17.1          DirichletMultinomial_1.33.2  
 [69] generics_0.1.0                RSQLite_2.2.3                
 [71] broom_0.7.5                   evaluate_0.14                
 [73] stringr_1.4.0                 fastmap_1.1.0                
 [75] yaml_2.2.1                    processx_3.4.5               
 [77] knitr_1.31                    bit64_4.0.5                  
 [79] zip_2.1.1                     purrr_0.3.4                  
 [81] KEGGREST_1.31.1               nlme_3.1-152                 
 [83] sparseMatrixStats_1.3.6       mime_0.10                    
 [85] debugme_1.1.0                 compiler_4.1.0               
 [87] png_0.1-7                     interactiveDisplayBase_1.29.0
 [89] beeswarm_0.2.3                filelock_1.0.2               
 [91] curl_4.3                      ggsignif_0.6.1               
 [93] tibble_3.1.0                  bslib_0.2.4                  
 [95] stringi_1.5.3                 highr_0.8                    
 [97] ps_1.6.0                      forcats_0.5.1                
 [99] Matrix_1.3-2                  vctrs_0.3.6                  
[101] rhdf5filters_1.3.4            pillar_1.5.0                 
[103] lifecycle_1.0.0               jquerylib_0.1.3              
[105] BiocNeighbors_1.9.4           data.table_1.14.0            
[107] cowplot_1.1.1                 bitops_1.0-6                 
[109] irlba_2.3.3                   httpuv_1.5.5                 
[111] R6_2.5.0                      promises_1.2.0.1             
[113] bookdown_0.21                 gridExtra_2.3                
[115] rio_0.5.26                    vipor_0.4.5                  
[117] codetools_0.2-18              MASS_7.3-53.1                
[119] assertthat_0.2.1              rhdf5_2.35.2                 
[121] withr_2.4.1                   GenomeInfoDbData_1.2.4       
[123] mgcv_1.8-34                   hms_1.0.0                    
[125] grid_4.1.0                    beachmat_2.7.7               
[127] tidyr_1.1.3                   rmarkdown_2.7                
[129] DelayedMatrixStats_1.13.5     carData_3.0-4                
[131] Rtsne_0.15                    ggpubr_0.4.0                 
[133] shiny_1.6.0                   ggbeeswarm_0.6.0             
```
</div>
