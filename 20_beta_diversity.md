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




# Community similarity

Where alpha diversity focuses on community variation within a
community (sample), beta diversity quantifies (dis-)similarites
between communities (samples). Some of the most popular beta diversity
measures in microbiome research include Bray-Curtis index (for
compositional data), Jaccard index (for presence / absence data,
ignoring abundance information), Aitchison distance (Euclidean
distance for clr transformed abundances, aiming to avoid the
compositionality bias), and the Unifrac distances (that take into
account the phylogenetic tree information). Only some of the commonly
used beta diversity measures are actual _distances_; this is a
mathematically well-defined concept and many ecological beta diversity
measures, such as Bray-Curtis index, are not proper distances.
Therefore, the term dissimilarity or beta diversity is commonly used.

Technically, beta diversities are usually represented as `dist`
objects, which contain triangular data describing the distance between
each pair of samples. These distances can be further subjected to
ordination. Ordination is a common concept in ecology that aims to
reduce the dimensionality of the data for further evaluation or
visualization. Ordination techniques aim to capture as much of
essential information in the data as possible in a lower dimensional
representation.  Dimension reduction is bound to loose information but
the common ordination techniques aim to preserve relevant information
of sample similarities in an optimal way, which is defined in
different ways by different methods. [TODO add references and/or link
to ordination chapter instead?]

Some of the most common ordination methods in microbiome research
include Principal Component Analysis (PCA), metric and non-metric
multi-dimensional scaling (MDS, NMDS), The MDS methods are also known
as Principal Coordinates Analysis (PCoA). Other recently popular
techniques include t-SNE and UMAP. 


## Explained variance

The percentage of explained variance is typically shown for PCA
ordination plots. This quantifies the proportion of overall variance
in the data that is captured by the PCA axes, or how well the
ordination axes reflect the original distances.

Sometimes a similar measure is shown for MDS/PCoA. The interpretation
is generally different, however, and hence we do not recommend using
it. PCA is a special case of PCoA with Euclidean distances.  With
non-Euclidean dissimilarities PCoA uses a trick where the pointwise
dissimilarities are first cast into similarities in a Euclidean space
(with some information loss i.e. stress) and then projected to the
maximal variance axes. In this case, the maximal variance axes do not
directly reflect the correspondence of the projected distances and
original distances, as they do for PCA.

In typical use cases, we would like to know how well the ordination
reflects the original similarity structures; then the quantity of
interest is the so-called "stress" function, which measures the
difference in pairwise similarities between the data points in the
original (high-dimensional) vs. projected (low-dimensional) space.

Hence, we propose that for PCoA and other ordination methods, users
would report relative stress (varies in the unit interval; the smaller
the better). This can be calculated as shown below. For further
examples, check the [note from Huber
lab](https://www.huber.embl.de/users/klaus/Teaching/statisticalMethods-lab.pdf).



```r
# Example data
library(mia)
data(GlobalPatterns, package="mia")

# Data matrix (features x samples)
tse <- GlobalPatterns
tse <- transformCounts(tse, method = "relabundance")

# Add group information Feces yes/no
colData(tse)$Group <- colData(tse)$SampleType=="Feces"

# Quantify dissimilarities in the original feature space
library(vegan)
x <- assay(tse, "relabundance") # Pick relabunance assay separately
d0 <- as.matrix(vegdist(t(x), "bray"))

# PCoA Ordination
pcoa <- as.data.frame(cmdscale(d0, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

# Quantify dissimilarities in the ordination space
dp <- as.matrix(dist(pcoa))

# Calculate stress i.e. relative difference in the original and
# projected dissimilarities
stress <- sum((dp - d0)^2)/sum(d0^2)
```


Shepard plot visualizes the original versus projected (ordination)
dissimilarities between the data points:


```r
ord <- order(as.vector(d0))
df <- data.frame(d0 = as.vector(d0)[ord],
                  dmds = as.vector(dp)[ord])

library(ggplot2)
ggplot(aes(x = d0, y = dmds), data=df) + 
       geom_smooth() +
       geom_point() +       
       labs(title = "Shepard plot",
       x = "Original distance",
       y = "MDS distance",       
            subtitle = paste("Stress:", round(stress, 2))) +
  theme_bw()
```

<img src="20_beta_diversity_files/figure-html/shepard-1.png" width="3000" />


## Community comparisons by beta diversity analysis

A typical comparison of community composition starts with a visual
comparison of the groups on a 2D ordination.

Then we estimate relative abundances and MDS ordination based on
Bray-Curtis (BC) dissimilarity between the groups, and visualize the
results.

In the following examples dissimilarities are calculated by 
functions supplied to the `FUN` argument. This function can be defined by
the user. It must return a `dist` function, which can then be used to
calculate reduced dimensions either via ordination methods (such as MDS
or NMDS), and the results can be stored in the `reducedDim`.

This entire process is wrapped in the `runMDS` and `runNMDS`
functions.


```r
library(scater)
tse <- runMDS(tse, FUN = vegan::vegdist, name = "PCoA_BC", exprs_values = "counts")
```

Sample similarities can be visualized on a lower-dimensional display
(typically 2D) using the `plotReducedDim` function in the `scater`
package. This provides also further tools to incorporate additional
information using variations in color, shape or size. Are there
visible differences between the groups?


```r
# Create ggplot object
p <- plotReducedDim(tse, "PCoA_BC", colour_by = "Group")

# Add explained variance for each axis
e <- attr(reducedDim(tse, "PCoA_BC"), "eig");
rel_eig <- e/sum(e[e>0])		  
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))

print(p)
```

<div class="figure">
<img src="20_beta_diversity_files/figure-html/plot-mds-bray-curtis-1.png" alt="MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset." width="3000" />
<p class="caption">(\#fig:plot-mds-bray-curtis)MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset.</p>
</div>




With additional tools from the `ggplot2` universe, comparisons can be 
performed informing on the applicability to visualize sample similarities in a 
meaningful way.


```r
tse <- runMDS(tse, FUN = vegan::vegdist, name = "MDS_euclidean",
             method = "euclidean", exprs_values = "counts")
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "NMDS_BC")
```

```
## initial  value 47.733208 
## iter   5 value 33.853364
## iter  10 value 32.891200
## final  value 32.823570 
## converged
```

```r
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "NMDS_euclidean",
               method = "euclidean")
```

```
## initial  value 31.882673 
## final  value 31.882673 
## converged
```

```r
plots <- lapply(c("PCoA_BC","MDS_euclidean","NMDS_BC","NMDS_euclidean"),
                plotReducedDim, object = tse, colour_by = "Group")
ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol = 2, common.legend = TRUE,
                  legend = "right")
```

<div class="figure">
<img src="20_beta_diversity_files/figure-html/plot-mds-nmds-comparison-1.png" alt="Comparison of MDS and NMDS plots based on the Bray-Curtis or euclidean distances on the GlobalPattern dataset." width="3000" />
<p class="caption">(\#fig:plot-mds-nmds-comparison)Comparison of MDS and NMDS plots based on the Bray-Curtis or euclidean distances on the GlobalPattern dataset.</p>
</div>

The _Unifrac_ method is a special case, as it requires data on the
relationship of features in form on a `phylo` tree. `calculateUnifrac`
performs the calculation to return a `dist` object, which can again be
used within `runMDS`.



```r
library(scater)
tse <- runMDS(tse, FUN = mia::calculateUnifrac, name = "Unifrac",
              tree = rowTree(tse),
              ntop = nrow(tse),
             exprs_values = "counts")
```


```r
plotReducedDim(tse, "Unifrac", colour_by = "Group")
```

<div class="figure">
<img src="20_beta_diversity_files/figure-html/plot-unifrac-1.png" alt="Unifrac distances scaled by MDS of the GlobalPattern dataset." width="3000" />
<p class="caption">(\#fig:plot-unifrac)Unifrac distances scaled by MDS of the GlobalPattern dataset.</p>
</div>

## Other ordination methods

Other dimension reduction methods, such as `PCA`, `t-SNE` and `UMAP` are 
inherited directly from the `scater` package.


```r
tse <- runPCA(tse, name = "PCA", exprs_values = "counts", ncomponents = 10)
```


```r
plotReducedDim(tse, "PCA", colour_by = "Group")
```

<div class="figure">
<img src="20_beta_diversity_files/figure-html/plot-pca-1.png" alt="PCA plot on the GlobalPatterns data set containing sample from different sources." width="3000" />
<p class="caption">(\#fig:plot-pca)PCA plot on the GlobalPatterns data set containing sample from different sources.</p>
</div>

As mentioned before, applicability of the different methods depends on your
sample set.

FIXME: let us switch to UMAP for the examples?


```r
tse <- runTSNE(tse, name = "TSNE", exprs_values = "counts", ncomponents = 3)
```


```r
plotReducedDim(tse, "TSNE", colour_by = "Group", ncomponents = c(1:3))
```

<div class="figure">
<img src="20_beta_diversity_files/figure-html/plot-tsne-1.png" alt="t-SNE plot on the GlobalPatterns data set containing sample from different sources." width="3000" />
<p class="caption">(\#fig:plot-tsne)t-SNE plot on the GlobalPatterns data set containing sample from different sources.</p>
</div>



## Visualizing the most dominant genus on PCoA

In this section we visualize most dominant genus on PCoA. A similar visualization was proposed in [Taxonomic signatures of cause-specific mortality risk in human gut microbiome](https://www.nature.com/articles/s41467-021-22962-y), Salosensaari et al. (2021).


Let us agglomerate the data at a Genus level and getting the dominant taxa per sample.


```r
# Agglomerate to genus level
tse_Genus <- agglomerateByRank(tse, rank="Genus")
# Convert to relative abundances
tse_Genus <- transformSamples(tse, method = "relabundance", abund_values="counts")
# Add info on dominant genus per sample
tse_Genus <- addPerSampleDominantTaxa(tse_Genus, abund_values="relabundance", name = "dominant_taxa")
```


Performing PCoA with Bray-Curtis dissimilarity.

```r
tse_Genus <- runMDS(tse_Genus, FUN = vegan::vegdist,
              name = "PCoA_BC", exprs_values = "relabundance")
```


Getting top taxa and visualizing the abundance on PCoA.


```r
# Getting the top taxa
top_taxa <- getTopTaxa(tse_Genus,top = 6, abund_values = "relabundance")

# Naming all the rest of non top-taxa as "Other"
most_abundant <- lapply(colData(tse_Genus)$dominant_taxa,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})

# Storing the previous results as a new column within colData
colData(tse_Genus)$most_abundant <- as.character(most_abundant)

# Calculating percentage of the most abundant
most_abundant_freq <- table(as.character(most_abundant))
most_abundant_percent <- round(most_abundant_freq/sum(most_abundant_freq)*100, 1)

# Retrieving the explained variance
e <- attr(reducedDim(tse_Genus, "PCoA_BC"), "eig");
var_explained <- e/sum(e[e>0])*100

# Visualization
plot <-plotReducedDim(tse_Genus,"PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(most_abundant_percent),"(",most_abundant_percent,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       color="")
plot
```

<img src="20_beta_diversity_files/figure-html/unnamed-chunk-7-1.png" width="3000" />

Note: A 3D interactive version of the earlier plot can be found from [here](https://microbiome.github.io/OMA/interactive-3d-plots.html).

Similarly let's visualize and compare the sub-population.


```r
# Calculating the frequencies and percentages for both categories
freq_TRUE <- table(as.character(most_abundant[colData(tse_Genus)$Group==TRUE]))
freq_FALSE <- table(as.character(most_abundant[colData(tse_Genus)$Group==FALSE]))
percent_TRUE <- round(freq_TRUE/sum(freq_TRUE)*100, 1)
percent_FALSE <- round(freq_FALSE/sum(freq_FALSE)*100, 1)

# Visualization
plotReducedDim(tse_Genus[,colData(tse_Genus)$Group==TRUE],
                          "PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(percent_TRUE),"(",percent_TRUE,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       title = "Group = TRUE", color="")
```

<img src="20_beta_diversity_files/figure-html/unnamed-chunk-8-1.png" width="3000" />

```r
plotReducedDim(tse_Genus[,colData(tse_Genus)$Group==FALSE],
                          "PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(percent_FALSE),"(",percent_FALSE,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       title = "Group = FALSE", color="")
```

<img src="20_beta_diversity_files/figure-html/unnamed-chunk-8-2.png" width="3000" />




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
permanova <- vegan::adonis(t(assay(tse,"relabundance")) ~ Group,
                           data = colData(tse),
                           permutations = 9999)

# P-value
print(as.data.frame(permanova$aov.tab)["Group", "Pr(>F)"])
```

```
## [1] 5e-04
```

In this case, the community composition is not significantly different
between the groups.

Let us visualize the model coefficients for species that exhibit the
largest differences between the groups. This gives some insights into
how the groups tend to differ from each other in terms of community
composition.



```r
coef <- coefficients(permanova)["Group1",]
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

<img src="20_beta_diversity_files/figure-html/plot-top-coef-anova-1.png" width="3000" />

In the above example, the largest differences between the two groups
can be attributed to _Bacteroides intestinalis_ (elevated in the first
group) and _Faecalibacterium prausnitzii_ (elevated in the second
group), and many other co-varying species.



### Checking the homogeneity condition 

It is important to note that the application of PERMANOVA assumes
homogeneous group dispersions (variances). This can be tested with the
PERMDISP2 method [@Anderson2006] by using the same assay and distance
method than in PERMANOVA.



```r
anova(vegan::betadisper(vegan::vegdist(t(assay(tse, "counts"))),
                        colData(tse)$Group))
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df Sum Sq Mean Sq F value  Pr(>F)    
## Groups     1 0.1607  0.1607     125 5.5e-11 ***
## Residuals 24 0.0309  0.0013                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

If the groups have similar dispersion, PERMANOVA can be seen as an
appropriate choice for comparing community compositions.


## Further reading


 - [How to extract information from clusters](http://bioconductor.org/books/release/OSCA/clustering.html)

 - Chapter \@ref(community-typing) on community typing

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
 [1] scater_1.22.0                  scuttle_1.4.0                 
 [3] ggplot2_3.3.5                  vegan_2.5-7                   
 [5] lattice_0.20-45                permute_0.9-5                 
 [7] mia_1.3.13                     MultiAssayExperiment_1.20.0   
 [9] TreeSummarizedExperiment_2.1.4 Biostrings_2.62.0             
[11] XVector_0.34.0                 SingleCellExperiment_1.16.0   
[13] SummarizedExperiment_1.24.0    Biobase_2.54.0                
[15] GenomicRanges_1.46.1           GenomeInfoDb_1.30.0           
[17] IRanges_2.28.0                 S4Vectors_0.32.3              
[19] BiocGenerics_0.40.0            MatrixGenerics_1.6.0          
[21] matrixStats_0.61.0-9001        BiocStyle_2.22.0              
[23] rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] Rtsne_0.15                  ggbeeswarm_0.6.0           
  [3] colorspace_2.0-2            ggsignif_0.6.3             
  [5] ellipsis_0.3.2              BiocNeighbors_1.12.0       
  [7] ggpubr_0.4.0                farver_2.1.0               
  [9] ggrepel_0.9.1               bit64_4.0.5                
 [11] fansi_1.0.0                 decontam_1.14.0            
 [13] splines_4.1.2               codetools_0.2-18           
 [15] sparseMatrixStats_1.6.0     cachem_1.0.6               
 [17] knitr_1.37                  jsonlite_1.7.2             
 [19] broom_0.7.11                cluster_2.1.2              
 [21] graph_1.72.0                BiocManager_1.30.16        
 [23] compiler_4.1.2              backports_1.4.1            
 [25] assertthat_0.2.1            Matrix_1.4-0               
 [27] fastmap_1.1.0               lazyeval_0.2.2             
 [29] BiocSingular_1.10.0         htmltools_0.5.2            
 [31] tools_4.1.2                 rsvd_1.0.5                 
 [33] gtable_0.3.0                glue_1.6.0                 
 [35] GenomeInfoDbData_1.2.7      reshape2_1.4.4             
 [37] dplyr_1.0.7                 Rcpp_1.0.7                 
 [39] carData_3.0-5               jquerylib_0.1.4            
 [41] vctrs_0.3.8                 ape_5.6-1                  
 [43] nlme_3.1-153                DECIPHER_2.22.0            
 [45] DelayedMatrixStats_1.16.0   xfun_0.29                  
 [47] stringr_1.4.0               beachmat_2.10.0            
 [49] lifecycle_1.0.1             irlba_2.3.5                
 [51] rstatix_0.7.0               XML_3.99-0.8               
 [53] zlibbioc_1.40.0             MASS_7.3-54                
 [55] scales_1.1.1                parallel_4.1.2             
 [57] yaml_2.2.1                  memoise_2.0.1              
 [59] gridExtra_2.3               yulab.utils_0.0.4          
 [61] sass_0.4.0                  stringi_1.7.6              
 [63] RSQLite_2.2.9               highr_0.9                  
 [65] ScaledMatrix_1.2.0          tidytree_0.3.6             
 [67] filelock_1.0.2              BiocParallel_1.28.3        
 [69] rlang_0.4.12                pkgconfig_2.0.3            
 [71] bitops_1.0-7                evaluate_0.14              
 [73] purrr_0.3.4                 labeling_0.4.2             
 [75] treeio_1.18.1               CodeDepends_0.6.5          
 [77] cowplot_1.1.1               bit_4.0.4                  
 [79] tidyselect_1.1.1            plyr_1.8.6                 
 [81] magrittr_2.0.1              bookdown_0.24              
 [83] R6_2.5.1                    generics_0.1.1             
 [85] DelayedArray_0.20.0         DBI_1.1.2                  
 [87] withr_2.4.3                 mgcv_1.8-38                
 [89] pillar_1.6.4                abind_1.4-5                
 [91] RCurl_1.98-1.5              tibble_3.1.6               
 [93] dir.expiry_1.2.0            car_3.0-12                 
 [95] crayon_1.4.2                utf8_1.2.2                 
 [97] rmarkdown_2.11              viridis_0.6.2              
 [99] grid_4.1.2                  blob_1.2.2                 
[101] digest_0.6.29               tidyr_1.1.4                
[103] munsell_0.5.0               DirichletMultinomial_1.36.0
[105] beeswarm_0.4.0              viridisLite_0.4.0          
[107] vipor_0.4.5                 bslib_0.3.1                
```
</div>
