# Community composition {#microbiome-community}

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

```
## Loading required package: ecodist
```


```r
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
```

## Visualizing taxonomic composition {#visual-composition}

### Composition barplot

A typical way to visualize microbiome composition is by using
composition barplot. In the following, relative abundance is
calculated and top taxa are retrieved for the Phylum rank. Thereafter,
the barplot is visualized ordering rank by abundance values and
samples by "Bacteroidetes":


```r
library(miaViz)
# Computing relative abundance
tse <- relAbundanceCounts(tse)

# Getting top taxa on a Phylum level
tse_phylum <- agglomerateByRank(tse, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopTaxa(tse_phylum,top = 5, abund_values = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$Phylum,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse)$Phylum <- as.character(phylum_renamed)

# Visualizing the composition barplot, with samples order by "Bacteroidetes"
plotAbundance(tse, abund_values="relabundance", rank = "Phylum",
              order_rank_by="abund", order_sample_by = "Bacteroidetes")
```

![](21_microbiome_community_files/figure-latex/unnamed-chunk-1-1.pdf)<!-- --> 

### Composition heatmap 

Community composition can be visualized with heatmap, where the
horizontal axis represents samples and the vertical axis the
taxa. Color of each intersection point represents abundance of a taxon
in a specific sample.

Here,  abundances are  first CLR  (centered log-ratio)  transformed to
remove  compositionality bias. Then  Z  transformation  is applied  to
CLR-transformed  data. This  shifts all  taxa  to zero  mean and  unit
variance, allowing visual comparison  between taxa that have different
absolute  abundance  levels.  After  these  rough  visual  exploration
techniques, we can visualize the abundances at Phylum level.


```r
library(ggplot2)
# Add clr-transformation on samples
tse_phylum <- transformSamples(tse_phylum, method = "relabundance", pseudocount = 1)
tse_phylum <- transformSamples(tse_phylum, abund_values = "relabundance", method = "clr")
# Add z-transformation on features (taxa)
tse_phylum <- transformFeatures(tse_phylum, abund_values = "clr", 
                                method = "z", name = "clr_z")
# Melts the assay
df <- meltAssay(tse_phylum, abund_values = "clr_z")

# Determines the scaling of colours
maxval <- round(max(abs(df$clr_z)))
limits <- c(-maxval, maxval)
breaks <- seq(from = min(limits), to = max(limits), by = 0.5)
colours <- c("darkblue", "blue", "white", "red", "darkred")

# Creates a ggplot object
ggplot(df, aes(x = SampleID, y = FeatureID, fill = clr_z)) +
  geom_tile() +
  scale_fill_gradientn(name = "CLR + Z transform", 
                       breaks = breaks, limits = limits, colours = colours) + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.size = unit(1, "cm")) +
  labs(x = "Samples", y = "Taxa")
```

![](21_microbiome_community_files/figure-latex/heatmap-1.pdf)<!-- --> 

In addition, there are also other packages that provide functions for more complex heatmaps,
such as [_iheatmapr_](https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html)
and [ComplexHeatmap](https://academic.oup.com/bioinformatics/article/32/18/2847/1743594?login=true).
[sechm](http://www.bioconductor.org/packages/release/bioc/vignettes/sechm/inst/doc/sechm.html)
package provides wrapper for _ComplexHeatmap_ and its usage is explained in chapter \@ref(viz-chapter)
along with the `pheatmap` package for clustered heatmaps.

# Community typing {#community-typing}





## Dirichlet Multinomial Mixtures (DMM)

This section focus on DMM analysis. 

One technique that allows to search for groups of samples that are
similar to each other is the [Dirichlet-Multinomial Mixture
Model](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030126). In
DMM, we first determine the number of clusters (k) that best fit the
data (model evidence) using Laplace approximation. After fitting the
model with k clusters, we obtain for each sample k probabilities that
reflect the probability that a sample belongs to the given cluster.

Let's cluster the data with DMM clustering. 


```r
# Runs model and calculates the most likely number of clusters from 1 to 7.
# Since this is a large dataset it takes long computational time.
# For this reason we use only a subset of the data; agglomerated by Phylum as a rank.
tse <- GlobalPatterns
tse <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree=TRUE)
```


```r
tse_dmn <- mia::runDMN(tse, name = "DMN", k = 1:7)
```


```r
# It is stored in metadata
tse_dmn
```

```
## class: TreeSummarizedExperiment 
## dim: 67 26 
## metadata(2): agglomerated_by_rank DMN
## assays(1): counts
## rownames(67): Phylum:Crenarchaeota Phylum:Euryarchaeota ...
##   Phylum:Synergistetes Phylum:SR1
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(7): X.SampleID Primer ... SampleType Description
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (67 rows)
## rowTree: 1 phylo tree(s) (66 leaves)
## colLinks: NULL
## colTree: NULL
```

Return information on metadata that the object contains.


```r
names(metadata(tse_dmn))
```

```
## [1] "agglomerated_by_rank" "DMN"
```

This returns a list of DMN objects for a closer investigation.


```r
getDMN(tse_dmn)
```

```
## [[1]]
## class: DMN 
## k: 1 
## samples x taxa: 26 x 67 
## Laplace: 7715 BIC: 7802 AIC: 7760 
## 
## [[2]]
## class: DMN 
## k: 2 
## samples x taxa: 26 x 67 
## Laplace: 7673 BIC: 7927 AIC: 7842 
## 
## [[3]]
## class: DMN 
## k: 3 
## samples x taxa: 26 x 67 
## Laplace: 7690 BIC: 8076 AIC: 7948 
## 
## [[4]]
## class: DMN 
## k: 4 
## samples x taxa: 26 x 67 
## Laplace: 7741 BIC: 8282 AIC: 8112 
## 
## [[5]]
## class: DMN 
## k: 5 
## samples x taxa: 26 x 67 
## Laplace: 7857 BIC: 8578 AIC: 8364 
## 
## [[6]]
## class: DMN 
## k: 6 
## samples x taxa: 26 x 67 
## Laplace: 7942 BIC: 8822 AIC: 8566 
## 
## [[7]]
## class: DMN 
## k: 7 
## samples x taxa: 26 x 67 
## Laplace: 8039 BIC: 9073 AIC: 8775
```


Show Laplace approximation (model evidence) for each model of the k models.


```r
library(miaViz)
plotDMNFit(tse_dmn, type = "laplace")
```

![](21_microbiome_community_files/figure-latex/unnamed-chunk-5-1.pdf)<!-- --> 

Return the model that has the best fit.


```r
getBestDMNFit(tse_dmn, type = "laplace")
```

```
## class: DMN 
## k: 2 
## samples x taxa: 26 x 67 
## Laplace: 7673 BIC: 7927 AIC: 7842
```

### PCoA for ASV-level data with Bray-Curtis; with DMM clusters shown with colors

Group samples and return DMNGroup object that contains a summary.
Patient status is used for grouping.


```r
dmn_group <- calculateDMNgroup(tse_dmn, variable = "SampleType",  exprs_values = "counts",
                               k = 2, seed=.Machine$integer.max)

dmn_group
```

```
## class: DMNGroup 
## summary:
##                    k samples taxa    NLE  LogDet Laplace    BIC  AIC
## Feces              2       4   67 1078.3 -106.26   901.1 1171.9 1213
## Freshwater         2       2   67  889.6  -97.20   716.9  936.4 1025
## Freshwater (creek) 2       3   67 1600.3  793.17  1872.8 1674.5 1735
## Mock               2       3   67  998.6  -70.65   839.2 1072.8 1134
## Ocean              2       3   67 1096.7  -56.66   944.3 1170.9 1232
## Sediment (estuary) 2       3   67 1195.5   18.63  1080.8 1269.7 1331
## Skin               2       3   67  992.6  -85.05   826.1 1066.8 1128
## Soil               2       3   67 1380.3   11.20  1261.8 1454.5 1515
## Tongue             2       2   67  783.0 -107.79   605.0  829.8  918
```

Mixture weights  (rough measure of the cluster size).



```r
DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
```

```
##       pi theta
## 1 0.5385 20.59
## 2 0.4615 15.32
```


Samples-cluster assignment probabilities / how probable it is that sample belongs
to each cluster


```r
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
```

```
##              [,1]      [,2]
## CL3     1.000e+00 4.444e-17
## CC1     1.000e+00 3.348e-22
## SV1     1.000e+00 1.703e-12
## M31Fcsw 7.367e-26 1.000e+00
## M11Fcsw 1.086e-16 1.000e+00
## M31Plmr 1.150e-13 1.000e+00
```

Contribution of each taxa to each component


```r
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))
```

```
##                          [,1]      [,2]
## Phylum:Crenarchaeota  0.30381 0.1354040
## Phylum:Euryarchaeota  0.23114 0.1468866
## Phylum:Actinobacteria 1.21370 1.0581361
## Phylum:Spirochaetes   0.21393 0.1318059
## Phylum:MVP-15         0.02983 0.0007634
## Phylum:Proteobacteria 6.84617 1.8113341
```
Get the assignment probabilities



```r
prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
# Add column names
colnames(prob) <- c("comp1", "comp2")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec <- colnames(prob)[max.col(prob,ties.method = "first")]
```

Computing the euclidean PCoA and storing it as a data frame


```r
# Does clr transformation. Pseudocount is added, because data contains zeros.
tse <- transformCounts(tse, method = "relabundance", pseudocount = 1)
tse <- transformCounts(tse, "relabundance", method = "clr")

library(scater)

# Does principal coordinate analysis
df <- calculateMDS(tse, exprs_values = "clr", method = "euclidean")

# Creates a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = df[,1], 
                                pcoa2 = df[,2])
```


```r
# Creates a data frame that contains principal coordinates and DMM information
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = vec)
# Creates a plot
euclidean_dmm_plot <- ggplot(data = euclidean_dmm_pcoa_df, 
                             aes(x=pcoa1, y=pcoa2,
                                 color = dmm_component)) +
  geom_point() +
  labs(x = "Coordinate 1",
       y = "Coordinate 2",
       title = "PCoA with Aitchison distances") +  
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_dmm_plot
```

![](21_microbiome_community_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> 

## Community Detection

Another approach for discovering communities within the samples of the
data, is to run community detection algorithms after building a
graph. The following demonstration builds a graph based on the k
nearest-neighbors and performs the community detection on the fly.

_`bluster`_ [@R-bluster] package offers several clustering methods,
among which graph-based are present, enabling the community detection
task.

Installing package:


```r
if(!require(bluster)){
  BiocManager::install("bluster")
}
```

The algorithm used is "short random walks" [@Pons2006]. Graph is
constructed using different k values (the number of nearest neighbors
to consider during graph construction) using the robust centered log
ratio (rclr) assay data. Then plotting the communities using UMAP
[@McInnes2018] ordination as a visual exploration aid.  In the
following demonstration we use the `enterotype` dataset from the
[@R-mia] package.


```r
library(bluster)
library(patchwork) # For arranging several plots as a grid
library(scater)

data("enterotype", package="mia")
tse <- enterotype
tse <- transformCounts(tse, method = "rclr")

# Performing and storing UMAP
tse <- runUMAP(tse, name="UMAP", exprs_values="rclr")

k <- c(2,3,5,10)
ClustAndPlot <- function(x) {
  # Creating the graph and running the short random walks algorithm  
  graph_clusters <- clusterRows(t(assays(tse)$rclr), NNGraphParam(k=x))
  
  # Results of the clustering as a color for each sample
  plotUMAP(tse, colour_by = I(graph_clusters)) +
    labs(title = paste0("k = ", x))
}

# Applying the function for different k values
plots <- lapply(k,ClustAndPlot)

# Displaying plots in a grid
(plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]])
```

![](21_microbiome_community_files/figure-latex/unnamed-chunk-15-1.pdf)<!-- --> 

Similarly, the _`bluster`_ [@R-bluster] package offers clustering
diagnostics that can be used for judging the clustering quality (see
[Assorted clustering
diagnostics](http://bioconductor.org/packages/release/bioc/vignettes/bluster/inst/doc/diagnostics.html)).
In the following, Silhouette width as a diagnostic tool is computed
and results are visualized for each case presented earlier. For more
about Silhouettes read [@Rousseeuw1987].


```r
ClustDiagPlot <- function(x) {
  # Getting the clustering results
  graph_clusters <- clusterRows(t(assays(tse)$rclr), NNGraphParam(k=x))
  
  # Computing the diagnostic info
  sil <- approxSilhouette(t(assays(tse)$rclr), graph_clusters)
  
  # Plotting as a boxlpot to observe cluster separation
  boxplot(split(sil$width, graph_clusters), main=paste0("k = ", x))
  
}
# Applying the function for different k values
res <- lapply(k,ClustDiagPlot)
```

![](21_microbiome_community_files/figure-latex/unnamed-chunk-16-1.pdf)<!-- --> ![](21_microbiome_community_files/figure-latex/unnamed-chunk-16-2.pdf)<!-- --> ![](21_microbiome_community_files/figure-latex/unnamed-chunk-16-3.pdf)<!-- --> ![](21_microbiome_community_files/figure-latex/unnamed-chunk-16-4.pdf)<!-- --> 

## Additional Community Typing

For more community typing techniques applied to the 'SprockettTHData' data set, see the attached .Rmd file.

Link:

   * [Rmd](add-comm-typing.Rmd)


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.2.0 (2022-04-22)
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
 [1] patchwork_1.1.1                bluster_1.6.0                 
 [3] scater_1.24.0                  scuttle_1.6.2                 
 [5] miaViz_1.3.3                   ggraph_2.0.5                  
 [7] ggplot2_3.3.6                  mia_1.3.26                    
 [9] MultiAssayExperiment_1.22.0    TreeSummarizedExperiment_2.1.4
[11] Biostrings_2.64.0              XVector_0.36.0                
[13] SingleCellExperiment_1.18.0    SummarizedExperiment_1.26.1   
[15] Biobase_2.56.0                 GenomicRanges_1.48.0          
[17] GenomeInfoDb_1.32.2            IRanges_2.30.0                
[19] S4Vectors_0.34.0               BiocGenerics_0.42.0           
[21] MatrixGenerics_1.8.0           matrixStats_0.62.0-9000       
[23] ecodist_2.0.9                  BiocStyle_2.24.0              
[25] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] plyr_1.8.7                  igraph_1.3.2               
  [3] lazyeval_0.2.2              splines_4.2.0              
  [5] BiocParallel_1.30.3         digest_0.6.29              
  [7] yulab.utils_0.0.4           htmltools_0.5.2            
  [9] viridis_0.6.2               fansi_1.0.3                
 [11] magrittr_2.0.3              memoise_2.0.1              
 [13] ScaledMatrix_1.4.0          cluster_2.1.3              
 [15] DECIPHER_2.24.0             graphlayouts_0.8.0         
 [17] colorspace_2.0-3            blob_1.2.3                 
 [19] ggrepel_0.9.1               xfun_0.31                  
 [21] dplyr_1.0.9                 crayon_1.5.1               
 [23] RCurl_1.98-1.7              jsonlite_1.8.0             
 [25] graph_1.74.0                ape_5.6-2                  
 [27] glue_1.6.2                  polyclip_1.10-0            
 [29] gtable_0.3.0                zlibbioc_1.42.0            
 [31] DelayedArray_0.22.0         BiocSingular_1.12.0        
 [33] scales_1.2.0                DBI_1.1.3                  
 [35] Rcpp_1.0.8.3                viridisLite_0.4.0          
 [37] decontam_1.16.0             gridGraphics_0.5-1         
 [39] tidytree_0.3.9              bit_4.0.4                  
 [41] rsvd_1.0.5                  FNN_1.1.3.1                
 [43] dir.expiry_1.4.0            ellipsis_0.3.2             
 [45] pkgconfig_2.0.3             XML_3.99-0.10              
 [47] farver_2.1.0                uwot_0.1.11                
 [49] CodeDepends_0.6.5           utf8_1.2.2                 
 [51] ggplotify_0.1.0             tidyselect_1.1.2           
 [53] labeling_0.4.2              rlang_1.0.2                
 [55] reshape2_1.4.4              munsell_0.5.0              
 [57] tools_4.2.0                 cachem_1.0.6               
 [59] cli_3.3.0                   DirichletMultinomial_1.38.0
 [61] generics_0.1.2              RSQLite_2.2.14             
 [63] evaluate_0.15               stringr_1.4.0              
 [65] fastmap_1.1.0               yaml_2.3.5                 
 [67] ggtree_3.4.0                knitr_1.39                 
 [69] bit64_4.0.5                 tidygraph_1.2.1            
 [71] purrr_0.3.4                 nlme_3.1-158               
 [73] sparseMatrixStats_1.8.0     aplot_0.1.6                
 [75] compiler_4.2.0              beeswarm_0.4.0             
 [77] filelock_1.0.2              treeio_1.20.0              
 [79] tibble_3.1.7                tweenr_1.0.2               
 [81] stringi_1.7.6               highr_0.9                  
 [83] lattice_0.20-45             Matrix_1.4-1               
 [85] vegan_2.6-2                 permute_0.9-7              
 [87] vctrs_0.4.1                 pillar_1.7.0               
 [89] lifecycle_1.0.1             BiocManager_1.30.18        
 [91] BiocNeighbors_1.14.0        cowplot_1.1.1              
 [93] bitops_1.0-7                irlba_2.3.5                
 [95] R6_2.5.1                    bookdown_0.27              
 [97] gridExtra_2.3               vipor_0.4.5                
 [99] codetools_0.2-18            MASS_7.3-57                
[101] assertthat_0.2.1            withr_2.5.0                
[103] GenomeInfoDbData_1.2.8      mgcv_1.8-40                
[105] parallel_4.2.0              grid_4.2.0                 
[107] ggfun_0.0.6                 beachmat_2.12.0            
[109] tidyr_1.2.0                 rmarkdown_2.14             
[111] DelayedMatrixStats_1.18.0   ggnewscale_0.4.7           
[113] ggforce_0.3.3               ggbeeswarm_0.6.0           
```
</div>

