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

```
## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
## logical.return = TRUE, : there is no package called 'ecodist'
```

```
## Installing package into '/__w/_temp/Library'
## (as 'lib' is unspecified)
```


```r
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
```

## Visualizing taxonomic composition

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

<img src="21_microbiome_community_files/figure-html/unnamed-chunk-1-1.png" width="672" />

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
tse_phylum <- transformSamples(tse_phylum, method = "clr", pseudocount = 1)
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

<img src="21_microbiome_community_files/figure-html/heatmap-1.png" width="672" />

_pheatmap_ is a package that provides methods to plot clustered heatmaps. 


```r
if(!require(pheatmap)){
    install.packages("pheatmap")
}
library(pheatmap)

# Takes subset: only samples from feces, skin, or tongue
tse_phylum_subset <- tse_phylum[ , colData(tse_phylum)$SampleType %in% c("Feces", "Skin", "Tongue") ]

# Does clr-transformation
tse_phylum_subset <- transformSamples(tse_phylum_subset, method = "clr", pseudocount = 1)
# Does z-transformation
tse_phylum_subset <- transformFeatures(tse_phylum_subset, abund_values = "clr", 
                                       method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_phylum_subset, top = 20)
tse_phylum_subset <- tse_phylum_subset[top_taxa, ]

# Gets the assay table
mat <- assay(tse_phylum_subset, "clr_z")

# Creates the heatmap
pheatmap(mat)
```

<img src="21_microbiome_community_files/figure-html/pheatmap1-1.png" width="672" />

We can create clusters by hierarchical clustering and 
visualize them with dendrogram.


```r
# Package for creating dendrograms
if(!require(dendextend)){
    install.packages("dendextend")
}
library(dendextend)

# Hierarchical clustering
taxa_clusters <- hclust(dist(mat), method = "complete")

# Creates a dendrogram
taxa_dendrogram <- as.dendrogram(taxa_clusters)

# Plots it
plot(taxa_dendrogram)
```

<img src="21_microbiome_community_files/figure-html/pheatmap2-1.png" width="672" />

Based on dendrogram, we decide to create three clusters. 


```r
# Creates clusters
taxa_clusters <- cutree(tree = taxa_dendrogram, k = 3)

# Prints taxa and their clusters
taxa_clusters 
```

```
##       Firmicutes    Bacteroidetes   Proteobacteria   Actinobacteria 
##                1                1                2                3 
##    Cyanobacteria     Fusobacteria      Tenericutes  Verrucomicrobia 
##                2                2                1                1 
##    Lentisphaerae    Euryarchaeota    Acidobacteria     Spirochaetes 
##                1                1                3                2 
##   Planctomycetes           Thermi      Chloroflexi              SR1 
##                3                3                3                2 
##    Synergistetes    Crenarchaeota Gemmatimonadetes           SAR406 
##                2                3                3                1
```


```r
# Creates clusters, and adds information to rowData
rowData(tse_phylum_subset)$clusters <- cutree(tree = taxa_dendrogram, k = 3)

# Prints taxa and their clusters
rowData(tse_phylum_subset)$clusters
```

```
##       Firmicutes    Bacteroidetes   Proteobacteria   Actinobacteria 
##                1                1                2                3 
##    Cyanobacteria     Fusobacteria      Tenericutes  Verrucomicrobia 
##                2                2                1                1 
##    Lentisphaerae    Euryarchaeota    Acidobacteria     Spirochaetes 
##                1                1                3                2 
##   Planctomycetes           Thermi      Chloroflexi              SR1 
##                3                3                3                2 
##    Synergistetes    Crenarchaeota Gemmatimonadetes           SAR406 
##                2                3                3                1
```

Now we can create heatmap with additional annotations.


```r
# Creates data frame that includes cluster data
taxa_clusters <- rowData(tse_phylum_subset)$clusters
taxa_clusters <- as.character(taxa_clusters)
taxa_clusters <- data.frame(cluster = taxa_clusters)
row.names(taxa_clusters) <- rownames(tse_phylum_subset)

# Creates data frame that includes data on sample type
sample_types <- unfactor(colData(tse_phylum_subset)$SampleType)
sample_types <- data.frame(sample_types = sample_types)
row.names(sample_types) <- colnames(tse_phylum_subset)

pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_types)
```

<img src="21_microbiome_community_files/figure-html/pheatmap5-1.png" width="672" />

In addition to _pheatmap_ package, there are also other packages that
provide functions for more complex heatmaps, such as [_iheatmapr_](https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html).



# Community typing {#community-typing}





### Dirichlet Multinomial Mixtures (DMM)

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
tse_dmn <- runDMN(tse, name = "DMN", k = 1:7)
```


```r
# It is stored in metadata
tse_dmn
```

```
## class: TreeSummarizedExperiment 
## dim: 67 26 
## metadata(1): DMN
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
## [1] "DMN"
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
## Laplace: 7689 BIC: 8076 AIC: 7948 
## 
## [[4]]
## class: DMN 
## k: 4 
## samples x taxa: 26 x 67 
## Laplace: 7751 BIC: 8274 AIC: 8103 
## 
## [[5]]
## class: DMN 
## k: 5 
## samples x taxa: 26 x 67 
## Laplace: 7858 BIC: 8578 AIC: 8364 
## 
## [[6]]
## class: DMN 
## k: 6 
## samples x taxa: 26 x 67 
## Laplace: 7899 BIC: 8753 AIC: 8497 
## 
## [[7]]
## class: DMN 
## k: 7 
## samples x taxa: 26 x 67 
## Laplace: 7991 BIC: 9021 AIC: 8722
```


Show Laplace approximation (model evidence) for each model of the k models.


```r
library(miaViz)
plotDMNFit(tse_dmn, type = "laplace")
```

<img src="21_microbiome_community_files/figure-html/unnamed-chunk-5-1.png" width="672" />

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
## Feces              2       4   67 1078.3 -106.14   901.2 1171.9 1213
## Freshwater         2       2   67  889.6  -97.17   717.0  936.4 1025
## Freshwater (creek) 2       3   67 1600.3  860.08  1906.3 1674.5 1735
## Mock               2       3   67  998.6  -70.60   839.3 1072.8 1134
## Ocean              2       3   67 1096.7  -56.21   944.6 1170.9 1232
## Sediment (estuary) 2       3   67 1195.5   18.63  1080.8 1269.7 1331
## Skin               2       3   67  992.6  -84.81   826.2 1066.8 1128
## Soil               2       3   67 1380.3   11.21  1261.8 1454.5 1515
## Tongue             2       2   67  783.0 -107.74   605.1  829.8  918
```

Mixture weights  (rough measure of the cluster size).



```r
DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))
```

```
##       pi theta
## 1 0.5385 20.58
## 2 0.4615 15.28
```


Samples-cluster assignment probabilities / how probable it is that sample belongs
to each cluster


```r
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
```

```
##              [,1]      [,2]
## CL3     1.000e+00 5.015e-17
## CC1     1.000e+00 3.864e-22
## SV1     1.000e+00 1.946e-12
## M31Fcsw 7.875e-26 1.000e+00
## M11Fcsw 1.133e-16 1.000e+00
## M31Plmr 1.123e-13 1.000e+00
```

Contribution of each taxa to each component


```r
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))
```

```
##                          [,1]      [,2]
## Phylum:Crenarchaeota  0.30381 0.1354661
## Phylum:Euryarchaeota  0.23114 0.1468621
## Phylum:Actinobacteria 1.21364 1.0601012
## Phylum:Spirochaetes   0.21393 0.1318418
## Phylum:MVP-15         0.02982 0.0007658
## Phylum:Proteobacteria 6.84499 1.8153496
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
tse <- transformCounts(tse, method = "clr", pseudocount = 1)

# Gets clr table
clr_assay <- assays(tse)$clr

# Transposes it to get taxa to columns
clr_assay <- t(clr_assay)

# Calculates Euclidean distances between samples. Because taxa is in columns,
# it is used to compare different samples.
euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")

# Does principal coordinate analysis
euclidean_pcoa <- ecodist::pco(euclidean_dist)

# Creates a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                pcoa2 = euclidean_pcoa$vectors[,2])
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

<img src="21_microbiome_community_files/figure-html/unnamed-chunk-13-1.png" width="672" />

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

<img src="21_microbiome_community_files/figure-html/unnamed-chunk-15-1.png" width="672" />

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

<img src="21_microbiome_community_files/figure-html/unnamed-chunk-16-1.png" width="672" /><img src="21_microbiome_community_files/figure-html/unnamed-chunk-16-2.png" width="672" /><img src="21_microbiome_community_files/figure-html/unnamed-chunk-16-3.png" width="672" /><img src="21_microbiome_community_files/figure-html/unnamed-chunk-16-4.png" width="672" />

## Additional Community Typing

For more community typing techniques applied to the 'SprockettTHData' data set, see the attached .Rmd file.

Link:

   * [Rmd](add-comm-typing.Rmd)


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
 [3] patchwork_1.1.1                bluster_1.4.0                 
 [5] dendextend_1.15.2              pheatmap_1.0.12               
 [7] miaViz_1.1.4                   ggraph_2.0.5                  
 [9] ggplot2_3.3.5                  mia_1.3.8                     
[11] MultiAssayExperiment_1.20.0    TreeSummarizedExperiment_2.1.4
[13] Biostrings_2.62.0              XVector_0.34.0                
[15] SingleCellExperiment_1.16.0    SummarizedExperiment_1.24.0   
[17] Biobase_2.54.0                 GenomicRanges_1.46.0          
[19] GenomeInfoDb_1.30.0            IRanges_2.28.0                
[21] S4Vectors_0.32.2               BiocGenerics_0.40.0           
[23] MatrixGenerics_1.6.0           matrixStats_0.61.0-9001       
[25] ecodist_2.0.7                  BiocStyle_2.22.0              
[27] rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] plyr_1.8.6                  igraph_1.2.8               
  [3] lazyeval_0.2.2              splines_4.1.2              
  [5] BiocParallel_1.28.0         digest_0.6.28              
  [7] yulab.utils_0.0.4           htmltools_0.5.2            
  [9] viridis_0.6.2               fansi_0.5.0                
 [11] magrittr_2.0.1              memoise_2.0.0              
 [13] ScaledMatrix_1.2.0          cluster_2.1.2              
 [15] DECIPHER_2.22.0             graphlayouts_0.7.1         
 [17] colorspace_2.0-2            blob_1.2.2                 
 [19] ggrepel_0.9.1               xfun_0.28                  
 [21] dplyr_1.0.7                 crayon_1.4.2               
 [23] RCurl_1.98-1.5              jsonlite_1.7.2             
 [25] graph_1.72.0                ape_5.5                    
 [27] glue_1.5.0                  polyclip_1.10-0            
 [29] gtable_0.3.0                zlibbioc_1.40.0            
 [31] DelayedArray_0.20.0         BiocSingular_1.10.0        
 [33] scales_1.1.1                DBI_1.1.1                  
 [35] Rcpp_1.0.7                  viridisLite_0.4.0          
 [37] decontam_1.14.0             gridGraphics_0.5-1         
 [39] tidytree_0.3.6              bit_4.0.4                  
 [41] rsvd_1.0.5                  FNN_1.1.3                  
 [43] RColorBrewer_1.1-2          dir.expiry_1.2.0           
 [45] ellipsis_0.3.2              pkgconfig_2.0.3            
 [47] XML_3.99-0.8                farver_2.1.0               
 [49] uwot_0.1.10                 CodeDepends_0.6.5          
 [51] sass_0.4.0                  utf8_1.2.2                 
 [53] ggplotify_0.1.0             tidyselect_1.1.1           
 [55] labeling_0.4.2              rlang_0.4.12               
 [57] reshape2_1.4.4              munsell_0.5.0              
 [59] tools_4.1.2                 cachem_1.0.6               
 [61] DirichletMultinomial_1.36.0 generics_0.1.1             
 [63] RSQLite_2.2.8               evaluate_0.14              
 [65] stringr_1.4.0               fastmap_1.1.0              
 [67] yaml_2.2.1                  ggtree_3.2.1               
 [69] knitr_1.36                  bit64_4.0.5                
 [71] tidygraph_1.2.0             purrr_0.3.4                
 [73] nlme_3.1-153                sparseMatrixStats_1.6.0    
 [75] aplot_0.1.1                 compiler_4.1.2             
 [77] beeswarm_0.4.0              filelock_1.0.2             
 [79] treeio_1.18.1               tibble_3.1.6               
 [81] tweenr_1.0.2                bslib_0.3.1                
 [83] stringi_1.7.5               highr_0.9                  
 [85] lattice_0.20-45             Matrix_1.3-4               
 [87] vegan_2.5-7                 permute_0.9-5              
 [89] vctrs_0.3.8                 pillar_1.6.4               
 [91] lifecycle_1.0.1             BiocManager_1.30.16        
 [93] jquerylib_0.1.4             BiocNeighbors_1.12.0       
 [95] cowplot_1.1.1               bitops_1.0-7               
 [97] irlba_2.3.3                 R6_2.5.1                   
 [99] bookdown_0.24               gridExtra_2.3              
[101] vipor_0.4.5                 codetools_0.2-18           
[103] MASS_7.3-54                 assertthat_0.2.1           
[105] withr_2.4.2                 GenomeInfoDbData_1.2.7     
[107] mgcv_1.8-38                 parallel_4.1.2             
[109] grid_4.1.2                  ggfun_0.0.4                
[111] beachmat_2.10.0             tidyr_1.1.4                
[113] rmarkdown_2.11              DelayedMatrixStats_1.16.0  
[115] ggnewscale_0.4.5            ggforce_0.3.3              
[117] ggbeeswarm_0.6.0           
```
</div>

