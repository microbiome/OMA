# Microbiome Community {#microbiome-community}

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

## Community composition

### Composition barplot

A typical way to visualize microbiome composition is by using composition barplot.
In the following, relative abundance is calculated and top taxa are retrieved for the
Phylum rank. Thereafter, the barplot is visualized ordering rank by abundance values
and samples by "Bacteroidetes":


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

<img src="16_microbiome_community_files/figure-html/unnamed-chunk-1-1.png" width="672" />

### Composition heatmap 

Community composition can be visualized with heatmap, where the horizontal axis represents
samples and the vertical axis the taxa. Color of each intersection point represents abundance
of a taxon in a specific sample. 

Here, abundances are first CLR (centered log ratio) transformed, and then 
Z transformation is applied to CLR-transformed data. After that, abundances are 
plotted at Phylum level. 


```r
library(ggplot2)
# Does clr-transformation
tse_phylum <- transformSamples(tse_phylum, method = "clr", pseudocount = 1)
# Does z-transformation
tse_phylum <- transformFeatures(tse_phylum, abund_values = "clr", 
                                method = "z", name = "clr_z")
# Melts the assay
df <- meltAssay(tse_phylum, assay_name = "clr_z")

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

<img src="16_microbiome_community_files/figure-html/heatmap-1.png" width="672" />

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

<img src="16_microbiome_community_files/figure-html/pheatmap1-1.png" width="672" />

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

<img src="16_microbiome_community_files/figure-html/pheatmap2-1.png" width="672" />

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

<img src="16_microbiome_community_files/figure-html/pheatmap5-1.png" width="672" />

In addition to _pheatmap_ package, there are also other packages that provide 
functions for more complex heatmaps. One example is _iheatmapr_ package. Examples
of using it you can find from its 
[vignette](https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html).

# Cross-correlation

With cross-correlation analysis, we can analyze how strongly and how differently
variables are associated between each other. For instance, we can analyze if 
higher presence of a specific taxon equals to higher levels of a biomolecule. 

Here, we analyze associations between taxa correlate and lipids. Data is from 
following publication Lahti _et al_. (2015) [Associations between the human intestinal 
microbiota, Lactobacillus rhamnosus GG and serum lipids indicated by 
integrated analysis of high-throughput profiling 
data](https://peerj.com/articles/32/).


```r
# Imports the data
tse <- microbiomeDataSets::peerj32()

# Microbiome data
tse[[1]] 
```

```
## class: TreeSummarizedExperiment 
## dim: 130 44 
## metadata(0):
## assays(1): counts
## rownames(130): Actinomycetaceae Aerococcus ... Xanthomonadaceae
##   Yersinia et rel.
## rowData names(3): Phylum Family Genus
## colnames(44): sample-1 sample-2 ... sample-43 sample-44
## colData names(0):
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```


```r
# Lipid data
tse[[2]]
```

```
## class: SummarizedExperiment 
## dim: 389 44 
## metadata(0):
## assays(1): counts
## rownames(389): Cer(d18:1/16:0).1 Cer(d18:1/16:0).2 ... TG(60:11)
##   TG(60:9)
## rowData names(0):
## colnames(44): sample-1 sample-2 ... sample-43 sample-44
## colData names(0):
```


```r
if(!require(microbiome)){
    BiocManager::install("microbiome")
}

# Does log10 transform for microbiome data
tse[[1]] <- transformSamples(tse[[1]], method = "log10", pseudocount = 1)

# Gets microbiome and lipid data to cross-correlate
x <- t(assay(tse[[1]], "log10"))
y <- t(assay(tse[[2]], "counts"))

# Cross correlates data sets
correlation_table <- microbiome::associate(x, y, method = "spearman", mode = "table", 
                                           p.adj.threshold = 0.05, n.signif = 1)

knitr::kable(head(correlation_table))
```



|    |X1                              |X2         | Correlation|  p.adj|
|:---|:-------------------------------|:----------|-----------:|------:|
|552 |Ruminococcus gnavus et rel.     |TG(54:5).2 |      0.7165| 0.0023|
|614 |Uncultured Bacteroidetes        |TG(56:2).1 |     -0.6964| 0.0039|
|100 |Lactobacillus plantarum et rel. |PC(40:3)   |     -0.6737| 0.0051|
|252 |Ruminococcus gnavus et rel.     |TG(50:4)   |      0.6912| 0.0051|
|357 |Ruminococcus gnavus et rel.     |TG(52:5)   |      0.6806| 0.0051|
|537 |Ruminococcus gnavus et rel.     |TG(54:4).2 |      0.6820| 0.0051|

Manipulates and reorders the table


```r
if(!require(reshape2)){
    install.packages("reshape2")
}
if(!require(dplyr)){
    install.packages("dplyr")
}
# Gets taxa that have at least one statistically significant correlation
taxa <- correlation_table %>% dplyr::filter(p.adj < 0.05 & abs(Correlation) > 0) %>% 
  dplyr::select(X1) %>% unique %>% sapply(as.character)
# Gets taxa that have at least one statistically significant correlation
lipids <- correlation_table %>% dplyr::filter(p.adj < 0.05 & abs(Correlation) > 0) %>% 
  dplyr::select(X2) %>% unique %>% sapply(as.character)

# Takes only those taxa and lipids that have statistically significant values
correlation_table <- correlation_table[correlation_table[["X1"]] %in% taxa & 
                                           correlation_table[["X2"]] %in% lipids, ]

# Converts data to matrix, correlations as values
mat <- reshape2::acast(correlation_table, X2 ~ X1, value.var = "Correlation")

# Hierarchical clustering, gets the order of taxa and lipids
taxa_indices <- hclust(as.dist(1 - cor(mat, use="pairwise.complete.obs")))$order
order_taxa <- colnames(mat)[taxa_indices]
lipids_indices <- hclust(as.dist(1 - cor(t(mat), use="pairwise.complete.obs")))$order
order_lipids <- rownames(mat)[lipids_indices]

# Converts taxa and lipids columns to factor so that they have the desired order
correlation_table[["X1"]] <- factor(correlation_table[["X1"]], levels = order_taxa)
correlation_table[["X2"]] <- factor(correlation_table[["X2"]], levels = order_lipids)
```

Creates the heatmap


```r
# Determines the scaling of colors
limits <- c(-1, 1)
breaks <- seq(from = min(limits), to = max(limits), by = 0.2)
colours <- c("darkblue", "blue", "white", "red", "darkred")

# Which observation have p-value under 0.05? --> creates a subset
cor_table_sub <- correlation_table[which(correlation_table[["Correlation"]] < 0.05), ]

# Creates a ggplot object
ggplot(correlation_table, aes(x = X1, y = X2, fill = Correlation)) +
  geom_tile() +
  
  scale_fill_gradientn(name = "Correlation", 
                       breaks = breaks, limits = limits, colours = colours) + 
  
  # Adds label to those observations that have p-value under 0.05
  geom_text(data = cor_table_sub, aes(x = X1, y = X2, label = "+")) +
  
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.size = unit(1, "cm")) +
  labs(x = "Taxa", y = "Lipids")
```

<img src="16_microbiome_community_files/figure-html/cross-correlation5-1.png" width="576" />

## Community typing

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
## Laplace: 7750 BIC: 8154 AIC: 8026 
## 
## [[4]]
## class: DMN 
## k: 4 
## samples x taxa: 26 x 67 
## Laplace: 7753 BIC: 8276 AIC: 8105 
## 
## [[5]]
## class: DMN 
## k: 5 
## samples x taxa: 26 x 67 
## Laplace: 7855 BIC: 8553 AIC: 8340 
## 
## [[6]]
## class: DMN 
## k: 6 
## samples x taxa: 26 x 67 
## Laplace: 7986 BIC: 8881 AIC: 8625 
## 
## [[7]]
## class: DMN 
## k: 7 
## samples x taxa: 26 x 67 
## Laplace: 8017 BIC: 9020 AIC: 8721
```


Show Laplace approximation (model evidence) for each model of the k models.


```r
library(miaViz)
plotDMNFit(tse_dmn, type = "laplace")
```

<img src="16_microbiome_community_files/figure-html/unnamed-chunk-5-1.png" width="672" />

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
## Feces              2       4   67 1078.3 -106.19   901.1 1171.9 1213
## Freshwater         2       2   67  889.6  -97.28   716.9  936.4 1025
## Freshwater (creek) 2       3   67 1600.3  860.08  1906.3 1674.5 1735
## Mock               2       3   67 1008.4  -55.37   856.6 1082.5 1143
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
## 1 0.5385 20.60
## 2 0.4615 15.28
```


Samples-cluster assignment probabilities / how probable it is that sample belongs
to each cluster


```r
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))
```

```
##              [,1]      [,2]
## CL3     1.000e+00 4.963e-17
## CC1     1.000e+00 3.749e-22
## SV1     1.000e+00 2.012e-12
## M31Fcsw 7.318e-26 1.000e+00
## M11Fcsw 1.062e-16 1.000e+00
## M31Plmr 9.985e-14 1.000e+00
```

Contribution of each taxa to each component


```r
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))
```

```
##                          [,1]      [,2]
## Phylum:Crenarchaeota  0.30431 0.1354654
## Phylum:Euryarchaeota  0.23143 0.1468629
## Phylum:Actinobacteria 1.21041 1.0600130
## Phylum:Spirochaetes   0.21410 0.1318414
## Phylum:MVP-15         0.02991 0.0007627
## Phylum:Proteobacteria 6.84236 1.8155208
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

<img src="16_microbiome_community_files/figure-html/unnamed-chunk-13-1.png" width="672" />

## Community Detection

Another approach for discovering communities within the samples of the data under question,
is to run community detection algorithms after building a graph. The following demonstration 
builds a graph based on the k nearest-neighbors and performs the community detection on the fly.

_`bluster`_ [@R-bluster] package offers several clustering methods, among which graph-based are
present, enabling the community detection task.

Installing package:


```r
if(!require(bluster)){
  BiocManager::install("bluster")
}
```

The algorithm used is "short random walks" [@Pons2006]. 
Graph is constructed using different k values (the number of nearest neighbors to consider during graph construction) 
using the robust centered log ratio (rclr) assay data. Then plotting the communities using UMAP [@McInnes2018] ordination as a visual exploration aid.
In the following demonstration we use the `enterotype` dataset from the [@R-mia] package.


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

<img src="16_microbiome_community_files/figure-html/unnamed-chunk-15-1.png" width="672" />

Similarly, the _`bluster`_ [@R-bluster] package offers clustering diagnostics
that can be used for judging the clustering quality (see [Assorted clustering diagnostics](http://bioconductor.org/packages/release/bioc/vignettes/bluster/inst/doc/diagnostics.html)).
In the following, Silhouette width as a diagnostic tool is computed and results are visualized
for each case presented earlier. For more about Silhouettes read [@Rousseeuw1987].


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

<img src="16_microbiome_community_files/figure-html/unnamed-chunk-16-1.png" width="672" /><img src="16_microbiome_community_files/figure-html/unnamed-chunk-16-2.png" width="672" /><img src="16_microbiome_community_files/figure-html/unnamed-chunk-16-3.png" width="672" /><img src="16_microbiome_community_files/figure-html/unnamed-chunk-16-4.png" width="672" />

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
 [1] scater_1.21.3                  scuttle_1.3.1                 
 [3] patchwork_1.1.1                bluster_1.3.0                 
 [5] dplyr_1.0.7                    reshape2_1.4.4                
 [7] microbiome_1.15.0              phyloseq_1.37.0               
 [9] microbiomeDataSets_1.1.1       MultiAssayExperiment_1.19.5   
[11] dendextend_1.15.1              pheatmap_1.0.12               
[13] miaViz_1.1.4                   ggraph_2.0.5                  
[15] ggplot2_3.3.5                  mia_1.1.9                     
[17] TreeSummarizedExperiment_2.1.4 Biostrings_2.61.2             
[19] XVector_0.33.0                 SingleCellExperiment_1.15.1   
[21] SummarizedExperiment_1.23.1    Biobase_2.53.0                
[23] GenomicRanges_1.45.0           GenomeInfoDb_1.29.3           
[25] IRanges_2.27.0                 S4Vectors_0.31.0              
[27] BiocGenerics_0.39.1            MatrixGenerics_1.5.3          
[29] matrixStats_0.60.0             ecodist_2.0.7                 
[31] BiocStyle_2.21.3               rebook_1.3.0                  

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                    tidyselect_1.1.1             
  [3] RSQLite_2.2.7                 AnnotationDbi_1.55.1         
  [5] grid_4.1.0                    BiocParallel_1.27.3          
  [7] Rtsne_0.15                    munsell_0.5.0                
  [9] ScaledMatrix_1.1.0            codetools_0.2-18             
 [11] withr_2.4.2                   colorspace_2.0-2             
 [13] filelock_1.0.2                highr_0.9                    
 [15] knitr_1.33                    labeling_0.4.2               
 [17] GenomeInfoDbData_1.2.6        polyclip_1.10-0              
 [19] bit64_4.0.5                   farver_2.1.0                 
 [21] rhdf5_2.37.0                  vctrs_0.3.8                  
 [23] treeio_1.17.2                 generics_0.1.0               
 [25] xfun_0.25                     BiocFileCache_2.1.1          
 [27] R6_2.5.0                      ggbeeswarm_0.6.0             
 [29] graphlayouts_0.7.1            rsvd_1.0.5                   
 [31] rhdf5filters_1.5.0            bitops_1.0-7                 
 [33] cachem_1.0.5                  DelayedArray_0.19.1          
 [35] assertthat_0.2.1              promises_1.2.0.1             
 [37] scales_1.1.1                  beeswarm_0.4.0               
 [39] gtable_0.3.0                  beachmat_2.9.1               
 [41] tidygraph_1.2.0               rlang_0.4.11                 
 [43] splines_4.1.0                 lazyeval_0.2.2               
 [45] BiocManager_1.30.16           yaml_2.2.1                   
 [47] httpuv_1.6.1                  tools_4.1.0                  
 [49] bookdown_0.22                 ellipsis_0.3.2               
 [51] decontam_1.13.0               jquerylib_0.1.4              
 [53] biomformat_1.21.0             RColorBrewer_1.1-2           
 [55] Rcpp_1.0.7                    plyr_1.8.6                   
 [57] sparseMatrixStats_1.5.2       zlibbioc_1.39.0              
 [59] purrr_0.3.4                   RCurl_1.98-1.3               
 [61] viridis_0.6.1                 cowplot_1.1.1                
 [63] ggrepel_0.9.1                 cluster_2.1.2                
 [65] DECIPHER_2.21.0               magrittr_2.0.1               
 [67] data.table_1.14.0             ggnewscale_0.4.5             
 [69] mime_0.11                     evaluate_0.14                
 [71] xtable_1.8-4                  XML_3.99-0.6                 
 [73] gridExtra_2.3                 compiler_4.1.0               
 [75] tibble_3.1.3                  crayon_1.4.1                 
 [77] htmltools_0.5.1.1             mgcv_1.8-36                  
 [79] later_1.2.0                   tidyr_1.1.3                  
 [81] aplot_0.0.6                   DBI_1.1.1                    
 [83] tweenr_1.0.2                  ExperimentHub_2.1.4          
 [85] dbplyr_2.1.1                  MASS_7.3-54                  
 [87] rappdirs_0.3.3                Matrix_1.3-4                 
 [89] ade4_1.7-17                   permute_0.9-5                
 [91] parallel_4.1.0                igraph_1.2.6                 
 [93] pkgconfig_2.0.3               rvcheck_0.1.8                
 [95] dir.expiry_1.1.0              foreach_1.5.1                
 [97] ggtree_3.1.3                  vipor_0.4.5                  
 [99] bslib_0.2.5.1                 DirichletMultinomial_1.35.0  
[101] multtest_2.49.0               stringr_1.4.0                
[103] digest_0.6.27                 vegan_2.5-7                  
[105] graph_1.71.2                  rmarkdown_2.10               
[107] tidytree_0.3.4                uwot_0.1.10                  
[109] DelayedMatrixStats_1.15.2     curl_4.3.2                   
[111] shiny_1.6.0                   lifecycle_1.0.0              
[113] nlme_3.1-152                  jsonlite_1.7.2               
[115] Rhdf5lib_1.15.2               BiocNeighbors_1.11.0         
[117] CodeDepends_0.6.5             viridisLite_0.4.0            
[119] fansi_0.5.0                   pillar_1.6.2                 
[121] lattice_0.20-44               survival_3.2-11              
[123] KEGGREST_1.33.0               fastmap_1.1.0                
[125] httr_1.4.2                    interactiveDisplayBase_1.31.2
[127] glue_1.4.2                    FNN_1.1.3                    
[129] png_0.1-7                     iterators_1.0.13             
[131] BiocVersion_3.14.0            bit_4.0.4                    
[133] ggforce_0.3.3                 stringi_1.7.3                
[135] sass_0.4.0                    blob_1.2.2                   
[137] BiocSingular_1.9.1            AnnotationHub_3.1.5          
[139] memoise_2.0.0                 irlba_2.3.3                  
[141] ape_5.5                      
```
</div>
