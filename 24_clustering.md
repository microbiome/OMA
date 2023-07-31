# Community Typing (Clustering) {#clustering}

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


```r
library(mia)
data("enterotype", package = "mia")
tse <- enterotype
```

Clustering is an unsupervised machine learning technique. The idea of
it is to find clusters in the data. A cluster is a group of
features/samples that share a pattern.  For example, with clustering, we
can find group of samples that share similar community
composition. There are multiple clustering algorithms available.

As mentioned before, clustering can be done either features or samples.
We will focus on the latter here. To learn about feature clustering,
check out chapter 6.3.

## Custom tools

*bluster* is a Bioconductor package providing tools for clustering data in 
in the `SummarizedExperiment` container. It offers multiple algorithms such 
as hierarchical clustering, DBSCAN, K-means, amongst others. The first thing 
to do when using this package is to load it, and transform the data if 
necessary, depending on your analysis goals.


```r
# Load dependencies
library(bluster)

# Apply transformation
tse <- transformAssay(tse, method = "relabundance")
```

The main focus here will be how to use mia's `cluster` function to 
cluster. It has multiple parameters that allow you to shape the result.
* The main new parameter allows you to choose the algorithm you want to use. 
In this example, we will use HclustParam which does the hierarchical 
clustering. This parameter itself has parameters on its own, you can check
them in the
[HclustParam documentation](https://rdrr.io/github/LTLA/bluster/man/HclustParam-class.html).
* Another parameter is `MARGIN`, which allows us to choose whether we want to
cluster the features or samples . Here we will cluster the latter.
* We will see the other parameters as we go along.


```r
# Simple use of the hierarchical clustering. Here, the default parameters
# set the cut height to half of the dendrogram height.
tse <- cluster(tse, assay.type = "relabundance", 
               MARGIN = "samples", HclustParam())

# Check the result contained in the clusters part of colData
colData(tse)$clusters
```

```
##    AM.AD.1    AM.AD.2  AM.F10.T1  AM.F10.T2    DA.AD.1   DA.AD.1T    DA.AD.2 
##          1          1          1          1          1          1          1 
##    DA.AD.3   DA.AD.3T    DA.AD.4    ES.AD.1    ES.AD.2    ES.AD.3    ES.AD.4 
##          1          1          1          2          3          2          1 
##    FR.AD.1    FR.AD.2    FR.AD.3    FR.AD.4    FR.AD.5    FR.AD.6    FR.AD.7 
##          1          1          2          1          1          1          1 
##    FR.AD.8    IT.AD.1    IT.AD.2    IT.AD.3    IT.AD.4    IT.AD.5    IT.AD.6 
##          2          2          1          2          3          1          1 
##    JP.AD.1    JP.AD.2    JP.AD.3    JP.AD.4    JP.AD.5    JP.AD.6    JP.AD.7 
##          2          1          1          2          1          2          2 
##    JP.AD.8    JP.AD.9    JP.IN.1    JP.IN.2    JP.IN.3    JP.IN.4     MH0001 
##          2          2          1          4          4          1          3 
##     MH0002     MH0003     MH0004     MH0005     MH0006     MH0007     MH0008 
##          1          1          1          1          1          1          2 
##     MH0009     MH0010     MH0011     MH0012     MH0013     MH0014     MH0015 
##          1          2          1          1          1          1          1 
##     MH0016     MH0017     MH0018     MH0019     MH0020     MH0021     MH0022 
##          2          1          3          2          2          1          1 
##     MH0023     MH0024     MH0025     MH0026     MH0027     MH0028     MH0030 
##          1          1          1          1          2          1          1 
##     MH0031     MH0032     MH0033     MH0034     MH0035     MH0036     MH0037 
##          1          3          1          3          1          1          1 
##     MH0038     MH0039     MH0040     MH0041     MH0042     MH0043     MH0044 
##          1          1          1          2          1          1          1 
##     MH0045     MH0046     MH0047     MH0048     MH0049     MH0050     MH0051 
##          2          1          1          1          1          1          1 
##     MH0052     MH0053     MH0054     MH0055     MH0056     MH0057     MH0058 
##          1          1          1          1          1          1          1 
##     MH0059     MH0060     MH0061     MH0062     MH0063     MH0064     MH0065 
##          1          1          2          1          1          1          1 
##     MH0066     MH0067     MH0068     MH0069     MH0070     MH0071     MH0072 
##          1          1          2          1          1          1          2 
##     MH0073     MH0074     MH0075     MH0076     MH0077     MH0078     MH0079 
##          2          1          1          1          1          1          1 
##     MH0080     MH0081     MH0082     MH0083     MH0084     MH0085     MH0086 
##          1          1          1          1          3          2          1 
##     TS1_V2    TS10_V2   TS100_V2 TS101.2_V2   TS103_V2   TS104_V2   TS105_V2 
##          5          6          5          7          8          5          7 
##   TS106_V2   TS107_V2   TS109_V2    TS11_V2   TS110_V2   TS111_V2   TS115_V2 
##          6          8          5          5          7          6          3 
##   TS116_V2   TS117_V2   TS118_V2   TS119_V2    TS12_V2   TS120_V2   TS124_V2 
##          8          6          5          6          5          3          6 
##   TS125_V2   TS126_V2   TS127_V2   TS128_V2   TS129_V2    TS13_V2   TS130_V2 
##          6          8          6          5          5          5          5 
##   TS131_V2   TS132_V2   TS133_V2   TS134_V2   TS135_V2   TS136_V2   TS137_V2 
##          6          5          5          6          5          7          5 
##   TS138_V2   TS139_V2    TS14_V2   TS140_V2   TS141_V2   TS142_V2   TS143_V2 
##          7          7          7          3          5          6          6 
##   TS144_V2   TS145_V2   TS146_V2   TS147_V2   TS148_V2   TS149_V2    TS15_V2 
##          3          5          8          6          5          5          3 
##   TS150_V2   TS151_V2   TS152_V2   TS153_V2 TS154.2_V2   TS155_V2   TS156_V2 
##          7          6          6          6          8          8          5 
##    TS16_V2   TS160_V2   TS161_V2   TS162_V2   TS163_V2   TS164_V2   TS165_V2 
##          5          3          6          5          5          5          8 
##   TS166_V2   TS167_V2   TS168_V2   TS169_V2    TS17_V2   TS170_V2   TS178_V2 
##          5          8          6          8          5          7          8 
##   TS179_V2   TS180_V2   TS181_V2   TS182_V2   TS183_V2   TS184_V2   TS185_V2 
##          7          5          5          5          5          5          6 
##   TS186_V2    TS19_V2   TS190_V2   TS191_V2   TS192_V2   TS193_V2   TS194_V2 
##          5          5          6          5          6          6          6 
##   TS195_V2     TS2_V2    TS20_V2    TS21_V2    TS22_V2    TS23_V2    TS25_V2 
##          5          8          3          6          7          7          6 
##    TS26_V2    TS27_V2    TS28_V2    TS29_V2     TS3_V2    TS30_V2    TS31_V2 
##          5          8          6          7          6          7          6 
##    TS32_V2    TS33_V2    TS34_V2    TS35_V2    TS37_V2    TS38_V2    TS39_V2 
##          3          7          6          3          5          7          5 
##     TS4_V2    TS43_V2    TS44_V2    TS49_V2     TS5_V2    TS50_V2    TS51_V2 
##          5          8          8          7          7          7          8 
##    TS55_V2    TS56_V2    TS57_V2     TS6_V2    TS61_V2    TS62_V2    TS63_V2 
##          6          5          9          3          6          5          6 
##    TS64_V2    TS65_V2    TS66_V2    TS67_V2    TS68_V2    TS69_V2     TS7_V2 
##          8          8          8          6          5          8          3 
##    TS70_V2    TS71_V2    TS72_V2    TS73_V2    TS74_V2    TS75_V2    TS76_V2 
##          8          6          6          5          8          6          5 
##    TS77_V2    TS78_V2     TS8_V2    TS82_V2    TS83_V2    TS84_V2    TS85_V2 
##          6          3          5          3          6          3          6 
##    TS86_V2    TS87_V2    TS88_V2    TS89_V2     TS9_V2    TS90_V2    TS91_V2 
##          5          8          5          6          5          5          7 
##    TS92_V2    TS94_V2    TS95_V2    TS96_V2    TS97_V2    TS98_V2  TS99.2_V2 
##          6          5          5          5          5          5          5 
## Levels: 1 2 3 4 5 6 7 8 9
```

Once the clustering on the samples is done, we can also plot the clusters.


```r
library(scater)

# Add the MDS dimensions for plotting
tse <- runMDS(tse, assay.type = "relabundance", 
              FUN = vegan::vegdist, method = "bray")

# Plot the clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")
```

![](24_clustering_files/figure-latex/bluster_sample_plot-1.pdf)<!-- --> 

We will now see different common algorithms and how to use them.

## Hierarchical clustering

The hierarchical clustering algorithm aims to find hierarchy between
samples/features. There are to approaches: agglomerative ("bottom-up")
and divisive ("top-down").

In agglomerative approach, each observation is first in a unique cluster.
The algorithm continues by agglomerating similar clusters. The divisive
approach starts with one cluster that contains all the observations. 
Clusters are split recursively to clusters that differ the most. 
The clustering ends when each cluster contains only one observation. 
In this algorithm, the similarity of two clusters is based on the distance
between them.

Hierarchical clustering can be visualized with a dendrogram tree. In each
splitting point, the tree is divided into two clusters leading to the
hierarchy.


```r
library(vegan)

# Load experimental data
tse <- enterotype
```

Hierarchical clustering requires 2 steps. 
1. Computation of the dissimilarities with a given distance.  
2. Clustering based on dissimilarities. 

Additionally, since sequencing data is compositional, we'll apply relative
transformation (as seen in the previous example).


```r
# Apply transformation
tse <- transformAssay(tse, method = "relabundance")

# Do the clustering
tse <- cluster(tse,
               assay.type = "relabundance",
               MARGIN = "samples",
               HclustParam(method = "complete",
                           dist.fun = vegdist,
                           metric = "bray"),
               full = TRUE,
               clust.col = "Hclust")
```

In this example, we wanted additional information on the clustering. To do so, 
we used the `full` parameter. We also computed the dissimilarities with the 
bray distance. Finally, the `clust.col` parameter allows us to choose the name
of the column in the colData (default name is `clusters`).

Next, we will plot the dendrogram, which is possible since we got the 
additional information from the clustering.


```r
library(dendextend)

# Get hclust data from metadata
hclust_data <- metadata(tse)$clusters$hclust

# Get the dendrogram object
dendro <- as.dendrogram(hclust_data)

# Plot dendrogram

dendro %>% set("labels", NULL) %>% plot()
```

![](24_clustering_files/figure-latex/hclust3-1.pdf)<!-- --> 

In our case, we cut the dendrogram in half by default. To know how many clusters
we have, we can check the colData.


```r
# Get the clusters
head(colData(tse)$Hclust)
```

```
##   AM.AD.1   AM.AD.2 AM.F10.T1 AM.F10.T2   DA.AD.1  DA.AD.1T 
##         1         1         2         1         3         3 
## 26 Levels: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ... 26
```

We can see that there are 26 clusters, but that probably isn't optimal since the 
the number of clusters was chosen arbitrarily. To determine the number of 
clusters, we can use the dendrogram. Usually the tree is split where the branch
length is the largest. However, as we can see from the dendrogram, clusters are
not clear. There are algorithms to identify the optimal number of clusters.

The `NbClust` library is useful to that end as it offers multiple methods to 
determine the optimal number of clusters. Here we will use the silhouette 
analysis to determine the optimal number of clusters. For each data point, this 
analysis measures the distance to other data points in the same cluster
(cohesion), and the distance to the other clusters (separation), establishing a
score. That score is then combined across the data points. `NbClust` does
this for multiple number of clusters and the best score corresponds to the
optimal number of clusters.


```r
library(NbClust)
diss <- metadata(tse)$clusters$dist

# Apply the silhouette analysis on the distance matrix
res <- NbClust(diss = diss, distance = NULL, method = "ward.D2",
               index = "silhouette")
```

```
## 
##  Only frey, mcclain, cindex, sihouette and dunn can be computed. To compute the other indices, data matrix is needed
```

```r
res$Best.nc
```

```
## Number_clusters     Value_Index 
##          2.0000          0.4783
```

Based on the result, let's divide observations into 2 clusters.


```r
library(dendextend)

# Get optimal number of clusters
k <- res$Best.nc[1]

# Making colors for 2 clusters
col_val_map <- randomcoloR::distinctColorPalette(k) %>%
    as.list() %>% 
    setNames(paste0("clust_", seq(k)))

dend <- color_branches(dendro, k = k, col = unlist(col_val_map))
labels(dend) <- NULL
plot(dend)
```

![](24_clustering_files/figure-latex/hclust6-1.pdf)<!-- --> 

## K-means clustering

Let's now try k-means clustering. Here observations are divided into clusters so
that the distances between observations and cluster centers are
minimized; an observation belongs to cluster whose center is the
nearest.

The algorithm starts by dividing observation to random clusters whose
number is defined by user. The centroids of the clusters are then
calculated. After that, observations' allocation to clusters are
updated so that the means are minimized. Again, the centroids are
calculated, and the algorithm continues iteratively until the assignments
do not change.

As an alternative to `NbClust` get the optimal number of clusters, we can 
visualize the silhouette analysis thanks to the library `factoextra`.


```r
library(factoextra)

# Convert dist object into matrix
diss <- as.matrix(diss)

# Perform silhouette analysis and plot the result
fviz_nbclust(diss, kmeans, method = "silhouette")
```

![](24_clustering_files/figure-latex/kmeans1-1.pdf)<!-- --> 

Based on the result of silhouette analysis, we confirm that 2 is the optimal
number of clusters in k-means clustering.


```r
# The first step is random, add seed for reproducibility
set.seed(15463)

# Perform k-means clustering with 2 clusters
km <- kmeans(diss, 2, nstart = 25)

# Add the result to colData
colData(tse)$clusters <- as.factor(km$cluster)

# Perform PCoA so that we can visualize clusters
tse <- runMDS(tse, assay.type = "relabundance", 
              FUN = vegan::vegdist, method = "bray")

# Plot PCoA and color clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")
```

![](24_clustering_files/figure-latex/kmeans2-1.pdf)<!-- --> 

## Dirichlet Multinomial Mixtures (DMM)

This section focus on DMM analysis. 

One technique that allows to search for groups of samples that are
similar to each other is the [Dirichlet-Multinomial Mixture
Model](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030126)
. In DMM, we first determine the number of clusters (k) that best fit the
data (model evidence) using Laplace approximation. After fitting the
model with k clusters, we obtain for each sample k probabilities that
reflect the probability that a sample belongs to the given cluster.

Let's cluster the data with DMM clustering. Since the dataset is large, the 
algorithm will take long computational time. Therefore, we use only a subset 
of the data; agglomerated by Phylum as a rank.


```r
# Get the data
data("GlobalPatterns", package = "mia")
tse <- GlobalPatterns

# Agglomerate by rank
tse <- mergeFeaturesByRank(tse, rank = "Phylum", agglomerateTree = TRUE)
```

Here we will further our use of `cluster` by renaming the clusters column in 
the metadata thanks to the `name` parameter.


```r
# Run the model and calculates the most likely number of clusters from 1 to 7
tse_dmm <- cluster(tse, name = "DMM", DmmParam(k = 1:7, type = "laplace"), 
                   MARGIN = "samples", full = TRUE)
```


```r
# The dmm info is stored in the metadata under the 'DMM' column
tse_dmm
```

```
## class: TreeSummarizedExperiment 
## dim: 67 26 
## metadata(2): agglomerated_by_rank DMM
## assays(1): counts
## rownames(67): Phylum:Crenarchaeota Phylum:Euryarchaeota ...
##   Phylum:Synergistetes Phylum:SR1
## rowData names(7): Kingdom Phylum ... Genus Species
## colnames(26): CL3 CC1 ... Even2 Even3
## colData names(8): X.SampleID Primer ... Description clusters
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (67 rows)
## rowTree: 1 phylo tree(s) (66 leaves)
## colLinks: NULL
## colTree: NULL
```

The following operation returns a list of DMM objects for closer investigation.


```r
metadata(tse_dmm)$DMM$dmm
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
## Laplace: 7792 BIC: 8357 AIC: 8187 
## 
## [[5]]
## class: DMN 
## k: 5 
## samples x taxa: 26 x 67 
## Laplace: 7844 BIC: 8548 AIC: 8335 
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
## Laplace: 8076 BIC: 9100 AIC: 8801
```

We can see the Laplace approximation (model evidence)
for each model of the k models.


```r
BiocManager::install("microbiome/miaViz")
library(miaViz)
plotDMNFit(tse_dmm, type = "laplace", name = "DMM")
```

![](24_clustering_files/figure-latex/dmm5-1.pdf)<!-- --> 

On the graph, we can see that the best number of clusters is 2. We can confirm
that with the following operation.


```r
# Get the model that has the best fit
bestFit <- metadata(tse_dmm)$DMM$dmm[[metadata(tse_dmm)$DMM$best]]
bestFit
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
dmm_group <- calculateDMNgroup(tse_dmm, variable = "SampleType", 
                               assay.type = "counts", k = 2, 
                               seed = .Machine$integer.max)

dmm_group
```

```
## class: DMNGroup 
## summary:
##                    k samples taxa    NLE  LogDet Laplace    BIC  AIC
## Feces              2       4   67 1078.3 -106.26   901.1 1171.9 1213
## Freshwater         2       2   67  889.6  -97.20   716.9  936.4 1025
## Freshwater (creek) 2       3   67 1600.3  862.19  1907.3 1674.5 1735
## Mock               2       3   67 1008.4  -55.40   856.6 1082.5 1143
## Ocean              2       3   67 1096.7  -56.66   944.3 1170.9 1232
## Sediment (estuary) 2       3   67 1195.5   18.63  1080.8 1269.7 1331
## Skin               2       3   67  992.6  -85.05   826.1 1066.8 1128
## Soil               2       3   67 1380.3   11.20  1261.8 1454.5 1515
## Tongue             2       2   67  783.0 -107.79   605.0  829.8  918
```

Mixture weights (rough measure of the cluster size).


```r
DirichletMultinomial::mixturewt(bestFit)
```

```
##       pi theta
## 1 0.5385 20.58
## 2 0.4615 15.32
```

It's also possible to get the samples-cluster assignment probabilities: how
probable it is that each sample belongs to each cluster


```r
prob <- metadata(tse_dmm)$DMM$prob
head(prob)
```

```
##                 1         2
## CL3     1.000e+00 4.501e-17
## CC1     1.000e+00 3.417e-22
## SV1     1.000e+00 1.712e-12
## M31Fcsw 7.425e-26 1.000e+00
## M11Fcsw 1.093e-16 1.000e+00
## M31Plmr 1.152e-13 1.000e+00
```

We can also know the contribution of each taxa to each component


```r
head(DirichletMultinomial::fitted(bestFit))
```

```
##                          [,1]      [,2]
## Phylum:Crenarchaeota  0.30381 0.1354045
## Phylum:Euryarchaeota  0.23114 0.1468879
## Phylum:Actinobacteria 1.21375 1.0581523
## Phylum:Spirochaetes   0.21393 0.1318064
## Phylum:MVP-15         0.02982 0.0007669
## Phylum:Proteobacteria 6.84470 1.8114033
```

Finally, to be able to visualize our data and clusters, we start by 
computing the euclidean PCoA and storing it as a data frame.


```r
# Do clr transformation. Pseudocount is added, because data contains zeros.
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformAssay(tse, assay.type = "pseudo", method = "relabundance")
tse <- transformAssay(tse, "relabundance", method = "clr")

# Do principal coordinate analysis
df <- calculateMDS(tse, assay.type = "clr", method = "euclidean")

# Create a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = df[, 1], pcoa2 = df[, 2])
```


```r
# Create a data frame that contains principal coordinates and DMM information
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = colData(tse_dmm)$clusters)

# Create a plot
euclidean_dmm_plot <- ggplot(data = euclidean_dmm_pcoa_df,
                             aes(x = pcoa1, y = pcoa2, color = dmm_component)) +
    geom_point() +
    labs(x = "Coordinate 1",y = "Coordinate 2", 
         title = "PCoA with Aitchison distances") +
    theme(plot.title = element_text(size = 12, # makes titles smaller
                                    hjust = 0.5)) 

euclidean_dmm_plot
```

![](24_clustering_files/figure-latex/dmm11-1.pdf)<!-- --> 

## Graph-based clustering

Another approach for discovering communities within the samples of the
data, is to run community detection algorithms after building a
graph. The following demonstration builds a graph based on the k
nearest-neighbors and performs the community detection on the fly.

Here, we will be using the `cluster` function with a graph-based clustering 
function, enabling the community detection task.

The algorithm used is "short random walks" [@Pons2006]. The graph is
constructed using different k values (the number of nearest neighbors
to consider during graph construction) using the robust centered log
ratio (rclr) assay data. Then plotting the communities using UMAP
[@McInnes2018] ordination as a visual exploration aid. Let us cluster the 
`enterotype` dataset.


```r
library(patchwork) # For arranging several plots as a grid

# Get enterotype dataset and transform data with rclr
tse <- enterotype
tse <- transformAssay(tse, method = "rclr")

# Perform and store UMAP
tse <- runUMAP(tse, name = "UMAP", assay.type = "rclr")

# Set different k values to test
k <- c(2, 3, 5, 10)

ClustAndPlot <- function(x) {
    # Add the clustering data from the short random walks algorithm to the TSE
    tse <- cluster(tse, assay.type = "rclr", 
                   MARGIN = "col", NNGraphParam(k = x))
    
    # Plot the results of the clustering as a color for each sample
    plotUMAP(tse, colour_by = I(colData(tse)$clusters)) +
    labs(title = paste0("k = ", x))
}

# Apply the function for different k values
plots <- lapply(k, ClustAndPlot)

# Display plots in a grid
plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plot_layout(ncol = 4)
```

![](24_clustering_files/figure-latex/NNGraph1-1.pdf)<!-- --> 

In this graph, we can clearly see the impact of the k choice on the quality
of the clustering

Similarly, the _`bluster`_ [@R_bluster] package offers clustering
diagnostics that can be used for judging the clustering quality (see
[Assorted clustering
diagnostics](http://bioconductor.org/packages/release/bioc/vignettes/bluster/inst/doc/diagnostics.html)).
In the following, Silhouette width as a diagnostic tool is computed
and results are visualized for each case presented earlier. For more
about Silhouettes read [@Rousseeuw1987].


```r
ClustDiagPlot <- function(x) {
  # Get the clustering results
    tse <- cluster(tse, assay.type = "rclr", 
                   MARGIN = "col", NNGraphParam(k = x))

    # Compute the diagnostic info
    sil <- approxSilhouette(t(assays(tse)$rclr), colData(tse)$clusters)

    # Plot as a boxplot to observe cluster separation
    boxplot(split(sil$width, colData(tse)$clusters), main = paste0("k = ", x))
}
# Apply the function for different k values
res <- lapply(k, ClustDiagPlot)
```

![](24_clustering_files/figure-latex/NNGraph2-1.pdf)<!-- --> ![](24_clustering_files/figure-latex/NNGraph2-2.pdf)<!-- --> ![](24_clustering_files/figure-latex/NNGraph2-3.pdf)<!-- --> ![](24_clustering_files/figure-latex/NNGraph2-4.pdf)<!-- --> 

## Biclustering

Biclustering methods cluster rows and columns simultaneously in order
to find subsets of correlated features/samples.

Here, we use following packages:

-   [biclust](https://cran.r-project.org/web/packages/biclust/index.html)
-   [cobiclust](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13582)

_cobiclust_ is especially developed for microbiome data whereas _biclust_ is more
general method. In this section, we show two different cases and example 
solutions to apply biclustering to them. 

1.   Taxa vs samples
2.   Taxa vs biomolecule/biomarker

Biclusters can be visualized using heatmap or boxplot, for
instance. For checking purposes, also scatter plot might be valid
choice.

Check more ideas for heatmaps from chapters \@ref(viz-chapter) and
\@ref(microbiome-community).

### Taxa vs samples

When you have microbial abundance matrices, we suggest to use
_cobiclust_ which is designed for microbial data.

Load example data

```r
library(cobiclust)
data("HintikkaXOData")
mae <- HintikkaXOData
```

Only the most prevalent taxa are included in analysis. 


```r
# Subset data in the first experiment
mae[[1]] <- subsetByPrevalentFeatures(mae[[1]], rank = "Genus", 
                                      prevalence = 0.2, 
                                      detection = 0.001)

# rclr-transform in the first experiment
mae[[1]] <- transformAssay(mae[[1]], method = "rclr")
```

_cobiclust_ takes counts table as an input and gives _cobiclust_ object as an output.
It includes clusters for taxa and samples. 


```r
# Do clustering using counts table
clusters <- cobiclust(assay(mae[[1]], "counts"))

# Get clusters
row_clusters <- clusters$classification$rowclass
col_clusters <- clusters$classification$colclass

# Add clusters to rowdata and coldata
rowData(mae[[1]])$clusters <- factor(row_clusters)
colData(mae[[1]])$clusters <- factor(col_clusters)

# Order data based on clusters
mae[[1]] <- mae[[1]][order(rowData(mae[[1]])$clusters),
                     order(colData(mae[[1]])$clusters)]

# Print clusters
clusters$classification
```

```
## $rowclass
##  [1] 1 1 1 1 2 2 1 1 1 1 1 1 2 2 2 2 1 2 1 1 2 1 2 2 1 1 2 1 1 1 1 1 2 1 1 2 1 1
## [39] 1 1 1 1 1 1 1 1 1 2 1 2 1 1 1 2 1 1 1
## 
## $colclass
##  C1  C2  C3  C4  C5  C6  C7  C8  C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 
##   1   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2 
## C21 C22 C23 C24 C25 C26 C27 C28 C29 C30 C31 C32 C33 C34 C35 C36 C37 C38 C39 C40 
##   2   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   1
```

Next we can plot clusters. Annotated heatmap is a common choice.


```r
library(pheatmap)
# z-transform for heatmap
mae[[1]] <- transformAssay(mae[[1]], assay.type = "rclr",
                            MARGIN = "features", method = "z", name = "rclr_z")

# Create annotations. When column names are equal, they should share levels.
# Here samples include 3 clusters, and taxa 2. That is why we have to make
# column names unique.
annotation_col <- data.frame(colData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_col) <- "col_clusters"

annotation_row <- data.frame(rowData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_row) <- "row_clusters"
```

Plot the heatmap.


```r
pheatmap(assay(mae[[1]], "rclr_z"), cluster_rows = F, cluster_cols = F,
         annotation_col = annotation_col, annotation_row = annotation_row)
```

![](24_clustering_files/figure-latex/cobiclust_3b-1.pdf)<!-- --> 

Boxplot is commonly used to summarize the results:


```r
library(ggplot2)
library(patchwork)

# ggplot requires data in melted format
melt_assay <- meltAssay(mae[[1]], assay.type = "rclr", 
                        add_col_data = T, add_row_data = T)

# patchwork two plots side-by-side
p1 <- ggplot(melt_assay) +
    geom_boxplot(aes(x = clusters.x, y = rclr)) +
    labs(x = "Taxa clusters")

p2 <- ggplot(melt_assay) +
    geom_boxplot(aes(x = clusters.y, y = rclr)) +
    labs(x = "Sample clusters")

p1 + p2
```

![](24_clustering_files/figure-latex/cobiclust_4-1.pdf)<!-- --> 

### Taxa vs biomolecules

Here, we analyze cross-correlation between taxa and metabolites. This
is a case, where we use _biclust_ method which is suitable for numeric
matrices in general. First we pre-process the data.


```r
# Samples must be in equal order
# (Only 1st experiment was ordered in cobiclust step leading to unequal order)
mae[[1]] <- mae[[1]][, colnames(mae[[2]])]

# Make rownames unique since it is required by other steps
rownames(mae[[1]]) <- make.unique(rownames(mae[[1]]))

# Transform the metabolites to be in log basis
mae[[2]] <- transformAssay(mae[[2]], assay.type = "nmr", method = "log10")

# Add missing data to the metabolites
replace_na <- function(row) {
    na_indices <- which(is.na(row))
    non_na_values <- row[!is.na(row)]
    row[na_indices] <- sample(non_na_values, length(na_indices), replace = TRUE)
    row
}
assay(mae[[2]], "log10") <- t(apply(assay(mae[[2]], "log10"), 1, replace_na))
```

Next, we compute the spearman correlation matrix.


```r
# Calculate correlations
corr <- getExperimentCrossCorrelation(mae, 1, 2, assay.type1 = "rclr",
                                      assay.type2 = "log10", mode = "matrix",
                                      correlation = "spearman")
```

_biclust_ takes a matrix as an input and returns a _biclust_ object. 


```r
library(biclust)
# Set seed for reproducibility
set.seed(3973)

# Find biclusters
bc <- biclust(corr, method = BCPlaid(), verbose = FALSE)

bc
```

```
## 
## An object of class Biclust 
## 
## call:
## 	biclust(x = corr, method = BCPlaid(), verbose = FALSE)
## 
## Number of Clusters found:  6 
## 
## First  5  Cluster sizes:
##                    BC 1 BC 2 BC 3 BC 4 BC 5
## Number of Rows:      11    9    6    3    2
## Number of Columns:   15   13    9    9   15
```

The object includes cluster information. However compared to
_cobiclust_, _biclust_ object includes only information about clusters
that were found, not general cluster.

Meaning that if one cluster size of 5 features was found out of 20 features, 
those 15 features do not belong to any cluster. That is why we have to create an
additional cluster for features/samples that are not assigned into any cluster.


```r
# Functions for obtaining biclust information

# Get clusters for rows and columns
.get_biclusters_from_biclust <- function(bc, assay) {
    # Get cluster information for columns and rows
    bc_columns <- t(bc@NumberxCol)
    bc_columns <- data.frame(bc_columns)
    bc_rows <- bc@RowxNumber
    bc_rows <- data.frame(bc_rows)

    # Get data into right format
    bc_columns <- .manipulate_bc_data(bc_columns, assay, "col")
    bc_rows <- .manipulate_bc_data(bc_rows, assay, "row")
    
    return(list(bc_columns = bc_columns, bc_rows = bc_rows))
}

# Input clusters, and how many observations there should be, i.e.,
# the number of samples or features
.manipulate_bc_data <- function(bc_clusters, assay, row_col) {
    # Get right dimension
    dim <- ifelse(row_col == "col", ncol(assay), nrow(assay))
    # Get column/row names
    if (row_col == "col") {
        names <- colnames(assay)
    } else {
        names <- rownames(assay)
    }

    # If no clusters were found, create one. Otherwise create additional
    # cluster which
    # contain those samples that are not included in clusters that were found.
    if (nrow(bc_clusters) != dim) {
        bc_clusters <- data.frame(cluster = rep(TRUE, dim))
    } else {
        # Create additional cluster that includes those samples/features that
        # are not included in other clusters.
        vec <- ifelse(rowSums(bc_clusters) > 0, FALSE, TRUE)

        # If additional cluster contains samples, then add it
        if (any(vec)) {
            bc_clusters <- cbind(bc_clusters, vec)
        }
    }
    
    # Adjust row and column names
    rownames(bc_clusters) <- names
    colnames(bc_clusters) <- paste0("cluster_", 1:ncol(bc_clusters))
    return(bc_clusters)
}
```


```r
# Get biclusters
bcs <- .get_biclusters_from_biclust(bc, corr)

bicluster_rows <- bcs$bc_rows
bicluster_columns <- bcs$bc_columns

# Print biclusters for rows
head(bicluster_rows)
```

```
##                           cluster_1 cluster_2 cluster_3 cluster_4 cluster_5
## D_5__Escherichia-Shigella     FALSE     FALSE     FALSE     FALSE     FALSE
## D_5__Ruminiclostridium 5       TRUE     FALSE      TRUE     FALSE     FALSE
## D_5__Lactobacillus            FALSE     FALSE     FALSE     FALSE     FALSE
## D_5__uncultured               FALSE     FALSE     FALSE     FALSE     FALSE
## D_5__uncultured bacterium     FALSE     FALSE     FALSE     FALSE     FALSE
## D_5__Lactococcus              FALSE     FALSE     FALSE     FALSE     FALSE
##                           cluster_6 cluster_7
## D_5__Escherichia-Shigella     FALSE      TRUE
## D_5__Ruminiclostridium 5      FALSE     FALSE
## D_5__Lactobacillus            FALSE      TRUE
## D_5__uncultured                TRUE     FALSE
## D_5__uncultured bacterium     FALSE      TRUE
## D_5__Lactococcus              FALSE      TRUE
```

Let's collect information for the scatter plot. 


```r
# Function for obtaining sample-wise sum, mean, median, and mean variance
# for each cluster

.sum_mean_median_var <- function(tse1, tse2, assay.type1, assay.type2, clusters1, clusters2) {
    list <- list()
    # Create a data frame that includes all the information
    for (i in 1:ncol(clusters1)) {
        # Subset data based on cluster
        tse_subset1 <- tse1[clusters1[, i], ]
        tse_subset2 <- tse2[clusters2[, i], ]
        # Get assay
        assay1 <- assay(tse_subset1, assay.type1)
        assay2 <- assay(tse_subset2, assay.type2)
        # Calculate sum, mean, median, and mean variance
        sum1 <- colSums2(assay1, na.rm = T)
        mean1 <- colMeans2(assay1, na.rm = T)
        median1 <- colMedians(assay1, na.rm = T)
        var1 <- colVars(assay1, na.rm = T)
        
        sum2 <- colSums2(assay2, na.rm = T)
        mean2 <- colMeans2(assay2, na.rm = T)
        median2 <- colMedians(assay2, na.rm = T)
        var2 <- colVars(assay2, na.rm = T)
        
        list[[i]] <- data.frame(sample = colnames(tse1), sum1, sum2, mean1, 
                                 mean2, median1, median2, var1, var2)
    }
    return(list)
}

# Calculate info
df <- .sum_mean_median_var(mae[[1]], mae[[2]], "rclr", "log10", bicluster_rows, bicluster_columns)
```

Now we can create a scatter plot. X-axis includes median clr abundance
of microbiome and y-axis median absolute concentration of each
metabolite. Each data point represents a single sample.

From the plots, we can see that there is low negative correlation in
both cluster 1 and 3.  This means that when abundance of bacteria
belonging to cluster 1 or 3 is higher, the concentration of
metabolites of cluster 1 or 3 is lower, and vice versa.


```r
pics <- list()
for (i in seq_along(df)) {
    pics[[i]] <- ggplot(df[[i]]) +
        geom_point(aes(x = median1, y = median2)) +
        labs(title = paste0("Cluster ", i), x = "Taxa (rclr median)",
             y = "Metabolites (abs. median)")
    print(pics[[i]])
}
```


\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-1} 
\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-2} 
\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-3} 
\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-4} 
\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-5} 
\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-6} 
\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-7} 

```r
pics[[1]] + pics[[2]] + pics[[3]]
```


\includegraphics[width=0.33\linewidth]{24_clustering_files/figure-latex/biclust_7-8} 

_pheatmap_ does not allow boolean values, so they must be converted into factors.


```r
bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))
```

Again, we can plot clusters with heatmap.


```r
# Adjust colors for all clusters
if (ncol(bicluster_rows) > ncol(bicluster_columns)) {
    cluster_names <- colnames(bicluster_rows)
} else {
    cluster_names <- colnames(bicluster_columns)
}
annotation_colors <- list()
for (name in cluster_names) {
    annotation_colors[[name]] <- c("TRUE" = "red", "FALSE" = "white")
}

# Create a heatmap
pheatmap(corr, cluster_cols = F, cluster_rows = F,
         annotation_col = bicluster_columns, annotation_row = bicluster_rows,
         annotation_colors = annotation_colors)
```

![](24_clustering_files/figure-latex/biclust_9-1.pdf)<!-- --> 

## Additional Community Typing

For more community typing techniques applied to the 'SprockettTHData'
data set, see the attached .Rmd file.

Link:

   * [Rmd](add-comm-typing.Rmd)
