# Visualization {#viz-chapter}

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

Whether a data set contains information on a microbial community or it originates from a different source, the way that data are visualized inevitably shapes how they will be interpreted, and motivates the next steps of the analysis.

A large variety of graphing methods belong to microbial analysis, but only few are the choices that will return useful answers about the data. Therefore, knowledge on the available tools and their possible applications plays an important role in selecting the most suitable method for the asked question.

This chapter introduces the reader to a number of visualization techniques found in this book, such as:

* bar plots
* box plots
* heatmaps
* ordination charts
* regression charts
* trees

The toolkit which provides the essential plotting functionality includes the following packages:

* _patchwork_, _cowplot_, _ggpubr_ and _gridExtra_: plot layout and multi-panel plotting
* _miaViz_: specific visualization tools for `TreeSummaizedExperiment` objects
* _scater_: specific visualization tools for `SingleCellExperiment` objects
* _ggplot2_, _pheatmap_, _ggtree_, _sechm_: composition heatmaps
* _ANCOMBC_, _ALDEx2_ and _Maaslin2_: visual differential abundance
* _fido_: tree-based methods for differential abundance
* _plotly_: animated and 3D plotting

For systematic and extensive tutorials on the visual tools available in _mia_, readers can refer to the following material:

* [microbiome tutorials](https://microbiome.github.io/tutorials/)



## Pre-analysis exploration

### Accessing row and column data

`SCE` and `TSE` objects contain multiple layers of information in the form of
rows, columns and meta data. The _scater_ package supports in accessing,
modifying and graphing the meta data related to features as well as samples.


```r
# list row meta data
names(rowData(tse))
```

```
## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
```

```r
# list column meta data
names(colData(tse))
```

```
## [1] "X.SampleID"               "Primer"                  
## [3] "Final_Barcode"            "Barcode_truncated_plus_T"
## [5] "Barcode_full_length"      "SampleType"              
## [7] "Description"
```

Such meta data can be directly plotted with the functions `plotRowData` and `plotColData`.


```r
# obtain QC data
tse <- addPerCellQC(tse)
tse <- addPerFeatureQC(tse)
# plot QC Mean against Species
plotRowData(tse, "mean", "Species") +
  theme(axis.text.x = element_blank()) +
  labs(x = "Species", y = "QC Mean")
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-3-1.pdf)<!-- --> 

```r
# plot QC Sum against Sample ID, colour-labeled by Sample Type
plotColData(tse, "sum", "X.SampleID", colour_by = "SampleType") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample ID", y = "QC Sum")
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-3-2.pdf)<!-- --> 

Alternatively, they can be converted to a `data.frame` object and passed to `ggplot`.


```r
# store colData into a data frame
coldata <- as.data.frame(colData(tse))
# plot Number of Samples against Sampling Site
ggplot(coldata, aes(x = SampleType)) +
  geom_bar(width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sampling Site",
       y = "Number of Samples")
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-4-1.pdf)<!-- --> 

Further methods of application can be found in the chapters \@ref(qc) and
\@ref(richness) and in a few [external tutorials](https://github.com/davismcc/scater_tutorials_open_data)
with open data. Additionally, `rowData` and `colData` allow manipulation and
subsetting of large data sets into smaller units, as explained in chapter \@ref(datamanipulation).

### Viewing abundance and prevalence patterns

Prior-to-analysis exploration may involve questions such as how microorganisms
are distributed across samples (abundance) and what microorganisms are present
in most of the samples (prevalence). The information on abundance and prevalence
can be summarized into a **jitter** or **density plot** and a **tree**,
respectively, with the _miaViz_ package.

Specifically, the functions `plotAbundance`, `plotAbundanceDensity` and
`plotRowTree` are used, and examples on their usage are discussed throughout
chapter \@ref(quality-control).

## Diversity estimation

Alpha diversity is commonly measured as one of the diversity indices explained
in chapter \@ref(community-diversity). Because the focus lies on each sample
separately, one-dimensional plots, such as **scatter**, **violin** and
**box plots**, are suitable. 

Beta diversity is generally evaluated as one of the dissimilarity indices
reported in chapter \@ref(community-similarity). Unlike alpha diversity,
samples are compared collectively to estimate the heterogeneity across them,
therefore multidimensional plots, such as **Shepard** and **ordination plots**
are suitable.

|                         | alpha diversity            | beta diversity            |
|:-----------------------:|:--------------------------:|:-------------------------:|
| used metrics            | diversity indices          | dissimilarity indices     |
|                         |                            |                           | 
| metric dimensionality   | one-dimensional            | multidimensional          |
|                         |                            |                           | 
| suitable visualization  | scatter, violin, box plots | Shepard, ordination plots |

As a conclusion, visualization techniques for alpha and beta diversity significantly differ from one another.

### Alpha diversity with scatter, violin and box plots

The basic method to visualize the diversity values assigned to the different
samples in a `TSE` object includes the following, where each data point represents one sample:


```r
# estimate shannon diversity index
tse <- mia::estimateDiversity(tse, 
                              abund_values = "counts",
                              index = "shannon", 
                              name = "shannon")
# plot shannon diversity index, colour-labeled by Sample Type
plotColData(tse, "shannon", colour_by = "SampleType")
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-5-1.pdf)<!-- --> 

The several indices available for the evaluation of alpha diversity often return
slightly divergent results, which can be visually compared with a multiple violin
or box plot. For this purpose, `plotColData` (for violin plots) or `ggplot`
(for box plots) are recursively applied to a number of diversity indices with
the function `lapply` and the multi-panel plotting functionality of the
_patchwork_ package is then exploited.


```r
# estimate faith diversity index
tse <- mia::estimateFaith(tse,
                          abund_values = "counts")
# store colData into a data frame
coldata <- as.data.frame(colData(tse))
# generate plots for shannon and faith indices
# and store them into a list
plots <- lapply(c("shannon", "faith"),
                function(i) ggplot(coldata, aes_string(y = i)) +
                  geom_boxplot() +
                  theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank()))
# combine plots with patchwork
plots[[1]] + plots[[2]]
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 

The analogous output in the form of a violin plot is obtained in chapter
\@ref(faith-diversity). In addition, box plots that group samples according to
certain information, such as origin, sex, age and health condition, can be
labeled with p-values for significant differences with the package _ggsignif_
package, as shown in chapter \@ref(estimate-diversity).

### Beta diversity with Shepard and coordination plots

The _scater_ package offers the general function `plotReducedDim`. In its basic
form, it takes a `TSE` object and the results on sample similarity stored in the
same object, which can be evaluated with the following coordination methods:

* `runMDS`
* `runNMDS`
* `runPCA`
* `runTSNE`
* `runUMAP`

Since these clustering techniques allow for multiple coordinates or components,
**coordination plots** can also span multiple dimensions, which is explained in chapter \@ref{extras}.


```r
# perform NMDS coordination method
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               name = "NMDS")
```

```
## initial  value 47.733208 
## iter   5 value 33.853364
## iter  10 value 32.891200
## final  value 32.823570 
## converged
```

```r
# plot results of a 2-component NMDS on tse,
# coloured-scaled by shannon diversity index
plotReducedDim(tse, "NMDS", colour_by = "shannon")
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-7-1.pdf)<!-- --> 

Multiple combinations of coordinates or dimensions can also be integrated into a multi-panel arrangement.


```r
# perform MDS coordination method
tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "bray",
              name = "MDS",
              exprs_values = "counts",
              ncomponents = 3)
# plot results of a 3-component MDS on tse,
# coloured-scaled by faith diversity index
plotReducedDim(tse, "MDS", ncomponents = c(1:3), colour_by = "faith")
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-8-1.pdf)<!-- --> 

Similarly to iterating `plotColData` over indices of alpha diversity, `lapply`
can be used in combination with _patchwork_ to recursively apply `plotReducedDim`
and visually compare results among various coordination methods.


```r
# generate plots for MDS and NMDS methods
# and store them into a list
plots <- lapply(c("MDS", "NMDS"),
                plotReducedDim,
                object = tse,
                colour_by = "shannon")
# combine plots with patchwork
plots[[1]] + plots[[2]] +
  plot_layout(guides = "collect")
```

![](19_visualization_techniques_files/figure-latex/unnamed-chunk-9-1.pdf)<!-- --> 

For similar examples, readers are referred to chapter \@ref(community-similarity).
Further material on the graphic capabilities of _patchwork_ is available in its
[official package tutorial](https://patchwork.data-imaginist.com/articles/patchwork.html).

## Statistical analysis

### Heatmaps

As described in chapter \@ref(visual-composition), bar plots and heatmaps can
offer a useful insight into the composition of a community. Simple methods involve
the functions `plotAbundance` and `geom_tile` in combination with `scale_fill_gradientn`
from the packages _miaViz_ and _ggplot2_, respectively.

For instance, below the composition of multiple samples (x axis) is reported
in terms of relative abundances (y axis) for the top 10 taxa at the Order rank.
Bar plots and heatmaps with analogous information at the Phylum level are
available in the aforementioned chapter.


```r
# agglomerate tse by Order
tse_order <- agglomerateByRank(tse,
                                rank = "Order",
                                onRankOnly = TRUE)
# transform counts into relative abundance
tse_order <- transformCounts(tse_order,
                              abund_values = "counts",
                              method = "relabundance")
# get top orders
top_taxa <- getTopTaxa(tse_order,
                       top = 10,
                       abund_values = "relabundance")
# leave only names for top 10 orders and label the rest with "Other"
order_renamed <- lapply(rowData(tse_order)$Order,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_order)$Order <- as.character(order_renamed)
# plot composition as a bar plot
plotAbundance(tse_order,
              abund_values = "relabundance",
              rank = "Order",
              order_rank_by = "abund",
              order_sample_by = "Clostridiales")
```

![](19_visualization_techniques_files/figure-latex/plotAbundance1-1.pdf)<!-- --> 

To add a sample annotation, you can combine plots that you get from the output
pf _plotAbundance_.


```r
# Create plots
plots <- plotAbundance(tse_order, abund_values = "relabundance", rank = "Order",
                       order_rank_by = "abund", order_sample_by = "Clostridiales",
                       features = "SampleType")

# Modify the legend of the first plot to be smaller 
plots[[1]] <- plots[[1]] +
    theme(legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))

# Modify the legend of the second plot to be smaller 
plots[[2]] <- plots[[2]] +
    theme(legend.key.height = unit(0.3, 'cm'),
          legend.key.width = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.direction = "vertical")

# Load required packages
if( !require("ggpubr") ){
    install.packages("ggpubr")
    library("ggpubr")
}
# Load required packages
if( !require("patchwork") ){
    install.packages("patchwork")
    library("patchwork")
}

# Combine legends
legend <- wrap_plots(as_ggplot(get_legend(plots[[1]])), as_ggplot(get_legend(plots[[2]])), ncol = 1) 

# Remove legends from the plots
plots[[1]] <- plots[[1]] + theme(legend.position = "none")
plots[[2]] <- plots[[2]] + theme(legend.position = "none", axis.title.x=element_blank()) 

# Combine plots
plot <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend
wrap_plots(plot, legend, nrow = 1, widths = c(2, 1))
```

![](19_visualization_techniques_files/figure-latex/plotAbundance2-1.pdf)<!-- --> 

For more sophisticated visualizations than those produced with `plotAbundance`
and _ggplot2_, the packages _pheatmap_ and _sechm_ provide methods to include
feature and sample clusters in a heatmap, along with further functionality.


```r
# Agglomerate tse by phylum
tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)
# Add clr-transformation on samples
tse_phylum <- transformSamples(tse_phylum, method = "relabundance", pseudocount = 1)
tse_phylum <- transformSamples(tse_phylum, method = "clr", abund_values = "relabundance")
# Add z-transformation on features (taxa)
tse_phylum <- transformFeatures(tse_phylum, abund_values = "clr", 
                                method = "z", name = "clr_z")

# Takes subset: only samples from feces, skin, or tongue
tse_phylum_subset <- tse_phylum[ , colData(tse_phylum)$SampleType %in% c("Feces", "Skin", "Tongue") ]

# Does clr-transformation
tse_phylum_subset <- transformSamples(tse_phylum_subset, method = "relabundance", pseudocount = 1)
tse_phylum_subset <- transformSamples(tse_phylum_subset, method = "clr", abund_values = "relabundance")
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

![](19_visualization_techniques_files/figure-latex/pheatmap1-1.pdf)<!-- --> 

We can cluster both samples and features hierarchically and add them to the
x and y axes of the heatmap, respectively.


```r
# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

# to view the tree, run
# taxa_tree
```

Based on phylo tree, we decide to create three clusters. 


```r
# Creates clusters
taxa_clusters <- cutree(tree = taxa_hclust, k = 3)

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters
```

```
##                  clusters
## Chloroflexi             3
## Actinobacteria          3
## Crenarchaeota           3
## Planctomycetes          3
## Gemmatimonadetes        3
## Thermi                  3
## Acidobacteria           3
## Spirochaetes            2
## Fusobacteria            2
## SR1                     2
## Cyanobacteria           2
## Proteobacteria          2
## Synergistetes           2
## Lentisphaerae           1
## Bacteroidetes           1
## Verrucomicrobia         1
## Tenericutes             1
## Firmicutes              1
## Euryarchaeota           1
## SAR406                  1
```

The information on the clusters is then added to the feature meta data.


```r
# Adds information to rowData
rowData(tse_phylum_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), 
                                                                 rownames(tse_phylum_subset))), ]

# Prints taxa and their clusters
rowData(tse_phylum_subset)$clusters
```

```
##  [1] 1 1 2 3 2 2 1 1 1 1 3 2 3 3 3 2 2 3 3 1
## Levels: 1 2 3
```

Similarly, samples are hierarchically grouped into clusters, the most suitable
number of clusters for the plot is selected and the new information is stored
into the sample meta data.


```r
# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
# sample_tree

# Creates clusters
sample_clusters <- factor(cutree(tree = sample_hclust, k = 3))

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 

# Order data based on 
tse_phylum_subset <- tse_phylum_subset[ , rownames(sample_data)]

# Add sample type data
sample_data$sample_types <- unfactor(colData(tse_phylum_subset)$SampleType)

sample_data
```

```
##         clusters sample_types
## M11Plmr        2         Skin
## M31Plmr        2         Skin
## F21Plmr        2         Skin
## M31Fcsw        1        Feces
## M11Fcsw        1        Feces
## TS28           3        Feces
## TS29           3        Feces
## M31Tong        3       Tongue
## M11Tong        3       Tongue
```

Now we can create heatmap with additional annotations.


```r
# Determines the scaling of colorss
# Scale colors
breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)

pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         breaks = breaks,
         color = colors)
```

![](19_visualization_techniques_files/figure-latex/pheatmap6-1.pdf)<!-- --> 

The package _sechm_ allows for further visual capabilities and flexibility.
In this case, the clustering step is automatically performed by the plotting
function and does not need to be executed in advance.


```r
# Stores annotation colros to metadata
metadata(tse_phylum_subset)$anno_colors$SampleType <- c(Feces = "blue", 
                                                        Skin = "red", 
                                                        Tongue = "gray")

# Create a plot
sechm(tse_phylum_subset, 
      features = rownames(tse_phylum_subset), 
      assayName = "clr", 
      do.scale = TRUE, 
      top_annotation = c("SampleType"), 
      gaps_at = "SampleType",
      cluster_cols = TRUE, cluster_rows = TRUE)
```

![](19_visualization_techniques_files/figure-latex/sechm-1.pdf)<!-- --> 

It is also possible to create an analogous heatmap by just using the _ggplot2_
package. However, a relatively long code is required to generate an identical output.


```r
# Add feature names to column as a factor
taxa_clusters$Feature <- rownames(taxa_clusters)
taxa_clusters$Feature <- factor(taxa_clusters$Feature, levels = taxa_clusters$Feature)

# Create annotation plot
row_annotation <- ggplot(taxa_clusters) + 
  geom_tile(aes(x = NA, y = Feature, fill = clusters)) + coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin=margin(0,0,0,0),
        ) +
      labs(fill = "Clusters", x = "Clusters")

# to view the notation, run
# row_annotation

# Add sample names to one of the columns
sample_data$sample <- factor(rownames(sample_data), levels = rownames(sample_data))

# Create annotation plot
sample_types_annotation <- ggplot(sample_data) + scale_y_discrete(position = "right", expand = c(0,0)) +
  geom_tile(aes(y = NA, x = sample, fill = sample_types)) + coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(0,0,0,0),
        axis.title.y.right = element_text(angle=0, vjust = 0.5)
        ) +
      labs(fill = "Sample types", y = "Sample types")
# to view the notation, run
# sample_types_annotation

# Create annotation plot
sample_clusters_annotation <- ggplot(sample_data) + scale_y_discrete(position = "right", expand = c(0,0)) +
  geom_tile(aes(y = NA, x = sample, fill = clusters)) + coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(0,0,0,0),
        axis.title.y.right = element_text(angle=0, vjust = 0.5)
        ) +
      labs(fill = "Clusters", y = "Clusters")
# to view the notation, run
# sample_clusters_annotation

# Order data based on clusters and sample types
mat <- mat[unfactor(taxa_clusters$Feature), unfactor(sample_data$sample)]

# ggplot requires data in melted format
melted_mat <- melt(mat)
colnames(melted_mat) <- c("Taxa", "Sample", "clr_z")

# Determines the scaling of colorss
maxval <- round(max(abs(melted_mat$clr_z)))
limits <- c(-maxval, maxval)
breaks <- seq(from = min(limits), to = max(limits), by = 0.5)
colours <- c("darkblue", "blue", "white", "red", "darkred")

heatmap <- ggplot(melted_mat) + 
  geom_tile(aes(x = Sample, y = Taxa, fill = clr_z)) +
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    
    plot.margin=margin(0,0,0,0), # removes margins
    legend.key.height= unit(1, 'cm')
    ) +
  scale_fill_gradientn(name = "CLR + Z transform", 
                       breaks = breaks, 
                       limits = limits, 
                       colours = colours) + 
  scale_y_discrete(position = "right")

heatmap
```

![](19_visualization_techniques_files/figure-latex/more_complex_heatmap-1.pdf)<!-- --> 


```r
# Create layout
design <- c(
  area(3, 1, 4, 1),
  area(1, 2, 1, 3),
  area(2, 2, 2, 3),
  area(3, 2, 4, 3)
)
# to view the design, run
# plot(design)

# Combine plots
plot <- row_annotation + sample_clusters_annotation + sample_types_annotation + heatmap  +
    plot_layout(design = design, guides = "collect", # Specify layout, collect legends
                
                # Adjust widths and heights to align plots.
                # When annotation plot is larger, it might not fit into its column/row.
                # Then you need to make column/row larger.
                
                # Relative widths and heights of each column and row:
                # Currently, the width of the first column is 15 % and the height of
                # first two rows are 30 % the size of others
                
                # To get this work most of the times, you can adjust all sizes to be 1, i.e. equal, 
                # but then the gaps between plots are larger.
                widths = c(0.15, 1, 1),
                heights = c(0.3, 0.3, 1, 1))

# plot
```


```r
# Create layout
design <- c(
  area(4, 1, 5, 1),
  area(4, 2, 5, 2),
  area(1, 3, 1, 4),
  area(2, 3, 2, 4),
  area(3, 3, 3, 4),
  area(4, 3, 5, 4)
)

# to view the design, run
# plot(design)

# Combine plots
plot <- taxa_tree + 
  row_annotation +
  sample_tree + 
  sample_clusters_annotation +
  sample_types_annotation +
  heatmap +
    plot_layout(design = design, guides = "collect", # Specify layout, collect legends
                widths = c(0.2, 0.15, 1, 1, 1),
                heights = c(0.1, 0.15, 0.15, 0.25, 1, 1))

plot
```

![](19_visualization_techniques_files/figure-latex/more_complex_heatmap3-1.pdf)<!-- --> 

Heatmaps find several other applications in biclustering and multi-assay
analyses, that are discussed in chapters \@ref(biclustering) and
\@ref(multi-assay_analyses), where the packages _cobiclust_ and _MOFA2_ are of interest.


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
 [1] ggpubr_0.4.0                   ggtree_3.4.0                  
 [3] ape_5.6-2                      pheatmap_1.0.12               
 [5] reshape2_1.4.4                 sechm_1.4.1                   
 [7] miaViz_1.3.3                   ggraph_2.0.5                  
 [9] patchwork_1.1.1                scater_1.24.0                 
[11] scuttle_1.6.2                  mia_1.3.27                    
[13] MultiAssayExperiment_1.22.0    TreeSummarizedExperiment_2.1.4
[15] Biostrings_2.64.0              XVector_0.36.0                
[17] SingleCellExperiment_1.18.0    SummarizedExperiment_1.26.1   
[19] Biobase_2.56.0                 GenomicRanges_1.48.0          
[21] GenomeInfoDb_1.32.2            IRanges_2.30.0                
[23] S4Vectors_0.34.0               BiocGenerics_0.42.0           
[25] MatrixGenerics_1.8.1           matrixStats_0.62.0-9000       
[27] ggplot2_3.3.6                  BiocStyle_2.24.0              
[29] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                  tidyselect_1.1.2           
  [3] RSQLite_2.2.14              grid_4.2.0                 
  [5] TSP_1.2-0                   BiocParallel_1.30.3        
  [7] Rtsne_0.16                  munsell_0.5.0              
  [9] ScaledMatrix_1.4.0          codetools_0.2-18           
 [11] withr_2.5.0                 colorspace_2.0-3           
 [13] filelock_1.0.2              highr_0.9                  
 [15] knitr_1.39                  ggsignif_0.6.3             
 [17] labeling_0.4.2              GenomeInfoDbData_1.2.8     
 [19] polyclip_1.10-0             bit64_4.0.5                
 [21] farver_2.1.0                vctrs_0.4.1                
 [23] treeio_1.20.0               generics_0.1.2             
 [25] xfun_0.31                   R6_2.5.1                   
 [27] doParallel_1.0.17           ggbeeswarm_0.6.0           
 [29] clue_0.3-61                 graphlayouts_0.8.0         
 [31] rsvd_1.0.5                  seriation_1.3.5            
 [33] bitops_1.0-7                cachem_1.0.6               
 [35] gridGraphics_0.5-1          DelayedArray_0.22.0        
 [37] assertthat_0.2.1            scales_1.2.0               
 [39] beeswarm_0.4.0              gtable_0.3.0               
 [41] beachmat_2.12.0             Cairo_1.5-15               
 [43] tidygraph_1.2.1             rlang_1.0.3                
 [45] GlobalOptions_0.1.2         splines_4.2.0              
 [47] rstatix_0.7.0               lazyeval_0.2.2             
 [49] broom_1.0.0                 BiocManager_1.30.18        
 [51] yaml_2.3.5                  abind_1.4-5                
 [53] backports_1.4.1             tools_4.2.0                
 [55] bookdown_0.27               ggplotify_0.1.0            
 [57] ellipsis_0.3.2              decontam_1.16.0            
 [59] RColorBrewer_1.1-3          Rcpp_1.0.8.3               
 [61] plyr_1.8.7                  sparseMatrixStats_1.8.0    
 [63] zlibbioc_1.42.0             purrr_0.3.4                
 [65] RCurl_1.98-1.7              GetoptLong_1.0.5           
 [67] viridis_0.6.2               cowplot_1.1.1              
 [69] ggrepel_0.9.1               cluster_2.1.3              
 [71] DECIPHER_2.24.0             magrittr_2.0.3             
 [73] circlize_0.4.15             ggnewscale_0.4.7           
 [75] randomcoloR_1.1.0.1         evaluate_0.15              
 [77] XML_3.99-0.10               gridExtra_2.3              
 [79] shape_1.4.6                 compiler_4.2.0             
 [81] tibble_3.1.7                V8_4.2.0                   
 [83] crayon_1.5.1                htmltools_0.5.2            
 [85] ggfun_0.0.6                 mgcv_1.8-40                
 [87] tidyr_1.2.0                 aplot_0.1.6                
 [89] DBI_1.1.3                   tweenr_1.0.2               
 [91] ComplexHeatmap_2.12.0       MASS_7.3-57                
 [93] Matrix_1.4-1                car_3.1-0                  
 [95] permute_0.9-7               cli_3.3.0                  
 [97] parallel_4.2.0              igraph_1.3.2               
 [99] pkgconfig_2.0.3             dir.expiry_1.4.0           
[101] registry_0.5-1              foreach_1.5.2              
[103] vipor_0.4.5                 DirichletMultinomial_1.38.0
[105] yulab.utils_0.0.5           stringr_1.4.0              
[107] digest_0.6.29               vegan_2.6-2                
[109] graph_1.74.0                rmarkdown_2.14             
[111] tidytree_0.3.9              DelayedMatrixStats_1.18.0  
[113] curl_4.3.2                  rjson_0.2.21               
[115] lifecycle_1.0.1             nlme_3.1-158               
[117] jsonlite_1.8.0              carData_3.0-5              
[119] BiocNeighbors_1.14.0        CodeDepends_0.6.5          
[121] viridisLite_0.4.0           fansi_1.0.3                
[123] pillar_1.7.0                lattice_0.20-45            
[125] fastmap_1.1.0               glue_1.6.2                 
[127] png_0.1-7                   iterators_1.0.14           
[129] bit_4.0.4                   ggforce_0.3.3              
[131] stringi_1.7.6               blob_1.2.3                 
[133] BiocSingular_1.12.0         memoise_2.0.1              
[135] dplyr_1.0.9                 irlba_2.3.5                
```
</div>
