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
top_taxa <- getTopTaxa(tse_phylum,top = 5, assay_name = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$Phylum,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse)$Phylum <- as.character(phylum_renamed)

# Visualizing the composition barplot, with samples order by "Bacteroidetes"
plotAbundance(tse, assay_name="relabundance", rank = "Phylum",
              order_rank_by="abund", 
              order_sample_by = "Bacteroidetes")
```

![](21_microbiome_community_files/figure-latex/unnamed-chunk-1-1.pdf)<!-- --> 

### Composition heatmap {#composition-heatmap}

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
assay(tse_phylum, "pseudo") <- assay(tse_phylum, "counts") + 1
tse_phylum <- transformCounts(tse_phylum, assay_name = "pseudo",
                              method = "relabundance")

tse_phylum <- transformCounts(tse_phylum,
                  assay_name = "relabundance",
		  method = "clr")

# Add z-transformation on features (taxa)
tse_phylum <- transformCounts(tse_phylum, assay_name = "clr", 
                              MARGIN = "features",
                              method = "z", name = "clr_z")
```

Visualize as heatmap.


```r
# Melt the assay for plotting purposes
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

![](21_microbiome_community_files/figure-latex/heatmapvisu-1.pdf)<!-- --> 


_pheatmap_ is a package that provides methods to plot clustered heatmaps. 


```r
if(!require(pheatmap)){install.packages("pheatmap"); library(pheatmap)}

# Takes subset: only samples from feces, skin, or tongue
tse_phylum_subset <- tse_phylum[ , colData(tse_phylum)$SampleType %in% c("Feces", "Skin", "Tongue") ]

# Add clr-transformation
tse_phylum_subset <- transformCounts(tse_phylum_subset,
                         method = "clr",
    			 pseudocount = 1)

tse_phylum_subset <- transformCounts(tse_phylum_subset, assay_name = "clr",
                                     MARGIN = "features", 
                                     method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopTaxa(tse_phylum_subset, top = 20)
tse_phylum_subset <- tse_phylum_subset[top_taxa, ]

# Gets the assay table
mat <- assay(tse_phylum_subset, "clr_z")

# Creates the heatmap
pheatmap(mat)
```

![](21_microbiome_community_files/figure-latex/pheatmap1-1.pdf)<!-- --> 

We can create clusters by hierarchical clustering and add them to the plot.


```r
if(!require(ape)){
    install.packages("ape")
    library(ape)
}

# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)
```


```r
if(!require(ggtree)){
    install.packages("ggtree")
    library(ggtree)
}

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

taxa_tree
```

![](21_microbiome_community_files/figure-latex/pheatmap3-1.pdf)<!-- --> 

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


```r
# Adds information to rowData
rowData(tse_phylum_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_phylum_subset))), ]

# Prints taxa and their clusters
rowData(tse_phylum_subset)$clusters
```

```
##  [1] 1 1 2 3 2 2 1 1 1 1 3 2 3 3 3 2 2 3 3 1
## Levels: 1 2 3
```


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

sample_tree
```

![](21_microbiome_community_files/figure-latex/pheatmap6-1.pdf)<!-- --> 


```r
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

![](21_microbiome_community_files/figure-latex/pheatmap8-1.pdf)<!-- --> 

In addition, there are also other packages that provide functions for
more complex heatmaps, such as
[_iheatmapr_](https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html)
and ComplexHeatmap [@ComplexHeatmap].
[sechm](http://www.bioconductor.org/packages/release/bioc/vignettes/sechm/inst/doc/sechm.html)
package provides wrapper for _ComplexHeatmap_ and its usage is
explained in chapter \@ref(viz-chapter) along with the `pheatmap`
package for clustered heatmaps.

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
 [1] ggtree_3.8.0                   ape_5.7-1                     
 [3] pheatmap_1.0.12                miaViz_1.7.5                  
 [5] ggraph_2.1.0                   ggplot2_3.4.2                 
 [7] mia_1.7.11                     MultiAssayExperiment_1.26.0   
 [9] TreeSummarizedExperiment_2.1.4 Biostrings_2.68.0             
[11] XVector_0.40.0                 SingleCellExperiment_1.22.0   
[13] SummarizedExperiment_1.30.0    Biobase_2.60.0                
[15] GenomicRanges_1.52.0           GenomeInfoDb_1.36.0           
[17] IRanges_2.34.0                 S4Vectors_0.38.0              
[19] BiocGenerics_0.46.0            MatrixGenerics_1.12.0         
[21] matrixStats_0.63.0-9003        BiocStyle_2.28.0              
[23] rebook_1.9.0                  

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3          jsonlite_1.8.4             
  [3] CodeDepends_0.6.5           magrittr_2.0.3             
  [5] ggbeeswarm_0.7.1            farver_2.1.1               
  [7] rmarkdown_2.21              zlibbioc_1.46.0            
  [9] vctrs_0.6.2                 memoise_2.0.1              
 [11] DelayedMatrixStats_1.22.0   RCurl_1.98-1.12            
 [13] htmltools_0.5.5             BiocNeighbors_1.18.0       
 [15] gridGraphics_0.5-1          plyr_1.8.8                 
 [17] DECIPHER_2.28.0             cachem_1.0.7               
 [19] igraph_1.4.2                lifecycle_1.0.3            
 [21] pkgconfig_2.0.3             rsvd_1.0.5                 
 [23] Matrix_1.5-4                R6_2.5.1                   
 [25] fastmap_1.1.1               GenomeInfoDbData_1.2.10    
 [27] digest_0.6.31               aplot_0.1.10               
 [29] colorspace_2.1-0            ggnewscale_0.4.8           
 [31] patchwork_1.1.2             scater_1.28.0              
 [33] irlba_2.3.5.1               RSQLite_2.3.1              
 [35] vegan_2.6-4                 beachmat_2.16.0            
 [37] labeling_0.4.2              filelock_1.0.2             
 [39] fansi_1.0.4                 polyclip_1.10-4            
 [41] mgcv_1.8-42                 compiler_4.3.0             
 [43] bit64_4.0.5                 withr_2.5.0                
 [45] BiocParallel_1.34.0         viridis_0.6.2              
 [47] DBI_1.1.3                   highr_0.10                 
 [49] ggforce_0.4.1               MASS_7.3-59                
 [51] DelayedArray_0.25.0         permute_0.9-7              
 [53] tools_4.3.0                 vipor_0.4.5                
 [55] beeswarm_0.4.0              glue_1.6.2                 
 [57] nlme_3.1-162                grid_4.3.0                 
 [59] cluster_2.1.4               reshape2_1.4.4             
 [61] generics_0.1.3              gtable_0.3.3               
 [63] tidyr_1.3.0                 BiocSingular_1.16.0        
 [65] tidygraph_1.2.3             ScaledMatrix_1.7.1         
 [67] utf8_1.2.3                  ggrepel_0.9.3              
 [69] pillar_1.9.0                stringr_1.5.0              
 [71] yulab.utils_0.0.6           splines_4.3.0              
 [73] dplyr_1.1.2                 tweenr_2.0.2               
 [75] treeio_1.24.0               lattice_0.21-8             
 [77] bit_4.0.5                   tidyselect_1.2.0           
 [79] DirichletMultinomial_1.42.0 scuttle_1.10.0             
 [81] knitr_1.42                  gridExtra_2.3              
 [83] bookdown_0.33               xfun_0.39                  
 [85] graphlayouts_0.8.4          stringi_1.7.12             
 [87] lazyeval_0.2.2              ggfun_0.0.9                
 [89] yaml_2.3.7                  evaluate_0.20              
 [91] codetools_0.2-19            tibble_3.2.1               
 [93] BiocManager_1.30.20         graph_1.78.0               
 [95] ggplotify_0.1.0             cli_3.6.1                  
 [97] munsell_0.5.0               Rcpp_1.0.10                
 [99] dir.expiry_1.8.0            XML_3.99-0.14              
[101] parallel_4.3.0              blob_1.2.4                 
[103] sparseMatrixStats_1.12.0    bitops_1.0-7               
[105] decontam_1.20.0             viridisLite_0.4.1          
[107] tidytree_0.4.2              scales_1.2.1               
[109] purrr_1.0.1                 crayon_1.5.2               
[111] rlang_1.1.0                
```
</div>

