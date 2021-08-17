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

In addition to _pheatmap_ package, there are also other packages that provide 
functions for more complex heatmaps. One example is _iheatmapr_ package. Examples
of using it you can find from its 
[vignette](https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html).



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
 [1] dendextend_1.15.1              pheatmap_1.0.12               
 [3] miaViz_1.1.4                   ggraph_2.0.5                  
 [5] ggplot2_3.3.5                  mia_1.1.9                     
 [7] TreeSummarizedExperiment_2.1.4 Biostrings_2.61.2             
 [9] XVector_0.33.0                 SingleCellExperiment_1.15.1   
[11] SummarizedExperiment_1.23.1    Biobase_2.53.0                
[13] GenomicRanges_1.45.0           GenomeInfoDb_1.29.3           
[15] IRanges_2.27.0                 S4Vectors_0.31.0              
[17] BiocGenerics_0.39.1            MatrixGenerics_1.5.3          
[19] matrixStats_0.60.0             ecodist_2.0.7                 
[21] BiocStyle_2.21.3               rebook_1.3.0                  

loaded via a namespace (and not attached):
  [1] plyr_1.8.6                  igraph_1.2.6               
  [3] lazyeval_0.2.2              splines_4.1.0              
  [5] BiocParallel_1.27.3         scater_1.21.3              
  [7] digest_0.6.27               htmltools_0.5.1.1          
  [9] viridis_0.6.1               fansi_0.5.0                
 [11] magrittr_2.0.1              memoise_2.0.0              
 [13] ScaledMatrix_1.1.0          cluster_2.1.2              
 [15] DECIPHER_2.21.0             graphlayouts_0.7.1         
 [17] colorspace_2.0-2            blob_1.2.2                 
 [19] ggrepel_0.9.1               xfun_0.25                  
 [21] dplyr_1.0.7                 crayon_1.4.1               
 [23] RCurl_1.98-1.3              jsonlite_1.7.2             
 [25] graph_1.71.2                ape_5.5                    
 [27] glue_1.4.2                  polyclip_1.10-0            
 [29] gtable_0.3.0                zlibbioc_1.39.0            
 [31] DelayedArray_0.19.1         BiocSingular_1.9.1         
 [33] scales_1.1.1                DBI_1.1.1                  
 [35] Rcpp_1.0.7                  viridisLite_0.4.0          
 [37] decontam_1.13.0             tidytree_0.3.4             
 [39] bit_4.0.4                   rsvd_1.0.5                 
 [41] dir.expiry_1.1.0            RColorBrewer_1.1-2         
 [43] ellipsis_0.3.2              pkgconfig_2.0.3            
 [45] XML_3.99-0.6                farver_2.1.0               
 [47] scuttle_1.3.1               CodeDepends_0.6.5          
 [49] sass_0.4.0                  utf8_1.2.2                 
 [51] tidyselect_1.1.1            labeling_0.4.2             
 [53] rlang_0.4.11                reshape2_1.4.4             
 [55] munsell_0.5.0               tools_4.1.0                
 [57] cachem_1.0.5                DirichletMultinomial_1.35.0
 [59] generics_0.1.0              RSQLite_2.2.7              
 [61] evaluate_0.14               stringr_1.4.0              
 [63] fastmap_1.1.0               yaml_2.2.1                 
 [65] ggtree_3.1.3                knitr_1.33                 
 [67] bit64_4.0.5                 tidygraph_1.2.0            
 [69] purrr_0.3.4                 nlme_3.1-152               
 [71] sparseMatrixStats_1.5.2     aplot_0.0.6                
 [73] compiler_4.1.0              beeswarm_0.4.0             
 [75] filelock_1.0.2              treeio_1.17.2              
 [77] tibble_3.1.3                tweenr_1.0.2               
 [79] bslib_0.2.5.1               stringi_1.7.3              
 [81] highr_0.9                   lattice_0.20-44            
 [83] Matrix_1.3-4                vegan_2.5-7                
 [85] permute_0.9-5               vctrs_0.3.8                
 [87] pillar_1.6.2                lifecycle_1.0.0            
 [89] BiocManager_1.30.16         jquerylib_0.1.4            
 [91] BiocNeighbors_1.11.0        bitops_1.0-7               
 [93] irlba_2.3.3                 patchwork_1.1.1            
 [95] R6_2.5.0                    bookdown_0.22              
 [97] gridExtra_2.3               vipor_0.4.5                
 [99] codetools_0.2-18            MASS_7.3-54                
[101] assertthat_0.2.1            withr_2.4.2                
[103] GenomeInfoDbData_1.2.6      mgcv_1.8-36                
[105] parallel_4.1.0              grid_4.1.0                 
[107] beachmat_2.9.1              tidyr_1.1.3                
[109] rmarkdown_2.10              DelayedMatrixStats_1.15.2  
[111] rvcheck_0.1.8               ggnewscale_0.4.5           
[113] ggforce_0.3.3               ggbeeswarm_0.6.0           
```
</div>

