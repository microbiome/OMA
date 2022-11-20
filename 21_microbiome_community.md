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
tse_phylum <- transformSamples(tse_phylum, assay_name = "relabundance", method = "clr")
# Add z-transformation on features (taxa)
tse_phylum <- transformFeatures(tse_phylum, assay_name = "clr", 
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

![](21_microbiome_community_files/figure-latex/heatmap-1.pdf)<!-- --> 

In addition, there are also other packages that provide functions for more complex heatmaps,
such as [_iheatmapr_](https://docs.ropensci.org/iheatmapr/articles/full_vignettes/iheatmapr.html)
and [ComplexHeatmap](https://academic.oup.com/bioinformatics/article/32/18/2847/1743594?login=true).
[sechm](http://www.bioconductor.org/packages/release/bioc/vignettes/sechm/inst/doc/sechm.html)
package provides wrapper for _ComplexHeatmap_ and its usage is explained in chapter \@ref(viz-chapter)
along with the `pheatmap` package for clustered heatmaps.

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.2.1 (2022-06-23)
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
 [1] miaViz_1.5.4                   ggraph_2.1.0                  
 [3] ggplot2_3.4.0                  mia_1.5.17                    
 [5] MultiAssayExperiment_1.24.0    TreeSummarizedExperiment_2.1.4
 [7] Biostrings_2.66.0              XVector_0.38.0                
 [9] SingleCellExperiment_1.20.0    SummarizedExperiment_1.28.0   
[11] Biobase_2.58.0                 GenomicRanges_1.50.1          
[13] GenomeInfoDb_1.34.3            IRanges_2.32.0                
[15] S4Vectors_0.36.0               BiocGenerics_0.44.0           
[17] MatrixGenerics_1.10.0          matrixStats_0.63.0-9003       
[19] BiocStyle_2.24.0               rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] ggtree_3.4.4                ggnewscale_0.4.8           
  [3] ggbeeswarm_0.6.0            colorspace_2.0-3           
  [5] ellipsis_0.3.2              scuttle_1.8.0              
  [7] BiocNeighbors_1.16.0        aplot_0.1.8                
  [9] farver_2.1.1                graphlayouts_0.8.3         
 [11] ggrepel_0.9.2               bit64_4.0.5                
 [13] fansi_1.0.3                 decontam_1.18.0            
 [15] codetools_0.2-18            splines_4.2.1              
 [17] sparseMatrixStats_1.10.0    cachem_1.0.6               
 [19] knitr_1.40                  scater_1.26.1              
 [21] polyclip_1.10-4             jsonlite_1.8.3             
 [23] cluster_2.1.4               graph_1.74.0               
 [25] ggforce_0.4.1               BiocManager_1.30.19        
 [27] compiler_4.2.1              assertthat_0.2.1           
 [29] Matrix_1.5-3                fastmap_1.1.0              
 [31] lazyeval_0.2.2              cli_3.4.1                  
 [33] tweenr_2.0.2                BiocSingular_1.14.0        
 [35] htmltools_0.5.3             tools_4.2.1                
 [37] igraph_1.3.5                rsvd_1.0.5                 
 [39] gtable_0.3.1                glue_1.6.2                 
 [41] GenomeInfoDbData_1.2.9      reshape2_1.4.4             
 [43] dplyr_1.0.10                Rcpp_1.0.9                 
 [45] vctrs_0.5.1                 ape_5.6-2                  
 [47] nlme_3.1-160                DECIPHER_2.26.0            
 [49] DelayedMatrixStats_1.20.0   xfun_0.35                  
 [51] stringr_1.4.1               beachmat_2.14.0            
 [53] lifecycle_1.0.3             irlba_2.3.5.1              
 [55] XML_3.99-0.12               zlibbioc_1.44.0            
 [57] MASS_7.3-58.1               scales_1.2.1               
 [59] tidygraph_1.2.2             parallel_4.2.1             
 [61] yaml_2.3.6                  memoise_2.0.1              
 [63] gridExtra_2.3               ggfun_0.0.8                
 [65] yulab.utils_0.0.5           stringi_1.7.8              
 [67] RSQLite_2.2.18              highr_0.9                  
 [69] ScaledMatrix_1.6.0          tidytree_0.4.1             
 [71] permute_0.9-7               filelock_1.0.2             
 [73] BiocParallel_1.32.1         rlang_1.0.6                
 [75] pkgconfig_2.0.3             bitops_1.0-7               
 [77] evaluate_0.18               lattice_0.20-45            
 [79] purrr_0.3.5                 labeling_0.4.2             
 [81] patchwork_1.1.2             treeio_1.22.0              
 [83] CodeDepends_0.6.5           bit_4.0.5                  
 [85] tidyselect_1.2.0            plyr_1.8.8                 
 [87] magrittr_2.0.3              bookdown_0.30              
 [89] R6_2.5.1                    generics_0.1.3             
 [91] DelayedArray_0.24.0         DBI_1.1.3                  
 [93] withr_2.5.0                 pillar_1.8.1               
 [95] mgcv_1.8-41                 RCurl_1.98-1.9             
 [97] tibble_3.1.8                dir.expiry_1.4.0           
 [99] crayon_1.5.2                utf8_1.2.2                 
[101] rmarkdown_2.18              viridis_0.6.2              
[103] grid_4.2.1                  blob_1.2.3                 
[105] vegan_2.6-4                 digest_0.6.30              
[107] tidyr_1.2.1                 gridGraphics_0.5-1         
[109] munsell_0.5.0               DirichletMultinomial_1.40.0
[111] ggplotify_0.1.0             beeswarm_0.4.0             
[113] viridisLite_0.4.1           vipor_0.4.5                
```
</div>

