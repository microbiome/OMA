# Biclustering {#biclustering}

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

Biclustering methods cluster rows and columns simultaneously in order
to find subsets of correlated features/samples.

Here, we use following packages:
-   [biclust](https://cran.r-project.org/web/packages/biclust/index.html)
-   [cobiclust](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13582)

_cobiclust_ is especially developed for microbiome data whereas _biclust_ is more
general method. In this section, we show three different cases and example 
solutions to apply biclustering to them. 

1.   Taxa vs samples
2.   Taxa vs biomolecule/biomarker
3.   Taxa vs taxa

Biclusters can be visualized using heatmap or boxplot, for instance. For checking purposes, 
also scatter plot might be valid choice.

Check more ideas for heatmaps from chapters \@ref(viz-chapter) and \@ref(microbiome-community.

## Taxa vs samples

When you have microbial abundance matrices, we suggest to use _cobiclust_ which is
designed for microbial data. 

Load example data

```r
library(mia)
data("HintikkaXOData")
mae <- HintikkaXOData
```

Only the most prevalent taxa are included in analysis. 


```r
# Subset data in the first experiment
mae[[1]] <- subsetByPrevalentTaxa(mae[[1]], rank = "Genus", prevalence = 0.2, detection = 0.001)
# clr-transform in the first experiment
mae[[1]] <- transformSamples(mae[[1]], method = "relabundance", pseudocount = 1)
mae[[1]] <- transformSamples(mae[[1]], "relabundance", method = "clr")
```

_cobiclust_ takes counts table as an input and gives _cobiclust_ object as an output.
It includes clusters for taxa and samples. 


```r
if(!require(cobiclust)){
    install.packages("cobiclust")
    library(cobiclust)
}

# Do clustering; use counts table´
clusters <- cobiclust(assay(mae[[1]], "counts"))

# Get clusters
row_clusters <- clusters$classification$rowclass
col_clusters <- clusters$classification$colclass

# Add clusters to rowdata and coldata
rowData(mae[[1]])$clusters <- factor(row_clusters)
colData(mae[[1]])$clusters <- factor(col_clusters)

# Order data based on clusters
mae[[1]] <- mae[[1]][order(rowData(mae[[1]])$clusters), order(colData(mae[[1]])$clusters)]

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

Next we can plot clusters. Commonly used plot is heatmap with annotations. 


```r
if(!require(pheatmap)){
    install.packages("pheatmap")
    library(pheatmap)
}
# z-transform for heatmap
mae[[1]] <- transformFeatures(mae[[1]], assay_name = "clr", method = "z", name = "clr_z")

# Create annotations. When column names are equal, they should share levels. 
# Here samples include 3 clusters, and taxa 2. That is why we have to make 
# column names unique. 
annotation_col <- data.frame(colData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_col) <- "col_clusters"

annotation_row <- data.frame(rowData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_row) <- "row_clusters"

# Create a heatmap
pheatmap(assay(mae[[1]], "clr_z"), cluster_rows = F, cluster_cols = F, 
         annotation_col = annotation_col,
         annotation_row = annotation_row)
```

![](24_biclustering_files/figure-latex/cobiclust_3-1.pdf)<!-- --> 

Boxplot is commonly used to summarize the results:


```r
if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
}
if(!require(patchwork)){
    install.packages("patchwork")
    library(patchwork)
}

# ggplot requires data in melted format
melt_assay <- meltAssay(mae[[1]], assay_name = "clr", add_col_data = T, add_row_data = T)

# patchwork two plots side-by-side
p1 <- ggplot(melt_assay) +
  geom_boxplot(aes(x = clusters.x, y = clr)) +
  labs(x = "Taxa clusters")

p2 <- ggplot(melt_assay) +
  geom_boxplot(aes(x = clusters.y, y = clr)) +
  labs(x = "Sample clusters")

p1 + p2
```

![](24_biclustering_files/figure-latex/cobiclust_4-1.pdf)<!-- --> 

## Taxa vs biomolecules

Here, we analyze cross-correlation between taxa and metabolites. This is a case, where
we use _biclust_ method which is suitable for numeric matrices in general.


```r
# Samples must be in equal order 
# (Only 1st experiment  was ordered in cobiclust step leading to unequal order)
mae[[1]] <- mae[[1]][ , colnames(mae[[2]]) ]

# Make rownames unique since it is require by other steps
rownames(mae[[1]]) <- make.unique(rownames(mae[[1]]))
# Calculate correlations
corr <- getExperimentCrossCorrelation(mae, 1, 2, 
                                      assay_name1 = "clr", 
                                      assay_name2 = "nmr", 
                                      mode = "matrix", 
                                      cor_threshold = 0.2)
```


_biclust_ takes matrix as an input and returns _biclust_ object. 


```r
# Load package
if(!require(biclust)){
    install.packages("biclust")
    library(biclust)
}

# Set seed for reproducibility
set.seed(3973)

# Find biclusters
bc <- biclust(corr, method=BCPlaid(), fit.model = y ~ m,
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10,
              iter.startup = 10, iter.layer = 100, verbose = FALSE)

bc
```

```
## 
## An object of class Biclust 
## 
## call:
## 	biclust(x = corr, method = BCPlaid(), fit.model = y ~ m, background = TRUE, 
## 	    shuffle = 100, back.fit = 0, max.layers = 10, iter.startup = 10, 
## 	    iter.layer = 100, verbose = FALSE)
## 
## There was no cluster found
```

The object includes cluster information. However compared to _cobiclust_, 
_biclust_ object includes only information about clusters that were found, not general cluster.

Meaning that if one cluster size of 5 features was found out of 20 features, 
those 15 features do not belong to any cluster. That is why we have to create an
additional cluster for features/samples that are not assigned into any cluster.


```r
# Functions for obtaining biclust information

# Get clusters for rows and columns
.get_biclusters_from_biclust <- function(bc, assay){
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

# Input clusters, and how many observations there should be, i.e., the number of samples or features
.manipulate_bc_data <- function(bc_clusters, assay, row_col){
  # Get right dimension
  dim <- ifelse(row_col == "col", ncol(assay), nrow(assay))
  # Get column/row names
  if( row_col == "col" ){
    names <- colnames(assay)
  } else{
    names <- rownames(assay)
  }
  
  # If no clusters were found, create one. Otherwise create additional cluster which
  # contain those samples that are not included in clusters that were found.
  if( nrow(bc_clusters) != dim ){
      bc_clusters <- data.frame(cluster = rep(TRUE, dim))
  } else {
      # Create additional cluster that includes those samples/features that
      # are not included in other clusters.
      vec <- ifelse(rowSums(bc_clusters) > 0, FALSE, TRUE)
      # If additional cluster contains samples, then add it
      if ( any(vec) ){
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
##                          cluster_1
## D_5__Ruminiclostridium 5      TRUE
## D_5__uncultured               TRUE
## D_5__Lactococcus              TRUE
## D_5__Lachnoclostridium        TRUE
## D_5__Holdemania               TRUE
## D_5__Anaerostipes             TRUE
```

Let's collect information for the scatter plot. 


```r
# Function for obtaining sample-wise sum, mean, median, and mean variance for each cluster
.sum_mean_median_var <- function(tse1, tse2, assay_name1, assay_name2, clusters1, clusters2){
  
  list <- list()
  # Create a data frame that includes all the information
  for(i in 1:ncol(clusters1) ){
    # Subset data based on cluster
    tse_subset1 <- tse1[clusters1[,i], ]
    tse_subset2 <- tse2[clusters2[,i], ]
    # Get assay
    assay1 <- assay(tse_subset1, assay_name1)
    assay2 <- assay(tse_subset2, assay_name2)
    # Calculate sum, mean, median, and mean variance
    sum1 <- colSums2(assay1, na.rm = T)
    mean1 <- colMeans2(assay1, na.rm = T)
    median1 <- colMedians(assay1, na.rm = T)
    var1 <- colVars(assay1, na.rm = T)
    
    sum2 <- colSums2(assay2, na.rm = T)
    mean2 <- colMeans2(assay2, na.rm = T)
    median2 <- colMedians(assay2, na.rm = T)
    var2 <- colVars(assay2, na.rm = T)
    
    list[[i]] <- data.frame(sample = colnames(tse1), sum1, sum2, mean1, mean2, 
                     median1, median2, var1, var2)
  }

  return(list)
}

# Calculate info
df <- .sum_mean_median_var(mae[[1]], mae[[2]], "clr", "nmr", bicluster_rows, bicluster_columns)
```

Now we can create a scatter plot. X-axis includes median clr abundance of microbiome
and y-axis median absolute concentration of each metabolite. Each data point represents
a single sample. 

From the plots, we can see that there is low negative correlation in both cluster 1 and 3.
This means that when abundance of bacteria belonging to cluster 1 or 3 is higher, 
the concentration of metabolites of cluster 1 or 3 is lower, and vice versa. 


```r
pics <- list()
for(i in seq_along(df)){
  pics[[i]] <- ggplot(df[[i]])  +
      geom_point(aes(x = median1, y = median2)) + 
      labs(title = paste0("Cluster ", i),
           x = "Taxa (clr median)",
           y = "Metabolites (abs. median)")
  print(pics[[i]])
}
```


\includegraphics[width=0.33\linewidth]{24_biclustering_files/figure-latex/biclust_6-1} 

```r
# pics[[1]] + pics[[2]] + pics[[3]]
```

_pheatmap_ does not allow boolean values, so they must be converted into factors.


```r
bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))
```

Again, we can plot clusters with heatmap.


```r
# Adjust colors for all clusters
if( ncol(bicluster_rows) > ncol(bicluster_columns) ){
  cluster_names <- colnames(bicluster_rows)
} else {
  cluster_names <- colnames(bicluster_columns)
}
annotation_colors <- list()
for(name in cluster_names){
  annotation_colors[[name]] <- c("TRUE" = "red", "FALSE" = "white")
}

# Create a heatmap
pheatmap(corr, cluster_cols = F, cluster_rows = F,
         annotation_col = bicluster_columns, 
         annotation_row = bicluster_rows,
         annotation_colors = annotation_colors)
```

![](24_biclustering_files/figure-latex/biclust_8-1.pdf)<!-- --> 

## Taxa vs taxa

Third and final example deals with situation where we want to analyze correlation
between taxa. _biclust_ is suitable for this. 


```r
# Calculate cross-correlation
corr <- getExperimentCrossCorrelation(mae, 1, 1, 
                                      assay_name1 = "clr", assay_name2 = "clr", 
                                      mode = "matrix",
                                      cor_threshold = 0.2, verbose = F, show_warning = F)

# Find biclusters
bc <- biclust(corr, method=BCPlaid(), fit.model = y ~ m,
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10,
              iter.startup = 10, iter.layer = 100, verbose = FALSE)
```


```r
# Get biclusters
bcs <- .get_biclusters_from_biclust(bc, corr)

bicluster_rows <- bcs$bc_rows
bicluster_columns <- bcs$bc_columns
```



```r
# Create a column that combines information
# If row/column includes in multiple clusters, cluster numbers are separated with "_&_"
bicluster_columns$clusters <- apply(bicluster_columns, 1, 
                                    function(x){paste(paste(which(x)), collapse = "_&_") })
bicluster_columns <- bicluster_columns[, "clusters", drop = FALSE]

bicluster_rows$clusters <- apply(bicluster_rows, 1, 
                                 function(x){paste(paste(which(x)), collapse = "_&_") })
bicluster_rows <- bicluster_rows[, "clusters", drop = FALSE]
```


```r
# Convert boolean values into factor
bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))

pheatmap(corr, cluster_cols = F, cluster_rows = F,
         annotation_col = bicluster_columns, 
         annotation_row = bicluster_rows)
```

![](24_biclustering_files/figure-latex/biclust_12-1.pdf)<!-- --> 

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
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] biclust_2.0.3                  lattice_0.20-45               
 [3] colorspace_2.0-3               MASS_7.3-58.1                 
 [5] patchwork_1.1.2                ggplot2_3.3.6                 
 [7] pheatmap_1.0.12                cobiclust_0.1.0               
 [9] mia_1.5.12                     MultiAssayExperiment_1.22.0   
[11] TreeSummarizedExperiment_2.1.4 Biostrings_2.64.1             
[13] XVector_0.36.0                 SingleCellExperiment_1.18.0   
[15] SummarizedExperiment_1.26.1    Biobase_2.56.0                
[17] GenomicRanges_1.48.0           GenomeInfoDb_1.32.3           
[19] IRanges_2.30.1                 S4Vectors_0.34.0              
[21] BiocGenerics_0.42.0            MatrixGenerics_1.8.1          
[23] matrixStats_0.62.0-9003        BiocStyle_2.24.0              
[25] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] ggbeeswarm_0.6.0            class_7.3-20               
  [3] modeltools_0.2-23           ellipsis_0.3.2             
  [5] scuttle_1.6.3               BiocNeighbors_1.14.0       
  [7] farver_2.1.1                ggrepel_0.9.1              
  [9] bit64_4.0.5                 fansi_1.0.3                
 [11] decontam_1.16.0             codetools_0.2-18           
 [13] splines_4.2.1               sparseMatrixStats_1.8.0    
 [15] cachem_1.0.6                knitr_1.40                 
 [17] scater_1.24.0               jsonlite_1.8.0             
 [19] cluster_2.1.4               graph_1.74.0               
 [21] BiocManager_1.30.18         compiler_4.2.1             
 [23] assertthat_0.2.1            Matrix_1.4-1               
 [25] fastmap_1.1.0               lazyeval_0.2.2             
 [27] cli_3.3.0                   BiocSingular_1.12.0        
 [29] htmltools_0.5.3             tools_4.2.1                
 [31] rsvd_1.0.5                  gtable_0.3.1               
 [33] glue_1.6.2                  GenomeInfoDbData_1.2.8     
 [35] reshape2_1.4.4              dplyr_1.0.10               
 [37] Rcpp_1.0.9                  vctrs_0.4.1                
 [39] ape_5.6-2                   nlme_3.1-159               
 [41] DECIPHER_2.24.0             DelayedMatrixStats_1.18.0  
 [43] xfun_0.32                   stringr_1.4.1              
 [45] beachmat_2.12.0             lifecycle_1.0.1            
 [47] irlba_2.3.5                 XML_3.99-0.10              
 [49] zlibbioc_1.42.0             scales_1.2.1               
 [51] parallel_4.2.1              RColorBrewer_1.1-3         
 [53] yaml_2.3.5                  memoise_2.0.1              
 [55] gridExtra_2.3               yulab.utils_0.0.5          
 [57] stringi_1.7.8               RSQLite_2.2.16             
 [59] highr_0.9                   ScaledMatrix_1.4.0         
 [61] tidytree_0.4.0              permute_0.9-7              
 [63] filelock_1.0.2              BiocParallel_1.30.3        
 [65] rlang_1.0.5                 pkgconfig_2.0.3            
 [67] bitops_1.0-7                evaluate_0.16              
 [69] purrr_0.3.4                 additivityTests_1.1-4.1    
 [71] labeling_0.4.2              treeio_1.20.2              
 [73] CodeDepends_0.6.5           bit_4.0.4                  
 [75] tidyselect_1.1.2            plyr_1.8.7                 
 [77] magrittr_2.0.3              bookdown_0.28              
 [79] R6_2.5.1                    generics_0.1.3             
 [81] DelayedArray_0.22.0         DBI_1.1.3                  
 [83] withr_2.5.0                 pillar_1.8.1               
 [85] mgcv_1.8-40                 RCurl_1.98-1.8             
 [87] tibble_3.1.8                dir.expiry_1.4.0           
 [89] crayon_1.5.1                utf8_1.2.2                 
 [91] rmarkdown_2.16              viridis_0.6.2              
 [93] flexclust_1.4-1             blob_1.2.3                 
 [95] vegan_2.6-2                 digest_0.6.29              
 [97] tidyr_1.2.0                 munsell_0.5.0              
 [99] DirichletMultinomial_1.38.0 beeswarm_0.4.0             
[101] viridisLite_0.4.1           vipor_0.4.5                
```
</div>
