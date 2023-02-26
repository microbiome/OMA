# Multi-assay analyses {#multi-assay-analyses}

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
```

Multi-omics means that we integrate data from multiple sources. For
example, we can integrate microbial abundances in the gut with
biomolecular profiling data from blood samples. This kind of
integrative multi-omic approaches can support the analysis of
microbiome dysbiosis and facilitate the discovery of novel biomarkers
for health and disease.

With cross-correlation analysis, we can analyze how strongly and how
differently variables are associated between each other. For instance,
we can analyze if higher presence of a specific taxon equals to higher
levels of a biomolecule.

The data containers that the _miaverse_ utilizes are scalable and they
can contain different types of data in a same container. Because of
that, the _miaverse_ is well-suitable for multi-assay microbiome data
which incorporates different types of complementary data sources in a
single reproducible workflow.

Another experiment can be stored in altExp slot of SE data container
or both experiments can be stored side-by-side in MAE data container
(see the sections \@ref(alt-exp) and \@ref(mae) to learn more about
altExp and MAE objects, respectively).

Different experiments are first imported into SE or TreeSE data
container similarly to the case when only one experiment is
present. After that different experiments are combined into the same
data container. Result is one TreeSE object with alternative
experiment in altExp slot, or MAE object with multiple experiment in
its experiment slot.

As an example data, we use data from following publication: Hintikka L
_et al._ (2021) Xylo-oligosaccharides in prevention of hepatic
steatosis and adipose tissue inflammation: associating taxonomic and
metabolomic patterns in fecal microbiotas with biclustering
[@Hintikka2021].

In this article, mice were fed with high-fat and low-fat diets with or
without prebiotics.  The purpose of this was to study if prebiotics
would reduce the negative impacts of high-fat diet.

This example data can be loaded from microbiomeDataSets. The data is
already in MAE format. It includes three different experiments:
microbial abundance data, metabolite concentrations, and data about
different biomarkers. Help for importing data into SE object you can
find from
[here](https://microbiome.github.io/OMA/containers.html#loading-experimental-microbiome-data).


```r
# Load the data
data(HintikkaXOData, package = "mia")
mae <- HintikkaXOData
mae
```

```
## A MultiAssayExperiment object of 3 listed
##  experiments with user-defined names and respective classes.
##  Containing an ExperimentList class object of length 3:
##  [1] microbiota: TreeSummarizedExperiment with 12706 rows and 40 columns
##  [2] metabolites: TreeSummarizedExperiment with 38 rows and 40 columns
##  [3] biomarkers: TreeSummarizedExperiment with 39 rows and 40 columns
## Functionality:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DataFrame
##  sampleMap() - the sample coordination DataFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DataFrame
##  assays() - convert ExperimentList to a SimpleList of matrices
##  exportClass() - save data to flat files
```


```r
if(!require(stringr)){
    install.packages("stringr")
    library(stringr)
}
# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]
# Clean taxonomy data, so that names do not include additional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))
# Microbiome data
mae[[1]]
```

```
## class: TreeSummarizedExperiment 
## dim: 12613 40 
## metadata(0):
## assays(1): counts
## rownames(12613): GAYR01026362.62.2014 CVJT01000011.50.2173 ...
##   JRJTB:03787:02429 JRJTB:03787:02478
## rowData names(7): Phylum Class ... Species OTU
## colnames(40): C1 C2 ... C39 C40
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
# Metabolite data
mae[[2]]
```

```
## class: TreeSummarizedExperiment 
## dim: 38 40 
## metadata(0):
## assays(1): nmr
## rownames(38): Butyrate Acetate ... Malonate 1,3-dihydroxyacetone
## rowData names(0):
## colnames(40): C1 C2 ... C39 C40
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
# Biomarker data
mae[[3]]
```

```
## class: TreeSummarizedExperiment 
## dim: 39 40 
## metadata(0):
## assays(1): signals
## rownames(39): Triglycerides_liver CLSs_epi ... NPY_serum Glycogen_liver
## rowData names(0):
## colnames(40): C1 C2 ... C39 C40
## colData names(0):
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: NULL
## rowTree: NULL
## colLinks: NULL
## colTree: NULL
```

## Cross-correlation Analysis

Next we can do the cross-correlation analysis.  Here we analyse if
individual bacteria genera correlates with concentrations of
individual metabolites. This helps as to answer the question: "If this
bacteria is present, is this metabolite's concentration then low or
high"?


```r
# Agglomerate microbiome data at family level
mae[[1]] <- agglomerateByPrevalence(mae[[1]], rank = "Family")
# Does log10 transform for microbiome data
mae[[1]] <- transformCounts(mae[[1]], method = "log10", pseudocount = 1)

# Give unique names so that we do not have problems when we are creating a plot
rownames(mae[[1]]) <- getTaxonomyLabels(mae[[1]])

# Cross correlates data sets
correlations <- testExperimentCrossCorrelation(mae, 
                                               experiment1 = 1,
                                               experiment2 = 2,
                                               assay_name1 = "log10", 
                                               assay_name2 = "nmr",
                                               method = "spearman", 
                                               p_adj_threshold = NULL,
                                               cor_threshold = NULL,
                                               # Remove when mia is fixed
                                               mode = "matrix",
                                               sort = TRUE,
                                               show_warnings = FALSE)
```

Creates the heatmap


```r
if( !require("ComplexHeatmap") ){
    BiocManager::install("ComplexHeatmap")
    library("ComplexHeatmap")
}

# Create a heatmap and store it
plot <- Heatmap(correlations$cor,
                # Print values to cells
                cell_fun = function(j, i, x, y, width, height, fill) {
                    # If the p-value is under threshold
                    if( !is.na(correlations$p_adj[i, j]) & correlations$p_adj[i, j] < 0.05 ){
                        # Print "X"
                        grid.text(sprintf("%s", "X"), x, y, gp = gpar(fontsize = 8, col = "black"))
                        }
                    },
                heatmap_legend_param = list(title = "", legend_height = unit(5, "cm"))
                )
plot
```

![](23_multi-assay_analyses_files/figure-latex/cross-correlation6-1.pdf)<!-- --> 

## Multi-Omics Factor Analysis

Multi-Omics Factor Analysis [@Argelaguet2018] (MOFA) is an
unsupervised method for integrating multi-omic data sets in a
downstream analysis.  It could be seen as a generalization of
principal component analysis. Yet, with the ability to infer a latent
(low-dimensional) representation, shared among the mutliple (-omics)
data sets in hand.

We use the R [MOFA2](https://biofam.github.io/MOFA2/index.html)
package for the analysis, and
[install](https://biofam.github.io/MOFA2/installation.html) the
corresponding dependencies.


```r
if(!require(MOFA2)){
    BiocManager::install("MOFA2")
}

# For inter-operability between Python and R, and setting Python dependencies,
# reticulate package is needed
if(!require(reticulate)){
    install.packages("reticulate")
}

# Let us assume that these have been installed already.
#reticulate::install_miniconda(force = TRUE)
#reticulate::use_miniconda(condaenv = "env1", required = FALSE)
#reticulate::py_install(packages = c("mofapy2"), pip = TRUE, python_version=3.6)
```

The `mae` object could be used straight to create the MOFA model. Yet,
we transform our assays since the model assumes normality per
default. Other distributions that can be used, include Poisson or
Bernoulli.


```r
library(MOFA2)
# For simplicity, classify all high-fat diets as high-fat, and all the low-fat 
# diets as low-fat diets
colData(mae)$Diet <- ifelse(colData(mae)$Diet == "High-fat" | 
                              colData(mae)$Diet == "High-fat + XOS", 
                            "High-fat", "Low-fat")

# Removing duplicates at the microbiome data
# which are also in form e.g. "Ambiguous" and "uncultured" taxa
mae[[1]] <- mae[[1]][!duplicated(rownames(assay(mae[[1]]))), ]

# Transforming microbiome data with rclr
mae[[1]] <- transformCounts(mae[[1]], method = "relabundance")
mae[[1]] <- transformCounts(mae[[1]], assay_name = "relabundance", method = "rclr")

# Transforming metabolomic data with log10
mae[[2]] <- transformCounts(mae[[2]], assay_name = "nmr",
                            MARGIN = "samples",
                            method = "log10")

# Transforming biomarker data with z-transform
mae[[3]] <- transformCounts(mae[[3]], assay_name = "signals",
                            MARGIN = "features",
                            method = "z", pseudocount = 1)

# Removing assays no longer needed
assay(mae[[1]], "counts") <- NULL
assay(mae[[1]], "log10") <- NULL
assay(mae[[2]], "nmr") <- NULL
assay(mae[[3]], "signals") <- NULL

# Building our mofa model
model <- create_mofa_from_MultiAssayExperiment(mae,
                                               groups = "Diet", 
                                               extract_metadata = TRUE)
model
```

```
## Untrained MOFA model with the following characteristics: 
##  Number of views: 3 
##  Views names: microbiota metabolites biomarkers 
##  Number of features (per view): 38 38 39 
##  Number of groups: 2 
##  Groups names: High-fat Low-fat 
##  Number of samples (per group): 20 20 
## 
```

Model options could be defined as follows:


```r
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 5
head(model_opts)
```

```
## $likelihoods
##  microbiota metabolites  biomarkers 
##  "gaussian"  "gaussian"  "gaussian" 
## 
## $num_factors
## [1] 5
## 
## $spikeslab_factors
## [1] FALSE
## 
## $spikeslab_weights
## [1] TRUE
## 
## $ard_factors
## [1] TRUE
## 
## $ard_weights
## [1] TRUE
```

Model's training options are defined with the following:


```r
train_opts <- get_default_training_options(model)
head(train_opts)
```

```
## $maxiter
## [1] 1000
## 
## $convergence_mode
## [1] "fast"
## 
## $drop_factor_threshold
## [1] -1
## 
## $verbose
## [1] FALSE
## 
## $startELBO
## [1] 1
## 
## $freqELBO
## [1] 5
```

Preparing and training the model:


```r
model.prepared <- prepare_mofa(
  object = model,
  model_options = model_opts
)
model.trained <- run_mofa(model.prepared)
```

Visualizing the variance explained:


```r
library(patchwork)
library(ggplot2)
wrap_plots(
    plot_variance_explained(model.trained, x="view", y="factor", plot_total = T),
    nrow = 2
) + plot_annotation(title = "Variance Explained per factor and assay",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
```

![](23_multi-assay_analyses_files/figure-latex/unnamed-chunk-5-1.pdf)<!-- --> 

The top weights for each assay using all 5 factors:


```r
plots <- lapply(c("microbiota", "metabolites","biomarkers"), function(name) {
    plot_top_weights(model.trained,
                     view = name,
                     factors = "all",
                     nfeatures = 10) +
        labs(title = paste0("Top weights of the ", name," assay"))
})
wrap_plots(plots, nrow = 3) & theme(text = element_text(size = 8))
```

![](23_multi-assay_analyses_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 

More tutorials and examples of using the package are found at: [link](https://biofam.github.io/MOFA2/tutorials.html)



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
 [1] ggplot2_3.4.1                  patchwork_1.1.2               
 [3] reticulate_1.28                MOFA2_1.6.0                   
 [5] ComplexHeatmap_2.12.1          stringr_1.5.0                 
 [7] mia_1.7.7                      MultiAssayExperiment_1.24.0   
 [9] TreeSummarizedExperiment_2.1.4 Biostrings_2.66.0             
[11] XVector_0.38.0                 SingleCellExperiment_1.20.0   
[13] SummarizedExperiment_1.28.0    Biobase_2.58.0                
[15] GenomicRanges_1.50.2           GenomeInfoDb_1.34.9           
[17] IRanges_2.32.0                 S4Vectors_0.36.1              
[19] BiocGenerics_0.44.0            MatrixGenerics_1.10.0         
[21] matrixStats_0.63.0-9003        BiocStyle_2.24.0              
[23] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] circlize_0.4.15             corrplot_0.92              
  [3] BiocBaseUtils_1.0.0         plyr_1.8.8                 
  [5] lazyeval_0.2.2              splines_4.2.1              
  [7] BiocParallel_1.32.5         scater_1.26.1              
  [9] digest_0.6.31               foreach_1.5.2              
 [11] yulab.utils_0.0.6           htmltools_0.5.4            
 [13] viridis_0.6.2               fansi_1.0.4                
 [15] magrittr_2.0.3              memoise_2.0.1              
 [17] ScaledMatrix_1.6.0          cluster_2.1.4              
 [19] doParallel_1.0.17           DECIPHER_2.26.0            
 [21] colorspace_2.1-0            rappdirs_0.3.3             
 [23] blob_1.2.3                  ggrepel_0.9.3              
 [25] xfun_0.37                   dplyr_1.1.0                
 [27] crayon_1.5.2                RCurl_1.98-1.10            
 [29] jsonlite_1.8.4              graph_1.74.0               
 [31] iterators_1.0.14            ape_5.7                    
 [33] glue_1.6.2                  gtable_0.3.1               
 [35] zlibbioc_1.44.0             GetoptLong_1.0.5           
 [37] DelayedArray_0.24.0         BiocSingular_1.14.0        
 [39] Rhdf5lib_1.18.2             shape_1.4.6                
 [41] HDF5Array_1.24.2            scales_1.2.1               
 [43] pheatmap_1.0.12             DBI_1.1.3                  
 [45] Rcpp_1.0.10                 viridisLite_0.4.1          
 [47] decontam_1.18.0             clue_0.3-64                
 [49] tidytree_0.4.2              bit_4.0.5                  
 [51] rsvd_1.0.5                  dir.expiry_1.4.0           
 [53] RColorBrewer_1.1-3          farver_2.1.1               
 [55] pkgconfig_2.0.3             XML_3.99-0.13              
 [57] scuttle_1.8.4               uwot_0.1.14                
 [59] CodeDepends_0.6.5           here_1.0.1                 
 [61] utf8_1.2.3                  labeling_0.4.2             
 [63] tidyselect_1.2.0            rlang_1.0.6                
 [65] reshape2_1.4.4              munsell_0.5.0              
 [67] tools_4.2.1                 cachem_1.0.6               
 [69] cli_3.6.0                   DirichletMultinomial_1.40.0
 [71] generics_0.1.3              RSQLite_2.3.0              
 [73] evaluate_0.20               fastmap_1.1.0              
 [75] yaml_2.3.7                  knitr_1.42                 
 [77] bit64_4.0.5                 purrr_1.0.1                
 [79] nlme_3.1-162                sparseMatrixStats_1.10.0   
 [81] compiler_4.2.1              beeswarm_0.4.0             
 [83] filelock_1.0.2              png_0.1-8                  
 [85] treeio_1.22.0               tibble_3.1.8               
 [87] stringi_1.7.12              highr_0.10                 
 [89] basilisk.utils_1.8.0        forcats_1.0.0              
 [91] lattice_0.20-45             Matrix_1.5-3               
 [93] vegan_2.6-4                 permute_0.9-7              
 [95] vctrs_0.5.2                 pillar_1.8.1               
 [97] lifecycle_1.0.3             rhdf5filters_1.8.0         
 [99] BiocManager_1.30.19         GlobalOptions_0.1.2        
[101] BiocNeighbors_1.16.0        cowplot_1.1.1              
[103] bitops_1.0-7                irlba_2.3.5.1              
[105] R6_2.5.1                    bookdown_0.32              
[107] gridExtra_2.3               vipor_0.4.5                
[109] codetools_0.2-19            MASS_7.3-58.2              
[111] rhdf5_2.40.0                rprojroot_2.0.3            
[113] rjson_0.2.21                withr_2.5.0                
[115] GenomeInfoDbData_1.2.9      mgcv_1.8-41                
[117] parallel_4.2.1              beachmat_2.14.0            
[119] tidyr_1.3.0                 basilisk_1.8.1             
[121] rmarkdown_2.20              DelayedMatrixStats_1.20.0  
[123] Rtsne_0.16                  Cairo_1.6-0                
[125] ggbeeswarm_0.7.1           
```
</div>

