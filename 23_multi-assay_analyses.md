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

Multi-omics approaches integrate data from multiple sources. For
example, we can integrate taxonomic abundance profiles with
metabolomic or other biomolecular profiling data to observe
associations, make predictions, or aim at causal
inferences. Integrating evidence across multiple sources can lead to
enhanced predictions, more holistic understanding, or facilitate the
discovery of novel biomarkers. In this section we demonstrate common
multi-assay data integration tasks.

Cross-correlation analysis is a straightforward approach that can
reveal assocation strengths and types between data sets. For instance,
we can analyze if higher presence of a specific taxon equals to higher
levels of a biomolecule.

The analyse can be facilitated by the multi-assay data containers,
_TreeSummarizedExperiment_ and _MultiAssayExperiment_. These are
scalable and contain different types of data in a single container,
making this framework particularly suited for multi-assay microbiome
data incorporating different types of complementary data sources in a
single reproducible workflow. Th different solutions for varying data
integration needs are discussed in more detail in Section
\@ref(containers). Another experiment can be stored in _altExp_ slot
of SE data container or both experiments can be stored side-by-side in
MAE data container (see the sections \@ref(alt-exp) and \@ref(mae) to
learn more about altExp and MAE objects, respectively). Different
experiments are first imported into these data containers similarly to
the case when only one experiment is present. After that, the
different experiments can be combined into the same multi-assay data
container. The result is one TreeSE object with alternative experiment in
altExp slot, or MAE object with multiple experiment in its experiment
slot, for instance.

As an example data, we use data from the following publication: Hintikka L
_et al._ (2021) Xylo-oligosaccharides in prevention of hepatic
steatosis and adipose tissue inflammation: associating taxonomic and
metabolomic patterns in fecal microbiota with biclustering

[@Hintikka2021].

In this study, mice were fed with high-fat and low-fat diets with or
without prebiotics.  The purpose of this was to study whether prebiotics
reduce negative impacts of a high-fat diet.

This example data can be loaded from microbiomeDataSets. The data is
already in MAE format. It includes three different experiments:
microbial abundance data, metabolite concentrations, and data about
different biomarkers.
=======
[@Hintikka2021]. In this article, mice were fed with high-fat and
low-fat diets with or without prebiotics.  The purpose of this was to
study if prebiotics would reduce the negative impact of high-fat diet.

This example data is readily available in the MultiAssayExperiment
format. It includes three different experiments: microbial abundance
data, metabolite concentrations, and data about different
biomarkers. If you like to construct the same data object from the
original files instead, [Here](https://microbiome.github.io/OMA/containers.html#loading-experimental-microbiome-data) 
you can find help for importing data into a SE object.


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
library(stringr)
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


Next we can do the cross-correlation analysis. Let us analyse if
individual bacteria genera correlates with concentrations of
individual metabolites. This helps to answer the question: "If this
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
                                               assay.type1 = "log10", 
                                               assay.type2 = "nmr",
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
library("ComplexHeatmap") 

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
(low-dimensional) representation, shared among the multiple (-omics)
data sets in hand.

We use the R [MOFA2](https://biofam.github.io/MOFA2/index.html)
package for the analysis, and
[install](https://biofam.github.io/MOFA2/installation.html) the
corresponding dependencies.


```r
library(MOFA2)

# For inter-operability between Python and R, and setting Python dependencies,
# reticulate package is needed
library(reticulate)
# Let us assume that these have been installed already.
#reticulate::install_miniconda(force = TRUE)
#reticulate::use_miniconda(condaenv = "env1", required = FALSE)
#reticulate::py_install(packages = c("mofapy2"), pip = TRUE, python_version=3.6)
```

The `mae` object could be used straight to create the MOFA model. Yet,
we transform our assays since the model assumes normality per
default. We can also use Poisson or Bernoulli distributions among others.


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
mae[[1]] <- transformCounts(mae[[1]], assay.type = "relabundance", method = "rclr")

# Transforming metabolomic data with log10
mae[[2]] <- transformCounts(mae[[2]], assay.type = "nmr",
                            MARGIN = "samples",
                            method = "log10")

# Transforming biomarker data with z-transform
mae[[3]] <- transformCounts(mae[[3]], assay.type = "signals",
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
## [1] FALSE
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

# Some systems may require the specification `use_basilisk = TRUE`
# so it has been added to the following code
model.trained <- run_mofa(model.prepared, use_basilisk = TRUE)
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

The top weights for each assay using all five factors:


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

More tutorials and examples of using the package are found at [link](https://biofam.github.io/MOFA2/tutorials.html)




