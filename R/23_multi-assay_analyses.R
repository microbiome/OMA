## ----setup, echo=FALSE, results="asis"-----------------
library(rebook)
chapterPreamble()


## ----load-pkg-data-------------------------------------
library(mia)


## ----cross-correlation1--------------------------------
# Load the data
data(HintikkaXOData, package = "mia")
mae <- HintikkaXOData


## ----cross-correlation2--------------------------------
library(stringr)
# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]
# Clean taxonomy data, so that names do not include additional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))


## ------------------------------------------------------
# Available alternative experiments
experiments(mae)


## ------------------------------------------------------
# Microbiome data
getWithColData(mae, "microbiota")


## ----cross-correlation3--------------------------------
# Metabolite data
getWithColData(mae, "metabolites")


## ----cross-correlation4--------------------------------
# Biomarker data
getWithColData(mae, "biomarkers")


## ----cross-correlation5--------------------------------
# Agglomerate microbiome data at family level
mae[[1]] <- mergeFeaturesByPrevalence(mae[[1]], rank = "Family")
# Does log10 transform for microbiome data
mae[[1]] <- transformAssay(mae[[1]], method = "log10", pseudocount = TRUE)

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


## ----cross-correlation6, fig.width=10, fig.height=8----
library(ComplexHeatmap) 

# Create a heatmap and store it
plot <- Heatmap(correlations$cor,
                # Print values to cells
                cell_fun = function(j, i, x, y, width, height, fill) {
                    # If the p-value is under threshold
                    if( !is.na(correlations$p_adj[i, j]) & correlations$p_adj[i, j] < 0.05 ){
                        # Print "X"
                        grid.text(sprintf("%s", "X"), x, y, gp = gpar(fontsize = 10, col = "#1dff00"))
                        }
                    },
                heatmap_legend_param = list(title = "", legend_height = unit(5, "cm"))
                )
plot


## ----MOFA2, message=FALSE, warning=FALSE---------------
library(MOFA2)

# For inter-operability between Python and R, and setting Python dependencies,
# reticulate package is needed
library(reticulate)
# Let us assume that these have been installed already.
#reticulate::install_miniconda(force = TRUE)
#reticulate::use_miniconda(condaenv = "env1", required = FALSE)
#reticulate::py_install(packages = c("mofapy2"), pip = TRUE, python_version=3.6)


## ---- message=FALSE, warning=FALSE---------------------
library(MOFA2)
# For simplicity, classify all high-fat diets as high-fat, and all the low-fat 
# diets as low-fat diets
colData(mae)$Diet <- ifelse(colData(mae)$Diet == "High-fat" | 
                              colData(mae)$Diet == "High-fat + XOS", 
                            "High-fat", "Low-fat")

# Transforming microbiome data with rclr
mae[[1]] <- transformAssay(mae[[1]], method = "relabundance")
mae[[1]] <- transformAssay(mae[[1]], assay.type = "relabundance", method = "rclr")

# Transforming metabolomic data with log10
mae[[2]] <- transformAssay(mae[[2]], assay.type = "nmr",
                            MARGIN = "samples",
                            method = "log10")

# Transforming biomarker data with z-transform
mae[[3]] <- transformAssay(mae[[3]], assay.type = "signals",
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


## ---- message=FALSE, warning=FALSE---------------------
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 5
head(model_opts)


## ---- message=FALSE, warning=FALSE---------------------
train_opts <- get_default_training_options(model)
head(train_opts)


## ---- message=FALSE, warning=FALSE---------------------
model.prepared <- prepare_mofa(
  object = model,
  model_options = model_opts
)

# Some systems may require the specification `use_basilisk = TRUE`
# so it has been added to the following code
model.trained <- run_mofa(model.prepared, use_basilisk = TRUE)


## ---- message=FALSE, warning=FALSE, fig.height=8, fig.width=10----
library(patchwork)
library(ggplot2)

plot_list <- plot_variance_explained(model.trained,
                                     x = "view", y = "factor",
                                     plot_total = T)

wrap_plots(plot_list, nrow = 2) +
  plot_annotation(title = "Variance Explained per factor and assay",
                  theme = theme(plot.title = element_text(hjust = 0.5)))


## ---- warning=FALSE, message=FALSE, fig.height=10, fig.width=10----
custom_plotter <- function(name) {
  
  p <- plot_top_weights(model.trained,
                        view = name,
                        factors = "all",
                        nfeatures = 10) +
    labs(title = paste0("Top weights of the ", name, " assay"))
  
}

plot_list <- lapply(c("microbiota", "metabolites", "biomarkers"), custom_plotter)

wrap_plots(plot_list, nrow = 3) & theme(text = element_text(size = 8))

