## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----load-pkg-data------------------------------------------------------------
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns


## ----dmm----------------------------------------------------------------------
# Runs model and calculates the most likely number of clusters from 1 to 7.
# Since this is a large dataset it takes long computational time.
# For this reason we use only a subset of the data; agglomerated by Phylum as a rank.
tse <- GlobalPatterns
tse <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree=TRUE)
tse_dmn <- runDMN(tse, name = "DMN", k = 1:7)


## -----------------------------------------------------------------------------
# It is stored in metadata
tse_dmn


## -----------------------------------------------------------------------------
names(metadata(tse_dmn))


## -----------------------------------------------------------------------------
getDMN(tse_dmn)


## -----------------------------------------------------------------------------
library(miaViz)
plotDMNFit(tse_dmn, type = "laplace")


## -----------------------------------------------------------------------------
getBestDMNFit(tse_dmn, type = "laplace")


## -----------------------------------------------------------------------------
dmn_group <- calculateDMNgroup(tse_dmn, variable = "SampleType",  exprs_values = "counts",
                               k = 2, seed=.Machine$integer.max)

dmn_group


## -----------------------------------------------------------------------------
DirichletMultinomial::mixturewt(getBestDMNFit(tse_dmn))


## -----------------------------------------------------------------------------
head(DirichletMultinomial::mixture(getBestDMNFit(tse_dmn)))


## -----------------------------------------------------------------------------
head(DirichletMultinomial::fitted(getBestDMNFit(tse_dmn)))


## -----------------------------------------------------------------------------
prob <- DirichletMultinomial::mixture(getBestDMNFit(tse_dmn))
# Add column names
colnames(prob) <- c("comp1", "comp2")

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec <- colnames(prob)[max.col(prob,ties.method = "first")]

# Add info to colData
colData(tse)$dmm_component <- vec 


## -----------------------------------------------------------------------------
# Does clr transformation. Pseudocount is added, because data contains zeros.
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformCounts(tse, assay_name = "pseudo", method = "relabundance")
tse <- transformCounts(tse, assay_name = "relabundance", method = "clr")

library(scater)
# Calculate PCoA
tse <- runMDS(tse, exprs_values = "clr", method = "euclidean")



## -----------------------------------------------------------------------------
# Create a plot
euclidean_dmm_plot <- plotReducedDim(tse, "MDS", colour_by = "dmm_component") +
    labs(x = "Coordinate 1",
       y = "Coordinate 2",
       title = "PCoA with Aitchison distances") +  
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_dmm_plot


## -----------------------------------------------------------------------------
if(!require(bluster)){
  BiocManager::install("bluster")
}


## ---- message=FALSE, warning=FALSE--------------------------------------------
library(bluster)
library(patchwork) # For arranging several plots as a grid
library(scater)

data("enterotype", package="mia")
tse <- enterotype
tse <- transformCounts(tse, method = "relabundance")
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


## ---- message=FALSE, warning=FALSE--------------------------------------------

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

