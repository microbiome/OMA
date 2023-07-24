## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----load-pkg-data1-----------------------------------------------------------
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns


## ----hclust1------------------------------------------------------------------
library(mia)
library(vegan)

# Load experimental data
data(peerj13075)
(tse <- peerj13075)


## ----hclust2------------------------------------------------------------------
if( !require(NbClust) ){install.packages("NbClust"); library(NbClust)}
if( !require(cobiclust) ){install.packages("cobiclust"); library(cobiclust)}

# Apply transformation
tse <- transformAssay(tse, method = "relabundance")
# Get the assay
assay <- assay(tse, "relabundance")
# Transpose assay --> samples are now in rows --> we are clustering samples
assay <- t(assay)

# Calculate distances
diss <- vegdist(assay, method = "bray")

# Perform hierarchical clustering
hc <- hclust(diss, method = "complete")

# To visualize, convert hclust object into dendrogram object
dendro <- as.dendrogram(hc)

# Plot dendrogram
plot(dendro)


## ----hclust3------------------------------------------------------------------
# Determine the optimal number of clusters
res <- NbClust(diss = diss, distance = NULL, method = "ward.D2",
               index = "silhouette")

res$Best.nc


## ----hclust4------------------------------------------------------------------
if( !require(dendextend) ){
    install.packages("dendextend")
    library(dendextend)
}

# Find clusters
cutree(hc, k = 15) 

# Making colors for 6 clusters
col_val_map <- randomcoloR::distinctColorPalette(15) %>%
     as.list() %>% setNames(paste0("clust_",seq(15)))

dend <- color_branches(dendro, k=15, col=unlist(col_val_map))
labels(dend) <- NULL
plot(dend)



## ----kmeans1------------------------------------------------------------------
if( !require(factoextra) ){
    install.packages("factoextra")
    library(factoextra)
}

# Convert dist object into matrix
diss <- as.matrix(diss)
# Perform silhouette analysis and plot the result
fviz_nbclust(diss, kmeans, method = "silhouette")


## ----kmeans2------------------------------------------------------------------
library(scater)

# The first step is random, add seed for reproducibility
set.seed(15463)
# Perform k-means clustering with 3 clusters
km <- kmeans(diss, 3, nstart = 25)
# Add the result to colData
colData(tse)$clusters <- as.factor(km$cluster)

# Perform PCoA so that we can visualize clusters
tse <- runMDS(tse, assay.type = "relabundance", FUN = vegan::vegdist, method = "bray")

# Plot PCoA and color clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")


## ----dmm1---------------------------------------------------------------------
# Runs model and calculates the most likely number of clusters from 1 to 7.
# Since this is a large dataset it takes long computational time.
# For this reason we use only a subset of the data; agglomerated by Phylum as a rank.
tse <- GlobalPatterns
tse <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree=TRUE)


## ----dmm2---------------------------------------------------------------------
tse_dmn <- mia::runDMN(tse, name = "DMN", k = 1:7)


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
dmn_group <- calculateDMNgroup(tse_dmn, variable = "SampleType",  assay.type = "counts",
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



## -----------------------------------------------------------------------------
# Does clr transformation. Pseudocount is added, because data contains zeros.
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformAssay(tse, assay.type = "pseudo", method = "relabundance")
tse <- transformAssay(tse, "relabundance", method = "clr")

library(scater)

# Does principal coordinate analysis
df <- calculateMDS(tse, assay.type = "clr", method = "euclidean")

# Creates a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = df[,1], 
                                pcoa2 = df[,2])



## -----------------------------------------------------------------------------
# Creates a data frame that contains principal coordinates and DMM information
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = vec)
# Creates a plot
euclidean_dmm_plot <- ggplot(data = euclidean_dmm_pcoa_df, 
                             aes(x=pcoa1, y=pcoa2,
                                 color = dmm_component)) +
  geom_point() +
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
tse <- transformAssay(tse, method = "rclr")

# Performing and storing UMAP
tse <- runUMAP(tse, name="UMAP", assay.type="rclr")

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


## ----load-pkg-data2-----------------------------------------------------------
library(mia)
data("HintikkaXOData")
mae <- HintikkaXOData


## ----cobiclust_1--------------------------------------------------------------
# Subset data in the first experiment
mae[[1]] <- subsetByPrevalentTaxa(mae[[1]], rank = "Genus", prevalence = 0.2, detection = 0.001)
# clr-transform in the first experiment
mae[[1]] <- transformAssay(mae[[1]], method = "relabundance")
mae[[1]] <- transformAssay(mae[[1]], "relabundance", method = "rclr")


## ----cobiclust_2--------------------------------------------------------------
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


## ----cobiclust_3, fig.width=14, fig.height=12---------------------------------

if(!require(pheatmap)){
    install.packages("pheatmap")
    library(pheatmap)
}
# z-transform for heatmap
mae[[1]] <- transformAssay(mae[[1]], assay.type = "rclr",
                            MARGIN = "features",
                            method = "z", name = "clr_z")

# Create annotations. When column names are equal, they should share levels. 
# Here samples include 3 clusters, and taxa 2. That is why we have to make 
# column names unique. 
annotation_col <- data.frame(colData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_col) <- "col_clusters"

annotation_row <- data.frame(rowData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_row) <- "row_clusters"


## ----cobiclust_3b, fig.width=14, fig.height=12--------------------------------
pheatmap(assay(mae[[1]], "clr_z"), cluster_rows = F, cluster_cols = F, 
         annotation_col = annotation_col,
         annotation_row = annotation_row)


## ----cobiclust_4--------------------------------------------------------------
if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
}
if(!require(patchwork)){
    install.packages("patchwork")
    library(patchwork)
}

# ggplot requires data in melted format
melt_assay <- meltAssay(mae[[1]], assay.type = "rclr", add_col_data = T, add_row_data = T)

# patchwork two plots side-by-side
p1 <- ggplot(melt_assay) +
  geom_boxplot(aes(x = clusters.x, y = rclr)) +
  labs(x = "Taxa clusters")

p2 <- ggplot(melt_assay) +
  geom_boxplot(aes(x = clusters.y, y = rclr)) +
  labs(x = "Sample clusters")

p1 + p2


## ----biclust_1----------------------------------------------------------------
# Samples must be in equal order 
# (Only 1st experiment  was ordered in cobiclust step leading to unequal order)
mae[[1]] <- mae[[1]][ , colnames(mae[[2]]) ]

# Make rownames unique since it is require by other steps
rownames(mae[[1]]) <- make.unique(rownames(mae[[1]]))
# Calculate correlations
corr <- getExperimentCrossCorrelation(mae, 1, 2, 
                                      assay.type1 = "rclr", 
                                      assay.type2 = "nmr", 
                                      mode = "matrix", 
                                      cor_threshold = 0.2)



## ----biclust_2----------------------------------------------------------------
# Set seed for reproducibility
set.seed(3973)

# Find biclusters
library(biclust)
bc <- biclust(corr, method=BCPlaid(), fit.model = y ~ m,
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10,
              iter.startup = 10, iter.layer = 100, verbose = FALSE)

bc


## ----biclust_3----------------------------------------------------------------
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

# Input clusters, and how many observations there should be, i.e.,
# the number of samples or features
.manipulate_bc_data <- function(bc_clusters, assay, row_col){
  # Get right dimension
  dim <- ifelse(row_col == "col", ncol(assay), nrow(assay))
  # Get column/row names
  if( row_col == "col" ){
    names <- colnames(assay)
  } else{
    names <- rownames(assay)
  }
  
  # If no clusters were found, create one. Otherwise create additional
  # cluster which
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


## ----biclust_4----------------------------------------------------------------
# Get biclusters
bcs <- .get_biclusters_from_biclust(bc, corr)

bicluster_rows <- bcs$bc_rows
bicluster_columns <- bcs$bc_columns

# Print biclusters for rows
head(bicluster_rows)


## ----biclust_5----------------------------------------------------------------
# Function for obtaining sample-wise sum, mean, median, and mean variance
# for each cluster

.sum_mean_median_var <- function(tse1, tse2, assay.type1, assay.type2, clusters1, clusters2){
  
  list <- list()
  # Create a data frame that includes all the information
  for(i in 1:ncol(clusters1) ){
    # Subset data based on cluster
    tse_subset1 <- tse1[clusters1[,i], ]
    tse_subset2 <- tse2[clusters2[,i], ]
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
    
    list[[i]] <- data.frame(sample = colnames(tse1), sum1, sum2, mean1, mean2, 
                     median1, median2, var1, var2)
  }

  return(list)
}

# Calculate info
df <- .sum_mean_median_var(mae[[1]], mae[[2]], "rclr", "nmr", bicluster_rows, bicluster_columns)


## ----biclust_6, fig.width=14, fig.height=6, fig.show="keep", out.width="33%"----
pics <- list()
for(i in seq_along(df)){
  pics[[i]] <- ggplot(df[[i]])  +
      geom_point(aes(x = median1, y = median2)) + 
      labs(title = paste0("Cluster ", i),
           x = "Taxa (rclr median)",
           y = "Metabolites (abs. median)")
  print(pics[[i]])
}
# pics[[1]] + pics[[2]] + pics[[3]]


## ----biclust_7----------------------------------------------------------------
bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))


## ----biclust_8, fig.width=10, fig.height=10-----------------------------------
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


## ----biclust_9----------------------------------------------------------------
# Calculate cross-correlation
corr <- getExperimentCrossCorrelation(mae, 1, 1, 
                                      assay.type1 = "rclr", assay.type2 = "rclr", 
                                      mode = "matrix",
                                      cor_threshold = 0.2, verbose = F, show_warning = F)

# Find biclusters
library(biclust)
bc <- biclust(corr, method=BCPlaid(), fit.model = y ~ m,
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10,
              iter.startup = 10, iter.layer = 100, verbose = FALSE)


## ----biclust_10---------------------------------------------------------------
# Get biclusters
bcs <- .get_biclusters_from_biclust(bc, corr)

bicluster_rows <- bcs$bc_rows
bicluster_columns <- bcs$bc_columns


## ----biclust_11---------------------------------------------------------------
# Create a column that combines information
# If row/column includes in multiple clusters, cluster numbers are separated with "_&_"
bicluster_columns$clusters <- apply(bicluster_columns, 1, 
                                    function(x){paste(paste(which(x)), collapse = "_&_") })
bicluster_columns <- bicluster_columns[, "clusters", drop = FALSE]

bicluster_rows$clusters <- apply(bicluster_rows, 1, 
                                 function(x){paste(paste(which(x)), collapse = "_&_") })
bicluster_rows <- bicluster_rows[, "clusters", drop = FALSE]


## ----biclust_12, fig.width=14, fig.height=12----------------------------------
# Convert boolean values into factor
bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))

pheatmap(corr, cluster_cols = F, cluster_rows = F,
         annotation_col = bicluster_columns, 
         annotation_row = bicluster_rows)

