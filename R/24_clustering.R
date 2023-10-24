## ----setup, echo=FALSE, results="asis"-----------------
library(rebook)
chapterPreamble()


## ----load-pkg-data1------------------------------------
library(mia)
data("enterotype", package = "mia")
tse <- enterotype


## ----load_bluster--------------------------------------
# Load dependencies
library(bluster)
library(kableExtra)

# Apply transformation
tse <- transformAssay(tse, method = "relabundance")


## ----bluster_sample_hclust1----------------------------
# Simple use of the hierarchical clustering. Here, the default parameters
# set the cut height to half of the dendrogram height.
tse <- cluster(tse, assay.type = "relabundance", 
               MARGIN = "samples", HclustParam())

# Check the result contained in the clusters part of colData

head(colData(tse)$clusters)


## ----bluster_sample_plot-------------------------------
library(scater)

# Add the MDS dimensions for plotting
tse <- runMDS(tse, assay.type = "relabundance", 
              FUN = vegan::vegdist, method = "bray")

# Plot the clusters
plotReducedDim(tse, "MDS", colour_by = "clusters")


## ----hclust1-------------------------------------------
library(vegan)

# Load experimental data
tse <- enterotype


## ----hclust2-------------------------------------------
# Apply transformation
tse <- transformAssay(tse, method = "relabundance")

# Do the clustering
tse <- cluster(tse,
               assay.type = "relabundance",
               MARGIN = "samples",
               HclustParam(method = "complete",
                           dist.fun = vegdist,
                           metric = "bray"),
               full = TRUE,
               clust.col = "Hclust")


## ----hclust3-------------------------------------------
library(dendextend)

# Get hclust data from metadata
hclust_data <- metadata(tse)$clusters$hclust

# Get the dendrogram object
dendro <- as.dendrogram(hclust_data)

# Plot dendrogram

dendro %>% set("labels", NULL) %>% plot()


## ----hclust4-------------------------------------------
# Get the clusters
head(colData(tse)$Hclust)


## ----hclust5-------------------------------------------
library(NbClust)
diss <- metadata(tse)$clusters$dist

# Apply the silhouette analysis on the distance matrix
res <- NbClust(diss = diss, distance = NULL, method = "ward.D2",
               index = "silhouette")

res$Best.nc


## ----hclust6-------------------------------------------
library(dendextend)

# Get optimal number of clusters
k <- res$Best.nc[1]

# Making colors for 2 clusters
col_val_map <- randomcoloR::distinctColorPalette(k) %>%
    as.list() %>% 
    setNames(paste0("clust_", seq(k)))

dend <- color_branches(dendro, k = k, col = unlist(col_val_map))
labels(dend) <- NULL
plot(dend)


## ----dmm1----------------------------------------------
# Get the data
data("GlobalPatterns", package = "mia")
tse <- GlobalPatterns

# Agglomerate by rank
tse <- mergeFeaturesByRank(tse, rank = "Phylum", agglomerateTree = TRUE)


## ----dmm2----------------------------------------------
# Run the model and calculates the most likely number of clusters from 1 to 7
tse <- cluster(tse, name = "DMM", DmmParam(k = 1:7, type = "laplace"), 
                   MARGIN = "samples", full = TRUE)


## ----dmm3----------------------------------------------
# The dmm info is stored in the metadata under the 'DMM' column
tse


## ----dmm4----------------------------------------------
head(metadata(tse)$DMM$dmm,3)


## ----dmm5----------------------------------------------
BiocManager::install("microbiome/miaViz")
library(miaViz)
plotDMNFit(tse, type = "laplace", name = "DMM")


## ----dmm6----------------------------------------------
# Get the model that has the best fit
bestFit <- metadata(tse)$DMM$dmm[[metadata(tse)$DMM$best]]
bestFit


## ----dmm_group-----------------------------------------
dmm_group <- calculateDMNgroup(tse, variable = "SampleType", 
                               assay.type = "counts", k = 2, 
                               seed = .Machine$integer.max)

dmm_group


## ----dmm7----------------------------------------------
DirichletMultinomial::mixturewt(bestFit)


## ----dmm8----------------------------------------------
prob <- metadata(tse)$DMM$prob
head(prob)


## ----dmm9----------------------------------------------
head(DirichletMultinomial::fitted(bestFit))


## ----dmm10---------------------------------------------
# Do clr transformation. Pseudocount is added, because data contains zeros.
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformAssay(tse, assay.type = "pseudo", method = "relabundance")
tse <- transformAssay(tse, "relabundance", method = "clr")

# Do principal coordinate analysis
df <- calculateMDS(tse, assay.type = "clr", method = "euclidean")

# Create a data frame from principal coordinates
euclidean_pcoa_df <- data.frame(pcoa1 = df[, 1], pcoa2 = df[, 2])


## ----dmm11---------------------------------------------
# Create a data frame that contains principal coordinates and DMM information
euclidean_dmm_pcoa_df <- cbind(euclidean_pcoa_df,
                               dmm_component = colData(tse)$clusters)

# Create a plot
euclidean_dmm_plot <- ggplot(data = euclidean_dmm_pcoa_df,
                             aes(x = pcoa1, y = pcoa2, color = dmm_component)) +
    geom_point() +
    labs(x = "Coordinate 1",y = "Coordinate 2", 
         title = "PCoA with Aitchison distances") +
    theme(plot.title = element_text(size = 12, # makes titles smaller
                                    hjust = 0.5)) 

euclidean_dmm_plot


## ----load-pkg-data2------------------------------------
library(cobiclust)
data("HintikkaXOData")
mae <- HintikkaXOData


## ----cobiclust_1---------------------------------------
# Subset data in the first experiment
mae[[1]] <- subsetByPrevalentFeatures(mae[[1]], rank = "Genus", 
                                      prevalence = 0.2, 
                                      detection = 0.001)

# rclr-transform in the first experiment
mae[[1]] <- transformAssay(mae[[1]], method = "rclr")


## ----cobiclust_2---------------------------------------
# Do clustering using counts table
clusters <- cobiclust(assay(mae[[1]], "counts"))

# Get clusters
row_clusters <- clusters$classification$rowclass
col_clusters <- clusters$classification$colclass

# Add clusters to rowdata and coldata
rowData(mae[[1]])$clusters <- factor(row_clusters)
colData(mae[[1]])$clusters <- factor(col_clusters)

# Order data based on clusters
mae[[1]] <- mae[[1]][order(rowData(mae[[1]])$clusters),
                     order(colData(mae[[1]])$clusters)]

# Print clusters
clusters$classification


## ----cobiclust_3a, fig.width=14, fig.height=12---------
library(pheatmap)
# z-transform for heatmap
mae[[1]] <- transformAssay(mae[[1]], assay.type = "rclr",
                            MARGIN = "features", method = "z", name = "rclr_z")

# Create annotations. When column names are equal, they should share levels.
# Here samples include 3 clusters, and taxa 2. That is why we have to make
# column names unique.
annotation_col <- data.frame(colData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_col) <- "col_clusters"

annotation_row <- data.frame(rowData(mae[[1]])[, "clusters", drop = F])
colnames(annotation_row) <- "row_clusters"


## ----cobiclust_3b, fig.width=14, fig.height=12---------
pheatmap(assay(mae[[1]], "rclr_z"), cluster_rows = F, cluster_cols = F,
         annotation_col = annotation_col, annotation_row = annotation_row)


## ----cobiclust_4---------------------------------------
library(ggplot2)
library(patchwork)

# ggplot requires data in melted format
melt_assay <- meltAssay(mae[[1]], assay.type = "rclr", 
                        add_col_data = T, add_row_data = T)

# patchwork two plots side-by-side
p1 <- ggplot(melt_assay) +
    geom_boxplot(aes(x = clusters.x, y = rclr)) +
    labs(x = "Taxa clusters")

p2 <- ggplot(melt_assay) +
    geom_boxplot(aes(x = clusters.y, y = rclr)) +
    labs(x = "Sample clusters")

p1 + p2


## ----biclust_1-----------------------------------------
# Samples must be in equal order
# (Only 1st experiment was ordered in cobiclust step leading to unequal order)
mae[[1]] <- mae[[1]][, colnames(mae[[2]])]

# Make rownames unique since it is required by other steps
rownames(mae[[1]]) <- make.unique(rownames(mae[[1]]))

# Transform the metabolites to be in log basis
mae[[2]] <- transformAssay(mae[[2]], assay.type = "nmr", method = "log10")

# Add missing data to the metabolites
replace_na <- function(row) {
    na_indices <- which(is.na(row))
    non_na_values <- row[!is.na(row)]
    row[na_indices] <- sample(non_na_values, length(na_indices), replace = TRUE)
    row
}
assay(mae[[2]], "log10") <- t(apply(assay(mae[[2]], "log10"), 1, replace_na))


## ----biclust_2-----------------------------------------
# Calculate correlations
corr <- getExperimentCrossCorrelation(mae, 1, 2, assay.type1 = "rclr",
                                      assay.type2 = "log10", mode = "matrix",
                                      correlation = "spearman")


## ----biclust_3-----------------------------------------
library(biclust)
# Set seed for reproducibility
set.seed(3973)

# Find biclusters
bc <- biclust(corr, method = BCPlaid(), verbose = FALSE)

bc


## ----biclust_4-----------------------------------------
# Functions for obtaining biclust information

# Get clusters for rows and columns
.get_biclusters_from_biclust <- function(bc, assay) {
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
.manipulate_bc_data <- function(bc_clusters, assay, row_col) {
    # Get right dimension
    dim <- ifelse(row_col == "col", ncol(assay), nrow(assay))
    # Get column/row names
    if (row_col == "col") {
        names <- colnames(assay)
    } else {
        names <- rownames(assay)
    }

    # If no clusters were found, create one. Otherwise create additional
    # cluster which
    # contain those samples that are not included in clusters that were found.
    if (nrow(bc_clusters) != dim) {
        bc_clusters <- data.frame(cluster = rep(TRUE, dim))
    } else {
        # Create additional cluster that includes those samples/features that
        # are not included in other clusters.
        vec <- ifelse(rowSums(bc_clusters) > 0, FALSE, TRUE)

        # If additional cluster contains samples, then add it
        if (any(vec)) {
            bc_clusters <- cbind(bc_clusters, vec)
        }
    }
    
    # Adjust row and column names
    rownames(bc_clusters) <- names
    colnames(bc_clusters) <- paste0("cluster_", 1:ncol(bc_clusters))
    return(bc_clusters)
}


## ----biclust_5-----------------------------------------
# Get biclusters
bcs <- .get_biclusters_from_biclust(bc, corr)

bicluster_rows <- bcs$bc_rows
bicluster_columns <- bcs$bc_columns

# Print biclusters for rows
head(bicluster_rows)


## ----biclust_6-----------------------------------------
# Function for obtaining sample-wise sum, mean, median, and mean variance
# for each cluster

.sum_mean_median_var <- function(tse1, tse2, assay.type1, assay.type2, clusters1, clusters2) {
    list <- list()
    # Create a data frame that includes all the information
    for (i in 1:ncol(clusters1)) {
        # Subset data based on cluster
        tse_subset1 <- tse1[clusters1[, i], ]
        tse_subset2 <- tse2[clusters2[, i], ]
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
        
        list[[i]] <- data.frame(sample = colnames(tse1), sum1, sum2, mean1, 
                                 mean2, median1, median2, var1, var2)
    }
    return(list)
}

# Calculate info
df <- .sum_mean_median_var(mae[[1]], mae[[2]], "rclr", "log10", bicluster_rows, bicluster_columns)


## ----biclust_7, fig.width=14, fig.height=6, fig.show="keep", out.width="33%"----
pics <- list()
for (i in seq_along(df)) {
    pics[[i]] <- ggplot(df[[i]]) +
        geom_point(aes(x = median1, y = median2)) +
        labs(title = paste0("Cluster ", i), x = "Taxa (rclr median)",
             y = "Metabolites (abs. median)")
    print(pics[[i]])
}
pics[[1]] + pics[[2]] + pics[[3]]


## ----biclust_8-----------------------------------------
bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))


## ----biclust_9, fig.width=10, fig.height=10------------
# Adjust colors for all clusters
if (ncol(bicluster_rows) > ncol(bicluster_columns)) {
    cluster_names <- colnames(bicluster_rows)
} else {
    cluster_names <- colnames(bicluster_columns)
}
annotation_colors <- list()
for (name in cluster_names) {
    annotation_colors[[name]] <- c("TRUE" = "red", "FALSE" = "white")
}

# Create a heatmap
pheatmap(corr, cluster_cols = F, cluster_rows = F,
         annotation_col = bicluster_columns, annotation_row = bicluster_rows,
         annotation_colors = annotation_colors)

