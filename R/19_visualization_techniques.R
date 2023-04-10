<<<<<<< HEAD
## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ---- include = FALSE---------------------------------------------------------
library(ggplot2)
theme_set(theme_classic())
library(mia)
library(scater)
library(patchwork)
library(miaViz)
library(sechm)
library(reshape2)
library(pheatmap)
library(ape)
library(ggtree)
# essential data
data("GlobalPatterns", package = "mia")
tse <- GlobalPatterns


## -----------------------------------------------------------------------------
# list row meta data
names(rowData(tse))
# list column meta data
names(colData(tse))


## ---- warning = FALSE---------------------------------------------------------
# obtain QC data
tse <- addPerCellQC(tse)
tse <- addPerFeatureQC(tse)
# plot QC Mean against Species
plotRowData(tse, "mean", "Species") +
  theme(axis.text.x = element_blank()) +
  labs(x = "Species", y = "QC Mean")
# plot QC Sum against Sample ID, colour-labeled by Sample Type
plotColData(tse, "sum", "X.SampleID", colour_by = "SampleType") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample ID", y = "QC Sum")


## -----------------------------------------------------------------------------
# store colData into a data frame
coldata <- as.data.frame(colData(tse))
# plot Number of Samples against Sampling Site
ggplot(coldata, aes(x = SampleType)) +
  geom_bar(width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sampling Site",
       y = "Number of Samples")


## -----------------------------------------------------------------------------
# estimate shannon diversity index
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "shannon", 
                              name = "shannon")
# plot shannon diversity index, colour-labeled by Sample Type
plotColData(tse, "shannon", colour_by = "SampleType")


## -----------------------------------------------------------------------------
# estimate faith diversity index
tse <- mia::estimateFaith(tse,
                          assay_name = "counts")
# store colData into a data frame
coldata <- as.data.frame(colData(tse))
# generate plots for shannon and faith indices
# and store them into a list
plots <- lapply(c("shannon", "faith"),
                function(i) ggplot(coldata, aes_string(y = i)) +
                  geom_boxplot() +
                  theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank()))
# combine plots with patchwork
plots[[1]] + plots[[2]]


## -----------------------------------------------------------------------------
# perform NMDS coordination method
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               name = "NMDS")
# plot results of a 2-component NMDS on tse,
# coloured-scaled by shannon diversity index
plotReducedDim(tse, "NMDS", colour_by = "shannon")


## -----------------------------------------------------------------------------
# perform MDS coordination method
tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "bray",
              name = "MDS",
              exprs_values = "counts",
              ncomponents = 3)
# plot results of a 3-component MDS on tse,
# coloured-scaled by faith diversity index
plotReducedDim(tse, "MDS", ncomponents = c(1:3), colour_by = "faith")


## -----------------------------------------------------------------------------
# generate plots for MDS and NMDS methods
# and store them into a list
plots <- lapply(c("MDS", "NMDS"),
                plotReducedDim,
                object = tse,
                colour_by = "shannon")
# combine plots with patchwork
plots[[1]] + plots[[2]] +
  plot_layout(guides = "collect")


## ----plotAbundance1-----------------------------------------------------------
# agglomerate tse by Order
tse_order <- agglomerateByRank(tse,
                                rank = "Order",
                                onRankOnly = TRUE)
# transform counts into relative abundance
tse_order <- transformCounts(tse_order,
                              assay_name = "counts",
                              method = "relabundance")
# get top orders
top_taxa <- getTopTaxa(tse_order,
                       top = 10,
                       assay_name = "relabundance")
# leave only names for top 10 orders and label the rest with "Other"
order_renamed <- lapply(rowData(tse_order)$Order,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_order)$Order <- as.character(order_renamed)
# plot composition as a bar plot
plotAbundance(tse_order,
              assay_name = "relabundance",
              rank = "Order",
              order_rank_by = "abund",
              order_sample_by = "Clostridiales")


## ----plotAbundance2-----------------------------------------------------------
# Create plots
plots <- plotAbundance(tse_order,
      	    assay_name = "relabundance",
	    rank = "Order",
            order_rank_by = "abund",
	    order_sample_by = "Clostridiales",
            features = "SampleType")

# Modify the legend of the first plot to be smaller 
plots[[1]] <- plots[[1]] +
    theme(legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))

# Modify the legend of the second plot to be smaller 
plots[[2]] <- plots[[2]] +
    theme(legend.key.height = unit(0.3, 'cm'),
          legend.key.width = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.direction = "vertical")

# Load required packages
if( !require("ggpubr") ){
    install.packages("ggpubr")
    library("ggpubr")
}
# Load required packages
if( !require("patchwork") ){
    install.packages("patchwork")
    library("patchwork")
}

# Combine legends
legend <- wrap_plots(as_ggplot(get_legend(plots[[1]])), as_ggplot(get_legend(plots[[2]])), ncol = 1) 

# Remove legends from the plots
plots[[1]] <- plots[[1]] + theme(legend.position = "none")
plots[[2]] <- plots[[2]] + theme(legend.position = "none", axis.title.x=element_blank()) 

# Combine plots
plot <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend
wrap_plots(plot, legend, nrow = 1, widths = c(2, 1))


## ----pheatmap1----------------------------------------------------------------
# Agglomerate tse by phylum
tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_phylum <- transformCounts(tse_phylum, MARGIN = "samples", method = "clr", assay_name = "counts", pseudocount=1)

# Add z-transformation on features (taxa)
tse_phylum <- transformCounts(tse_phylum, assay_name = "clr",
                              MARGIN = "features", 
                              method = "z", name = "clr_z")

# Take subset: only samples from feces, skin, or tongue
tse_phylum_subset <- tse_phylum[ , tse_phylum$SampleType %in% c("Feces", "Skin", "Tongue") ]

# Add clr-transformation
tse_phylum_subset <- transformCounts(tse_phylum_subset, method = "clr",
                                     MARGIN="samples",
                                     assay_name = "counts", pseudocount=1)
# Does z-transformation
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


## ----pheatmap2----------------------------------------------------------------
# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

# to view the tree, run
# taxa_tree


## ----pheatmap3----------------------------------------------------------------
# Creates clusters
taxa_clusters <- cutree(tree = taxa_hclust, k = 3)

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters


## ----pheatmap4----------------------------------------------------------------
# Adds information to rowData
rowData(tse_phylum_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_phylum_subset))), ]

# Prints taxa and their clusters
rowData(tse_phylum_subset)$clusters


## ----pheatmap5----------------------------------------------------------------
# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
# sample_tree

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


## ----pheatmap6----------------------------------------------------------------
# Determines the scaling of colorss
# Scale colors
breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)

pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         breaks = breaks,
         color = colors)


## ----sechm--------------------------------------------------------------------
# Stores annotation colros to metadata
metadata(tse_phylum_subset)$anno_colors$SampleType <- c(Feces = "blue", 
                                                        Skin = "red", 
                                                        Tongue = "gray")

# Create a plot
sechm(tse_phylum_subset, 
      features = rownames(tse_phylum_subset), 
      assayName = "clr", 
      do.scale = TRUE, 
      top_annotation = c("SampleType"), 
      gaps_at = "SampleType",
      cluster_cols = TRUE, cluster_rows = TRUE)


## ----more_complex_heatmap-----------------------------------------------------
# Add feature names to column as a factor
taxa_clusters$Feature <- rownames(taxa_clusters)
taxa_clusters$Feature <- factor(taxa_clusters$Feature, levels = taxa_clusters$Feature)

# Create annotation plot
row_annotation <- ggplot(taxa_clusters) + 
  geom_tile(aes(x = NA, y = Feature, fill = clusters)) +
  coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin=margin(0,0,0,0),
        ) +
      labs(fill = "Clusters", x = "Clusters")

# to view the notation, run
# row_annotation

# Add sample names to one of the columns
sample_data$sample <- factor(rownames(sample_data), levels = rownames(sample_data))

# Create annotation plot
sample_types_annotation <- ggplot(sample_data) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  geom_tile(aes(y = NA, x = sample, fill = sample_types)) +
  coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(0,0,0,0),
        axis.title.y.right = element_text(angle=0, vjust = 0.5)
        ) +
      labs(fill = "Sample types", y = "Sample types")
# to view the notation, run
# sample_types_annotation

# Create annotation plot
sample_clusters_annotation <- ggplot(sample_data) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  geom_tile(aes(y = NA, x = sample, fill = clusters)) +
  coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(0,0,0,0),
        axis.title.y.right = element_text(angle=0, vjust = 0.5)
        ) +
      labs(fill = "Clusters", y = "Clusters")
# to view the notation, run
# sample_clusters_annotation

# Order data based on clusters and sample types
mat <- mat[unfactor(taxa_clusters$Feature), unfactor(sample_data$sample)]

# ggplot requires data in melted format
melted_mat <- melt(mat)
colnames(melted_mat) <- c("Taxa", "Sample", "clr_z")

# Determines the scaling of colorss
maxval <- round(max(abs(melted_mat$clr_z)))
limits <- c(-maxval, maxval)
breaks <- seq(from = min(limits), to = max(limits), by = 0.5)
colours <- c("darkblue", "blue", "white", "red", "darkred")

heatmap <- ggplot(melted_mat) + 
  geom_tile(aes(x = Sample, y = Taxa, fill = clr_z)) +
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    
    plot.margin=margin(0,0,0,0), # removes margins
    legend.key.height= unit(1, 'cm')
    ) +
  scale_fill_gradientn(name = "CLR + Z transform", 
                       breaks = breaks, 
                       limits = limits, 
                       colours = colours) + 
  scale_y_discrete(position = "right")

heatmap


## ----more_complex_heatmap2, fig.width = 10, fig.height = 8, eval=FALSE--------
## library(patchwork)
## 
## # Create layout
## design <- c(
##   patchwork::area(3, 1, 4, 1),
##   patchwork::area(1, 2, 1, 3),
##   patchwork::area(2, 2, 2, 3),
##   patchwork::area(3, 2, 4, 3)
## )
## # to view the design, run
## # plot(design)
## 
## # Combine plots
## plot <- row_annotation + sample_clusters_annotation +
##                          sample_types_annotation +
## 			 heatmap  +
##     plot_layout(design = design, guides = "collect",
##                 # Specify layout, collect legends
## 
##                 # Adjust widths and heights to align plots.
##                 # When annotation plot is larger, it might not fit into
## 		# its column/row.
##                 # Then you need to make column/row larger.
## 
##                 # Relative widths and heights of each column and row:
##                 # Currently, the width of the first column is 15 % and the height of
##                 # first two rows are 30 % the size of others
## 
##                 # To get this work most of the times, you can adjust all sizes to be 1, i.e. equal,
##                 # but then the gaps between plots are larger.
##                 widths = c(0.15, 1, 1),
##                 heights = c(0.3, 0.3, 1, 1))
## 
## # plot


## ----more_complex_heatmap3, fig.width = 10, fig.height = 8, eval=FALSE--------
## # Create layout
## design <- c(
##   patchwork::area(4, 1, 5, 1),
##   patchwork::area(4, 2, 5, 2),
##   patchwork::area(1, 3, 1, 4),
##   patchwork::area(2, 3, 2, 4),
##   patchwork::area(3, 3, 3, 4),
##   patchwork::area(4, 3, 5, 4)
## )
## 
## # to view the design, run
## # plot(design)
## 
## # Combine plots
## plot <- taxa_tree +
##   row_annotation +
##   sample_tree +
##   sample_clusters_annotation +
##   sample_types_annotation +
##   heatmap +
##     plot_layout(design = design, guides = "collect", # Specify layout, collect legends
##                 widths = c(0.2, 0.15, 1, 1, 1),
##                 heights = c(0.1, 0.15, 0.15, 0.25, 1, 1))
## 
## plot


## ----sessionInfo, echo = FALSE, results = "asis"------------------------------
prettySessionInfo()

||||||| merged common ancestors
=======
## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ---- include = FALSE---------------------------------------------------------
library(ggplot2)
theme_set(theme_classic())
library(mia)
library(scater)
library(patchwork)
library(miaViz)
library(sechm)
library(reshape2)
library(pheatmap)
library(ape)
library(ggtree)
# essential data
data("GlobalPatterns", package = "mia")
tse <- GlobalPatterns


## -----------------------------------------------------------------------------
# list row meta data
names(rowData(tse))
# list column meta data
names(colData(tse))


## ---- warning = FALSE---------------------------------------------------------
# obtain QC data
tse <- addPerCellQC(tse)
tse <- addPerFeatureQC(tse)
# plot QC Mean against Species
plotRowData(tse, "mean", "Species") +
  theme(axis.text.x = element_blank()) +
  labs(x = "Species", y = "QC Mean")
# plot QC Sum against Sample ID, colour-labeled by Sample Type
plotColData(tse, "sum", "X.SampleID", colour_by = "SampleType") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sample ID", y = "QC Sum")


## -----------------------------------------------------------------------------
# store colData into a data frame
coldata <- as.data.frame(colData(tse))
# plot Number of Samples against Sampling Site
ggplot(coldata, aes(x = SampleType)) +
  geom_bar(width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Sampling Site",
       y = "Number of Samples")


## -----------------------------------------------------------------------------
# estimate shannon diversity index
tse <- mia::estimateDiversity(tse, 
                              assay_name = "counts",
                              index = "shannon", 
                              name = "shannon")
# plot shannon diversity index, colour-labeled by Sample Type
plotColData(tse, "shannon", colour_by = "SampleType")


## -----------------------------------------------------------------------------
# estimate faith diversity index
tse <- mia::estimateFaith(tse,
                          assay_name = "counts")
# store colData into a data frame
coldata <- as.data.frame(colData(tse))
# generate plots for shannon and faith indices
# and store them into a list
plots <- lapply(c("shannon", "faith"),
                function(i) ggplot(coldata, aes_string(y = i)) +
                  geom_boxplot() +
                  theme(axis.text.x = element_blank(),
                        axis.ticks.x = element_blank()))
# combine plots with patchwork
plots[[1]] + plots[[2]]


## -----------------------------------------------------------------------------
# perform NMDS coordination method
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               name = "NMDS")
# plot results of a 2-component NMDS on tse,
# coloured-scaled by shannon diversity index
plotReducedDim(tse, "NMDS", colour_by = "shannon")


## -----------------------------------------------------------------------------
# perform MDS coordination method
tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "bray",
              name = "MDS",
              exprs_values = "counts",
              ncomponents = 3)
# plot results of a 3-component MDS on tse,
# coloured-scaled by faith diversity index
plotReducedDim(tse, "MDS", ncomponents = c(1:3), colour_by = "faith")


## -----------------------------------------------------------------------------
# generate plots for MDS and NMDS methods
# and store them into a list
plots <- lapply(c("MDS", "NMDS"),
                plotReducedDim,
                object = tse,
                colour_by = "shannon")
# combine plots with patchwork
plots[[1]] + plots[[2]] +
  plot_layout(guides = "collect")


## ----plotAbundance1-----------------------------------------------------------
# agglomerate tse by Order
tse_order <- agglomerateByRank(tse,
                                rank = "Order",
                                onRankOnly = TRUE)
# transform counts into relative abundance
tse_order <- transformCounts(tse_order,
                              assay_name = "counts",
                              method = "relabundance")
# get top orders
top_taxa <- getTopTaxa(tse_order,
                       top = 10,
                       assay_name = "relabundance")
# leave only names for top 10 orders and label the rest with "Other"
order_renamed <- lapply(rowData(tse_order)$Order,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse_order)$Order <- as.character(order_renamed)
# plot composition as a bar plot
plotAbundance(tse_order,
              assay_name = "relabundance",
              rank = "Order",
              order_rank_by = "abund",
              order_sample_by = "Clostridiales")


## ----plotAbundance2-----------------------------------------------------------
# Create plots
plots <- plotAbundance(tse_order,
      	    assay_name = "relabundance",
	    rank = "Order",
            order_rank_by = "abund",
	    order_sample_by = "Clostridiales",
            features = "SampleType")

# Modify the legend of the first plot to be smaller 
plots[[1]] <- plots[[1]] +
    theme(legend.key.size = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8))

# Modify the legend of the second plot to be smaller 
plots[[2]] <- plots[[2]] +
    theme(legend.key.height = unit(0.3, 'cm'),
          legend.key.width = unit(0.3, 'cm'),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.direction = "vertical")

# Load required packages
if( !require("ggpubr") ){
    install.packages("ggpubr")
    library("ggpubr")
}
# Load required packages
if( !require("patchwork") ){
    install.packages("patchwork")
    library("patchwork")
}

# Combine legends
legend <- wrap_plots(as_ggplot(get_legend(plots[[1]])), as_ggplot(get_legend(plots[[2]])), ncol = 1) 

# Remove legends from the plots
plots[[1]] <- plots[[1]] + theme(legend.position = "none")
plots[[2]] <- plots[[2]] + theme(legend.position = "none", axis.title.x=element_blank()) 

# Combine plots
plot <- wrap_plots(plots[[2]], plots[[1]], ncol = 1, heights = c(2, 10))
# Combine the plot with the legend
wrap_plots(plot, legend, nrow = 1, widths = c(2, 1))


## ----pheatmap1----------------------------------------------------------------
# Agglomerate tse by phylum
tse_phylum <- agglomerateByRank(tse,
                                rank = "Phylum",
                                onRankOnly = TRUE)

# Add clr-transformation on samples
tse_phylum <- transformCounts(tse_phylum, MARGIN = "samples", method = "clr", assay_name = "counts", pseudocount=1)

# Add z-transformation on features (taxa)
tse_phylum <- transformCounts(tse_phylum, assay_name = "clr",
                              MARGIN = "features", 
                              method = "z", name = "clr_z")

# Take subset: only samples from feces, skin, or tongue
tse_phylum_subset <- tse_phylum[ , tse_phylum$SampleType %in% c("Feces", "Skin", "Tongue") ]

# Add clr-transformation
tse_phylum_subset <- transformCounts(tse_phylum_subset, method = "clr",
                                     MARGIN="samples",
                                     assay_name = "counts", pseudocount=1)
# Does z-transformation
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


## ----pheatmap2----------------------------------------------------------------
# Hierarchical clustering
taxa_hclust <- hclust(dist(mat), method = "complete")

# Creates a phylogenetic tree
taxa_tree <- as.phylo(taxa_hclust)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

# to view the tree, run
# taxa_tree


## ----pheatmap3----------------------------------------------------------------
# Creates clusters
taxa_clusters <- cutree(tree = taxa_hclust, k = 3)

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters


## ----pheatmap4----------------------------------------------------------------
# Adds information to rowData
rowData(tse_phylum_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_phylum_subset))), ]

# Prints taxa and their clusters
rowData(tse_phylum_subset)$clusters


## ----pheatmap5----------------------------------------------------------------
# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

# to view the tree, run
# sample_tree

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


## ----pheatmap6----------------------------------------------------------------
# Determines the scaling of colorss
# Scale colors
breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)

pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         breaks = breaks,
         color = colors)


## ----sechm--------------------------------------------------------------------
# Stores annotation colros to metadata
metadata(tse_phylum_subset)$anno_colors$SampleType <- c(Feces = "blue", 
                                                        Skin = "red", 
                                                        Tongue = "gray")

# Create a plot
sechm(tse_phylum_subset, 
      features = rownames(tse_phylum_subset), 
      assayName = "clr", 
      do.scale = TRUE, 
      top_annotation = c("SampleType"), 
      gaps_at = "SampleType",
      cluster_cols = TRUE, cluster_rows = TRUE)


## ----more_complex_heatmap-----------------------------------------------------
# Add feature names to column as a factor
taxa_clusters$Feature <- rownames(taxa_clusters)
taxa_clusters$Feature <- factor(taxa_clusters$Feature, levels = taxa_clusters$Feature)

# Create annotation plot
row_annotation <- ggplot(taxa_clusters) + 
  geom_tile(aes(x = NA, y = Feature, fill = clusters)) +
  coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin=margin(0,0,0,0),
        ) +
      labs(fill = "Clusters", x = "Clusters")

# to view the notation, run
# row_annotation

# Add sample names to one of the columns
sample_data$sample <- factor(rownames(sample_data), levels = rownames(sample_data))

# Create annotation plot
sample_types_annotation <- ggplot(sample_data) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  geom_tile(aes(y = NA, x = sample, fill = sample_types)) +
  coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(0,0,0,0),
        axis.title.y.right = element_text(angle=0, vjust = 0.5)
        ) +
      labs(fill = "Sample types", y = "Sample types")
# to view the notation, run
# sample_types_annotation

# Create annotation plot
sample_clusters_annotation <- ggplot(sample_data) +
  scale_y_discrete(position = "right", expand = c(0,0)) +
  geom_tile(aes(y = NA, x = sample, fill = clusters)) +
  coord_equal(ratio = 1) +
  theme(
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=margin(0,0,0,0),
        axis.title.y.right = element_text(angle=0, vjust = 0.5)
        ) +
      labs(fill = "Clusters", y = "Clusters")
# to view the notation, run
# sample_clusters_annotation

# Order data based on clusters and sample types
mat <- mat[unfactor(taxa_clusters$Feature), unfactor(sample_data$sample)]

# ggplot requires data in melted format
melted_mat <- melt(mat)
colnames(melted_mat) <- c("Taxa", "Sample", "clr_z")

# Determines the scaling of colorss
maxval <- round(max(abs(melted_mat$clr_z)))
limits <- c(-maxval, maxval)
breaks <- seq(from = min(limits), to = max(limits), by = 0.5)
colours <- c("darkblue", "blue", "white", "red", "darkred")

heatmap <- ggplot(melted_mat) + 
  geom_tile(aes(x = Sample, y = Taxa, fill = clr_z)) +
  theme(
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.y=element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    
    plot.margin=margin(0,0,0,0), # removes margins
    legend.key.height= unit(1, 'cm')
    ) +
  scale_fill_gradientn(name = "CLR + Z transform", 
                       breaks = breaks, 
                       limits = limits, 
                       colours = colours) + 
  scale_y_discrete(position = "right")

heatmap


## ----more_complex_heatmap2, fig.width = 10, fig.height = 8--------------------
library(patchwork)

# Create layout
design <- c(
  patchwork::area(3, 1, 4, 1),
  patchwork::area(1, 2, 1, 3),
  patchwork::area(2, 2, 2, 3),
  patchwork::area(3, 2, 4, 3)
)
# to view the design, run
# plot(design)

# Combine plots
plot <- row_annotation + sample_clusters_annotation +
                         sample_types_annotation +
			 heatmap  +
    plot_layout(design = design, guides = "collect",
                # Specify layout, collect legends
                
                # Adjust widths and heights to align plots.
                # When annotation plot is larger, it might not fit into
		# its column/row.
                # Then you need to make column/row larger.
                
                # Relative widths and heights of each column and row:
                # Currently, the width of the first column is 15 % and the height of
                # first two rows are 30 % the size of others
                
                # To get this work most of the times, you can adjust all sizes to be 1, i.e. equal, 
                # but then the gaps between plots are larger.
                widths = c(0.15, 1, 1),
                heights = c(0.3, 0.3, 1, 1))

# plot


## ----more_complex_heatmap3, fig.width = 10, fig.height = 8--------------------
# Create layout
design <- c(
  patchwork::area(4, 1, 5, 1),
  patchwork::area(4, 2, 5, 2),
  patchwork::area(1, 3, 1, 4),
  patchwork::area(2, 3, 2, 4),
  patchwork::area(3, 3, 3, 4),
  patchwork::area(4, 3, 5, 4)
)

# to view the design, run
# plot(design)

# Combine plots
plot <- taxa_tree + 
  row_annotation +
  sample_tree + 
  sample_clusters_annotation +
  sample_types_annotation +
  heatmap +
    plot_layout(design = design, guides = "collect", # Specify layout, collect legends
                widths = c(0.2, 0.15, 1, 1, 1),
                heights = c(0.1, 0.15, 0.15, 0.25, 1, 1))

plot


## ----sessionInfo, echo = FALSE, results = "asis"------------------------------
prettySessionInfo()

>>>>>>> a4695d89c40cde94395c726a99a6073cf631e43a
