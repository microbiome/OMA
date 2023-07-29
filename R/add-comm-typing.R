## ----data, message=FALSE, warning=FALSE---------------------------------------
# Load data
library(mia)
library(microbiomeDataSets)
tse <- SprockettTHData()


## ---- message=FALSE, warning=FALSE--------------------------------------------
library(miaViz)

# Only consider Forest samples
tse <- tse[, colData(tse)$Country == "Finland"]
# Agglomerate to phylum level
tse <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree=TRUE)
# Get top taxa
rel_taxa <- getTopFeatures(tse, top = 8, assay.type = "counts")
# Take only the top taxa
tse <- tse[is.element(rownames(tse), rel_taxa), ]

# Visualise composition barplot, with samples order by "Firmicutes"
plotAbundance(tse, rank = "Phylum", order_rank_by = "abund", order_sample_by = "Firmicutes")


## ---- message=FALSE, warning=FALSE--------------------------------------------
library(pheatmap)
library(grid)
library(RColorBrewer)
# Agglomerate to phylum level
tse <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree=TRUE)
# Take only the top taxa
tse <- tse[is.element(rownames(tse), rel_taxa),]

# CLR and Z transforms
tse <- transformAssay(tse, MARGIN = "samples", assay.type = "counts",
                       method = "clr", pseudocount=1)
tse <- transformAssay(tse, MARGIN = "features", assay.type = "clr", method = "z")

Countries <- data.frame("Country" = colData(tse)$Country)
rownames(Countries) <- colData(tse)$Sample_ID
# count matrix for pheatmap
mat <- assays(tse)$z

# Order community types
mat <- mat[, order(Countries$Country)]
colnames(mat) <- colnames(mat)[order(Countries$Country)]
rownames(mat) <- stringr::str_remove(rownames(mat), "Phylum:")
# Make grid for heatmap
breaks <- seq(-3, 3, length.out = 10)
setHook("grid.newpage", function() pushViewport(viewport(x = 1, y = 1, width = 0.9, 
                                                         height = 0.9, name = "vp", 
                                                         just = c("right","top"))), 
        action = "prepend")
pheatmap(mat, color = rev(brewer.pal(9, "RdBu")), breaks = breaks, main = "Countries", treeheight_row = 0, treeheight_col = 0, show_colnames = 0, annotation_col = Countries, cluster_cols = F)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample", x = 0.39, y = -0.04, gp = gpar(fontsize = 16))
grid.text("Phylum", x = -0.04, y = 0.47, rot = 90, gp = gpar(fontsize = 16))


## ----scree--------------------------------------------------------------------
library(ggplot2); th <- theme_bw()

# Only consider Finland samples
tse <- tse[, colData(tse)$Country == "Finland"]

# MDS analysis with the first 20 dimensions
tse  <- scater::runMDS(tse, FUN = vegan::vegdist, method = "bray", 
                       name = "MDS_BC", assay.type = "counts", ncomponents = 20)
ord  <- reducedDim(tse, "MDS_BC", withDimnames = TRUE)
# retrieve eigenvalues
eigs <- attr(ord, "eig")

# variance in each of the axes
var <- eigs / sum(eigs)
# first 12 values
df <- data.frame(x = c(1:12), y = var[1:12])
# create scree plot
p <- ggplot(df, aes(x, y)) +
     geom_bar(stat="identity") +
     xlab("Principal Component") +
     ylab("Variance") +
     ggtitle("Scree Plot")
p


## -----------------------------------------------------------------------------
# histogram of MDS eigenvalues from the fifth dimension onward
h <- hist(eigs[3:length(eigs)], 100)


## ----message = FALSE, warning = FALSE-----------------------------------------
plot(h$mids, h$count, log = "y", type = "h", lwd = 10, lend = 2)


## ----elbow, message = FALSE---------------------------------------------------
library(factoextra)

# take only first 5 dimensions
NDIM <- 5
x    <- ord[, 1:NDIM]

# Elbow Method
factoextra::fviz_nbclust(x, kmeans, method = "wss") +
                         geom_vline(xintercept = 3, linetype = 2) +
                         labs(subtitle = "Elbow Method") + th


## ----silhouette---------------------------------------------------------------
# Silhouette method
factoextra::fviz_nbclust(x, kmeans, method = "silhouette") +
                         labs(subtitle = "Silhouette method") + th


## ----gap-statistic------------------------------------------------------------
# Gap Statistic Method
factoextra::fviz_nbclust(x, kmeans, method = "gap_stat", nboot = 50)+
                         labs(subtitle = "Gap Statistic Method") + th


## ----create clusters----------------------------------------------------------
library(cluster)
tse2 <- agglomerateByRank(tse, rank = "Phylum", agglomerateTree=TRUE)

# assume 6 clusters
K <- 6
x <- ord[, 1:NDIM]

clust <- as.factor(pam(x, k = K, cluster.only = T))
# assign CSTs
colData(tse2)$CST <- clust
CSTs <- as.character(seq(K))


## ----message = FALSE, warning = FALSE-----------------------------------------
library(scater)
library(RColorBrewer)
library(patchwork)

# set up colours
CSTColors <- brewer.pal(6, "Paired")[c(2, 5, 3, 4, 1, 6)]
names(CSTColors) <- CSTs

CSTColorScale <- scale_colour_manual(name = "CST", values = CSTColors)
CSTFillScale <- scale_fill_manual(name = "CST", values = CSTColors)

# plot MDS with Bray-Curtis dimensions 1 and 2
p1 <- scater::plotReducedDim(tse2, "MDS_BC", colour_by = "CST", point_alpha = 1, 
                             percentVar = var[c(1, 2)]*100) + th + labs(title = "Ordination by Cluster") +
                             theme(plot.title = element_text(hjust = 0.5))
# plot MDS with Bray-Curtis dimensions 3 and 4
p2 <- scater::plotReducedDim(tse2, "MDS_BC", colour_by = "CST", point_alpha = 1, 
                             ncomponents = c(3, 4), percentVar = var[c(1, 2, 3, 4)]*100) + th
# show results
(p1 + CSTColorScale) / (p2 + CSTColorScale)


## ----message = FALSE, warning = FALSE-----------------------------------------
tse2  <- runNMDS(tse2, FUN = vegan::vegdist, method = "bray", 
                name = "NMDS_BC", assay.type = "counts", ncomponents = 20)
scater::plotReducedDim(tse2, "NMDS_BC", colour_by = "CST", point_alpha = 1) + th + 
        labs(title = "NMDS Bray-Curtis by Cluster") +
        theme(plot.title = element_text(hjust = 0.5)) + CSTColorScale


## ----message = FALSE, warning = FALSE, results = FALSE------------------------
# Z transform of CLR counts
tse2 <- transformAssay(tse2, MARGIN = "samples", assay.type = "counts",
                        method = "clr", pseudocount=1)
tse2 <- transformAssay(tse2, MARGIN = "features", assay.type = "clr", method = "z")
# get top taxa
tse2 <- tse2[is.element(rownames(tse2), rel_taxa), ]

mat <- assays(tse2)$z

# Order CSTs
mat <- mat[, order(clust)]
colnames(mat) <- names(sort(clust))
rownames(mat) <- stringr::str_remove(rownames(mat), "Phylum:")


## ----messages = FALSE, warning = FALSE----------------------------------------
# Plot
CSTs        <- as.data.frame(sort(clust))
names(CSTs) <- "CST"
breaks <- seq(-2, 2, length.out = 10)
# Make grid for heatmap
setHook("grid.newpage", function() pushViewport(viewport(x = 1, y = 1, width = 0.9, 
                                                         height = 0.9, name = "vp", 
                                                         just = c("right","top"))), 
        action = "prepend")
pheatmap(mat, color = rev(brewer.pal(9, "RdBu")), breaks = breaks, main = "All CSTs", treeheight_row = 0, treeheight_col = 0, show_colnames = 0, annotation_col = CSTs, cluster_cols = F)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample", x = 0.39, y = -0.04, gp = gpar(fontsize = 16))
grid.text("Phylum", x = -0.04, y = 0.47, rot = 90, gp = gpar(fontsize = 16))


## ----dmm----------------------------------------------------------------------
# Runs model and calculates the most likely number of clusters from 1 to 7.
# Since this is a large dataset it takes a long time to run
# For this reason we use only a subset of the data agglomerated to the phylum level
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
dmn_group <- calculateDMNgroup(tse_dmn, variable = "Country",  assay.type = "counts",
                               k = 3, seed=.Machine$integer.max)

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
colnames(prob) <- paste0("comp", seq_len(ncol(prob)))

# For each row, finds column that has the highest value. Then extract the column 
# names of highest values.
vec <- colnames(prob)[max.col(prob,ties.method = "first")]


## -----------------------------------------------------------------------------
# Calculate relative abundances
tse_dmn <- transformAssay(tse_dmn, method = "relabundance")
# Does principal coordinate analysis
bray_pcoa_df <- calculateMDS(tse_dmn, FUN = vegan::vegdist, method = "bray", 
                             assay.type = "relabundance")

# Convert to data.frame
bray_pcoa_df <- as.data.frame(bray_pcoa_df)
colnames(bray_pcoa_df) <- c("pcoa1", "pcoa2")
head(bray_pcoa_df)


## -----------------------------------------------------------------------------
# Creates a data frame that contains principal coordinates and DMM information
bray_dmm_pcoa_df <- bray_pcoa_df 
bray_dmm_pcoa_df$dmm_component <- vec
# Creates a plot
bray_dmm_plot <- ggplot(data = bray_dmm_pcoa_df, 
                        aes(x = pcoa1, y = pcoa2, color = dmm_component)) +
  geom_point() +
  labs(x = "Coordinate 1",
       y = "Coordinate 2",
       title = "PCoA with Bray-Curtis Dissimilarity") +  
  theme(title = element_text(size = 12)) + theme_bw() # makes titles smaller

bray_dmm_plot


## -----------------------------------------------------------------------------
# get clr + z-transformed counts
tse_dmn <- transformAssay(tse_dmn, MARGIN = "samples",
                           assay.type = "counts", method = "clr",
			   pseudocount = 1)
tse_dmn <- transformAssay(tse_dmn, MARGIN = "features",
                           assay.type = "clr", method = "z")
# objects containing dmm component information
clust <- factor(vec)
names(clust) <- colnames(tse_dmn)
# get top taxa
tse_dmn <- tse_dmn[is.element(rownames(tse_dmn), rel_taxa), ]
# get just counts
mat <- assays(tse_dmn)$z
# order according to dmm component
mat <- mat[, order(clust)]
colnames(mat) <- names(sort(clust))
rownames(mat) <- stringr::str_remove(rownames(mat), "Phylum:")
# Plot
CSTs        <- as.data.frame(sort(clust))
names(CSTs) <- "CST"
breaks <- seq(-2, 2, length.out = 10)
# Make grid for heatmap
setHook("grid.newpage", function() pushViewport(viewport(x = 1, y = 1, width = 0.9, 
                                                         height = 0.9, name = "vp", 
                                                         just = c("right","top"))), 
        action = "prepend")
pheatmap(mat, color = rev(brewer.pal(9, "RdBu")), breaks = breaks, main = "All CSTs", treeheight_row = 0, treeheight_col = 0, show_colnames = 0, annotation_col = CSTs, cluster_cols = F)
setHook("grid.newpage", NULL, "replace")
grid.text("Sample", x = 0.39, y = -0.04, gp = gpar(fontsize = 16))
grid.text("Phylum", x = -0.04, y = 0.47, rot = 90, gp = gpar(fontsize = 16))


## -----------------------------------------------------------------------------
sessionInfo()

