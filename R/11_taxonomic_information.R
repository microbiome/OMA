## ----setup, echo=FALSE, results="asis"---------------------------------------------------------------------------------
library(rebook)
chapterPreamble()


## ---- message=FALSE----------------------------------------------------------------------------------------------------
library(mia)
data("GlobalPatterns", package = "mia")
tse <- GlobalPatterns


## ----------------------------------------------------------------------------------------------------------------------
checkTaxonomy(tse)


## ----------------------------------------------------------------------------------------------------------------------
taxonomyRanks(tse)


## ----------------------------------------------------------------------------------------------------------------------
rowData(tse)[, taxonomyRanks(tse)]


## ----------------------------------------------------------------------------------------------------------------------
all(!taxonomyRankEmpty(tse, rank = "Kingdom"))
table(taxonomyRankEmpty(tse, rank = "Genus"))
table(taxonomyRankEmpty(tse, rank = "Species"))


## ----------------------------------------------------------------------------------------------------------------------
head(getTaxonomyLabels(tse))


## ----------------------------------------------------------------------------------------------------------------------
phylum <- !is.na(rowData(tse)$Phylum) &
    vapply(data.frame(apply(rowData(tse)[, taxonomyRanks(tse)[3:7]], 1L, is.na)), all, logical(1))
head(getTaxonomyLabels(tse[phylum,]))
head(getTaxonomyLabels(tse[phylum,], with_rank = TRUE))


## ----------------------------------------------------------------------------------------------------------------------
head(getTaxonomyLabels(tse[phylum,], with_rank = TRUE, make_unique = FALSE))


## ----------------------------------------------------------------------------------------------------------------------
head(getUniqueFeatures(tse, rank = "Phylum"))


## ----------------------------------------------------------------------------------------------------------------------
taxonomyTree(tse)


## ----------------------------------------------------------------------------------------------------------------------
tse <- addTaxonomyTree(tse)
tse


## ----------------------------------------------------------------------------------------------------------------------
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
altExp(tse, "Family") <- mergeFeaturesByRank(tse, rank = "Family",
                                           agglomerateTree = TRUE)
altExp(tse, "Family")


## ----------------------------------------------------------------------------------------------------------------------
assayNames(tse)
assayNames(altExp(tse, "Family"))


## ----------------------------------------------------------------------------------------------------------------------
assay(altExp(tse, "Family"), "relabundance")[1:5, 1:7]


## ----taxinfo_altexp_example--------------------------------------------------------------------------------------------
assay(altExp(tse, "Family"), "counts")[1:5, 1:7]


## ----------------------------------------------------------------------------------------------------------------------
altExp(tse, "Species_byPrevalence") <- mergeFeaturesByPrevalence(tse, 
                                                               rank = "Species", 
                                                               other_label = "Other", 
                                                               prevalence = 5 / 100, 
                                                               detection = 1 / 100, 
                                                               as_relative = T)
altExp(tse, "Species_byPrevalence")

assay(altExp(tse, "Species_byPrevalence"), "relabundance")[88:92, 1:7]


## ----------------------------------------------------------------------------------------------------------------------
# Saving the tse for later
tseGlobalPatterns <- tse


## ----bluster_dependence------------------------------------------------------------------------------------------------
library(bluster)


## ----taxa_clustering---------------------------------------------------------------------------------------------------
# Get the data
data("peerj13075", package = "mia")
tse <- peerj13075

# The result of the CLR transform is stored in the assay clr
tse <- transformAssay(tse, method = "clr", pseudocount = 1)

tse <- transformAssay(tse, assay.type = "clr", method = "z", 
                      MARGIN = "features")

# Cluster (with euclidean distance) on the features of the z assay
tse <- cluster(tse,
               assay.type = "z",
               clust.col = "hclustEuclidean",
	       MARGIN = "features",
               HclustParam(dist.fun = stats::dist, method = "ward.D2"))

# Declare the Kendall dissimilarity computation function
kendall_dissimilarity <- function(x) {
    as.dist(1 - cor(t(x), method = "kendall"))
}

# Cluster (with Kendall dissimilarity) on the features of the z assay
tse <- cluster(tse,
               assay.type = "z",
               clust.col = "hclustKendall",
       	       MARGIN = "features", 	       
               HclustParam(dist.fun = kendall_dissimilarity, method = "ward.D2"))


## ----taxa_clustering_result--------------------------------------------------------------------------------------------
# Checking the clusters
clusters_euclidean <- rowData(tse)$hclustEuclidean
head(clusters_euclidean, 10)

clusters_kendall <- rowData(tse)$hclustKendall
head(clusters_kendall, 10)


## ----taxa_clustering_histogram-----------------------------------------------------------------------------------------
library(ggplot2)
library(patchwork) # TO arrange several plots as a grid
plot1 <- ggplot(as.data.frame(rowData(tse)), aes(x = clusters_euclidean)) +
    geom_bar() +
    labs(title = "CAG size distribution (Euclidean distance)",
         x = "Clusters", y = "Feature count (n)")
plot2 <- ggplot(as.data.frame(rowData(tse)), aes(x = clusters_kendall)) +
    geom_bar() +
    labs(title = "CAG size distribution (1 - tau)",
         x = "Clusters", y = "Feature count (n)")
plot1 + plot2 + plot_layout(ncol = 2)


## ----taxa_clustering_row_merge-----------------------------------------------------------------------------------------
# Aggregate clusters as a sum of each cluster values
tse_merged <- mergeFeatures(tse, clusters_euclidean)
tse_merged


## ----------------------------------------------------------------------------------------------------------------------
tse <- tseGlobalPatterns
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance", pseudocount = 1)
tse <- transformAssay(x = tse, assay.type = "relabundance", method = "clr", 
                      pseudocount = 1, name = "clr")

head(assay(tse, "clr"))


## ----------------------------------------------------------------------------------------------------------------------
tse <- transformAssay(tse, method = "pa")

head(assay(tse, "pa"))


## ----------------------------------------------------------------------------------------------------------------------
# list of abundance tables that assays slot contains
assays(tse)

