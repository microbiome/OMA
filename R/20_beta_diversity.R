## ----setup, echo = FALSE, results = "asis"-------------
library(rebook)
chapterPreamble()


## ----include = FALSE-----------------------------------
# global knitr options
knitr::opts_chunk$set(
  message = FALSE,
  fig.width = 10,
  dpi = 300,
  dev = "png",
  dev.args = list(type = "cairo-png")
)


## ----prep-tse------------------------------------------
# Load mia and import sample dataset
library(mia)
data("GlobalPatterns", package = "mia")
tse <- GlobalPatterns

# Beta diversity metrics like Bray-Curtis are often applied to relabundances
tse <- transformAssay(tse,
                      assay.type = "counts",
                      method = "relabundance")

# Other metrics like Aitchison to clr-transformed data
tse <- transformAssay(tse,
                      assay.type = "relabundance",
                      method = "clr",
                      pseudocount = TRUE)

# Add group information Feces yes/no
tse$Group <- tse$SampleType == "Feces"


## ----runMDS--------------------------------------------
# Load package to plot reducedDim
library(scater)

# Run PCoA on relabundance assay with Bray-Curtis distances
tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "bray",
              assay.type = "relabundance",
              name = "MDS_bray")


## ----plot-mds-bray-curtis, fig.cap = "MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset."----
# Create ggplot object
p <- plotReducedDim(tse, "MDS_bray",
                    colour_by = "Group")

# Calculate explained variance
e <- attr(reducedDim(tse, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = ""))

p


## ----mds-nmds-comparison, results='hide'---------------
# Run NMDS on relabundance assay with Bray-Curtis distances
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               method = "bray",
               assay.type = "relabundance",
               name = "NMDS_bray")

# Run MDS on clr assay with Aitchison distances
tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "euclidean",
              assay.type = "clr",
              name = "MDS_aitchison")

# Run NMDS on clr assay with Euclidean distances
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               method = "euclidean",
               assay.type = "clr",
               name = "NMDS_aitchison")


## ---- fig.cap = "Comparison of MDS and NMDS plots based on the Bray-Curtis or Aitchison distances on the GlobalPattern dataset."----
# Load package for multi-panel plotting
library(patchwork)

# Generate plots for all 4 reducedDims
plots <- lapply(c("MDS_bray", "MDS_aitchison",
                  "NMDS_bray", "NMDS_aitchison"),
                plotReducedDim,
                object = tse,
                colour_by = "Group")

# Generate multi-panel plot
wrap_plots(plots) +
  plot_layout(guides = "collect")


## ----plot-unifrac, fig.cap = "Unifrac distances scaled by MDS of the GlobalPattern dataset."----
tse <- runMDS(tse,
              FUN = mia::calculateUnifrac,
              name = "Unifrac",
              tree = rowTree(tse),
              ntop = nrow(tse),
              assay.type = "counts")

plotReducedDim(tse, "Unifrac",
               colour_by = "Group")


## ----plot-pca, fig.cap = "PCA plot on the GlobalPatterns data set containing sample from different sources."----
tse <- runPCA(tse,
              name = "PCA",
              assay.type = "counts",
              ncomponents = 10)

plotReducedDim(tse, "PCA",
               colour_by = "Group")


## ----plot-umap, fig.cap = "UMAP plot on the GlobalPatterns data set containing sample from different sources."----
tse <- runUMAP(tse,
               name = "UMAP",
               assay.type = "counts",
               ncomponents = 3)

plotReducedDim(tse, "UMAP",
               colour_by = "Group",
               ncomponents = c(1:3))


## ----relstress-----------------------------------------
# Load vegan package
library(vegan)

# Quantify dissimilarities in the original feature space
x <- assay(tse, "relabundance") # Pick relabunance assay separately
d0 <- as.matrix(vegdist(t(x), "bray"))

# PCoA Ordination
pcoa <- as.data.frame(cmdscale(d0, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

# Quantify dissimilarities in the ordination space
dp <- as.matrix(dist(pcoa))

# Calculate stress i.e. relative difference in the original and
# projected dissimilarities
stress <- sum((dp - d0)^2) / sum(d0^2)


## ----shepard-------------------------------------------
ord <- order(as.vector(d0))
df <- data.frame(d0 = as.vector(d0)[ord],
                 dmds = as.vector(dp)[ord])

ggplot(df, aes(x = d0, y = dmds)) +
  geom_smooth() +
  geom_point() +    
  labs(title = "Shepard plot",
       x = "Original distance",
       y = "MDS distance",   
       subtitle = paste("Stress:", round(stress, 2))) +
  theme_bw()


## ----import-rda-dataset--------------------------------
# Load data
data("enterotype", package = "mia")
tse2 <- enterotype

# Apply relative transform
tse2 <- transformAssay(tse2,
                       method = "relabundance")


## ----run-rda-------------------------------------------
# Perform RDA
tse2 <- runRDA(tse2,
               assay.type = "relabundance",
               formula = assay ~ ClinicalStatus + Gender + Age,
               distance = "bray",
               na.action = na.exclude)

# Store results of PERMANOVA test
rda_info <- attr(reducedDim(tse2, "RDA"), "significance")


## ----rda-permanova-res---------------------------------
rda_info$permanova |>
  knitr::kable()


## ----rda-homogeneity-res-------------------------------
rda_info$homogeneity |>
  knitr::kable()


## ----plot-rda------------------------------------------
# Load packages for plotting function
library(miaViz)

# Generate RDA plot coloured by clinical status
plotRDA(tse2, "RDA", colour_by = "ClinicalStatus")


## ------------------------------------------------------
# Agglomerate to genus level
tse_genus <- mergeFeaturesByRank(tse,
                                 rank = "Genus")

# Convert to relative abundances
tse_genus <- transformAssay(tse,
                            method = "relabundance",
                            assay.type = "counts")

# Add info on dominant genus per sample
tse_genus <- addPerSampleDominantFeatures(tse_genus,
                                          assay.type = "relabundance",
                                          name = "dominant_taxa")
# Overview
countDominantFeatures(tse_genus, rank = "Genus", digits = 3, name = "dominant_taxa")


## ------------------------------------------------------
tse_genus <- runMDS(tse_genus,
                    FUN = vegan::vegdist,
                    name = "PCoA_BC",
                    assay.type = "relabundance")


## ------------------------------------------------------
# Getting the top taxa
top_taxa <- getTopFeatures(tse_genus,
                           top = 6,
                           assay.type = "relabundance")

# Naming all the rest of non top-taxa as "Other"
most_abundant <- lapply(colData(tse_genus)$dominant_taxa,
                        function(x) {if (x %in% top_taxa) {x} else {"Other"}})

# Storing the previous results as a new column within colData
colData(tse_genus)$most_abundant <- as.character(most_abundant)

# Calculating percentage of the most abundant
most_abundant_freq <- table(as.character(most_abundant))
most_abundant_percent <- round(most_abundant_freq / sum(most_abundant_freq) * 100, 1)

# Retrieving the explained variance
e <- attr(reducedDim(tse_genus, "PCoA_BC"), "eig")
var_explained <- e / sum(e[e > 0]) * 100

# Define colors for visualization
my_colors <- c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red")

# Visualization
plot <-plotReducedDim(tse_genus, "PCoA_BC",
                      colour_by = "most_abundant") +
  scale_colour_manual(values = my_colors,
                      labels = paste0(names(most_abundant_percent), "(", most_abundant_percent, "%)")) +
  labs(x = paste("PC 1 (", round(var_explained[1], 1), "%)"),
       y = paste("PC 2 (", round(var_explained[2], 1), "%)"),
       color = "")

plot


## ------------------------------------------------------
# Calculating the frequencies and percentages for both categories
freq_TRUE <- table(as.character(most_abundant[colData(tse_genus)$Group == TRUE]))
freq_FALSE <- table(as.character(most_abundant[colData(tse_genus)$Group == FALSE]))
percent_TRUE <- round(freq_TRUE / sum(freq_TRUE) * 100, 1)
percent_FALSE <- round(freq_FALSE / sum(freq_FALSE) * 100, 1)

# Visualization
plotReducedDim(tse_genus[ , colData(tse_genus)$Group == TRUE], "PCoA_BC",
               colour_by = "most_abundant") +
  scale_colour_manual(values = my_colors,
                      labels = paste0(names(percent_TRUE), "(", percent_TRUE, "%)")) +
  labs(x = paste("PC 1 (", round(var_explained[1], 1), "%)"),
       y = paste("PC 2 (", round(var_explained[2], 1), "%)"),
       title = "Group = TRUE", color = "")

plotReducedDim(tse_genus[ , colData(tse_genus)$Group == FALSE], "PCoA_BC",
               colour_by = "most_abundant") +
  scale_colour_manual(values = my_colors,
                      labels = paste0(names(percent_FALSE), "(", percent_FALSE, "%)")) +
  labs(x = paste("PC 1 (", round(var_explained[1], 1), "%)"),
       y = paste("PC 2 (", round(var_explained[2], 1), "%)"),
       title = "Group = FALSE", color = "")


## ------------------------------------------------------
# Agglomerate data to Species level
tse <- mergeFeaturesByRank(tse,
                           rank = "Species")

# Set seed for reproducibility
set.seed(1576)
# We choose 99 random permutations. Consider applying more (999 or 9999) in your
# analysis. 
permanova <- adonis2(t(assay(tse, "relabundance")) ~ Group,
                     by = "margin", # each term (here only 'Group') analyzed individually
                     data = colData(tse),
                     method = "euclidean",
                     permutations = 99)

# Set seed for reproducibility
set.seed(1576)
# Perform dbRDA
dbrda <- dbrda(t(assay(tse,"relabundance")) ~ Group, 
               data = colData(tse))
# Perform permutational analysis
permanova2 <- anova.cca(dbrda,
                        by = "margin", # each term (here only 'Group') analyzed individually
                        method = "euclidean",
                        permutations = 99)

# Get p-values
p_values <- c(permanova["Group", "Pr(>F)"], permanova2["Group", "Pr(>F)"])
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values


## ------------------------------------------------------
# Add taxa info
sppscores(dbrda) <- t(assay(tse, "relabundance"))
# Get coefficients
coef <- dbrda$CCA$v
# Get the taxa with biggest weights
top.coef <- head(coef[rev(order(abs(coef))), , drop = FALSE], 20)
# Sort weights in increasing order
top.coef <- top.coef[order(top.coef), ]
# Get top names
top_names <- names(top.coef)[order(abs(top.coef), decreasing = TRUE)]


## ----plot-top-coef-anova, fig.cap = ""-----------------
df <- data.frame(x = top.coef,
                 y = factor(names(top.coef), unique(names(top.coef))))

ggplot(df, aes(x = x, y = y)) +
  geom_bar(stat = "identity") +
  labs(x = "", y= "", title = "Top Taxa") +
  theme_bw()


## ------------------------------------------------------
anova(betadisper(vegdist(t(assay(tse, "counts"))), colData(tse)$Group))

