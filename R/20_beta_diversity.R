## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----include=FALSE------------------------------------------------------------
# global knitr options
knitr::opts_chunk$set(
  fig.width=10,
  dpi=300,
  dev = "png",
  dev.args = list(type = "cairo-png")
)


## ----relstress----------------------------------------------------------------
# Example data
library(mia)
data(GlobalPatterns, package="mia")

# Data matrix (features x samples)
tse <- GlobalPatterns
tse <- transformAssay(tse, method = "relabundance")

# Add group information Feces yes/no
colData(tse)$Group <- tse$SampleType=="Feces"

# Quantify dissimilarities in the original feature space
library(vegan)
x <- assay(tse, "relabundance") # Pick relabunance assay separately
d0 <- as.matrix(vegdist(t(x), "bray"))

# PCoA Ordination
pcoa <- as.data.frame(cmdscale(d0, k = 2))
names(pcoa) <- c("PCoA1", "PCoA2")

# Quantify dissimilarities in the ordination space
dp <- as.matrix(dist(pcoa))

# Calculate stress i.e. relative difference in the original and
# projected dissimilarities
stress <- sum((dp - d0)^2)/sum(d0^2)


## ----shepard------------------------------------------------------------------
ord <- order(as.vector(d0))
df <- data.frame(d0 = as.vector(d0)[ord],
                  dmds = as.vector(dp)[ord])

library(ggplot2)
ggplot(aes(x = d0, y = dmds), data=df) + 
       geom_smooth() +
       geom_point() +       
       labs(title = "Shepard plot",
       x = "Original distance",
       y = "MDS distance",       
            subtitle = paste("Stress:", round(stress, 2))) +
  theme_bw()


## ----runMDS, message=FALSE----------------------------------------------------
library(scater)

# Bray-Curtis is usually applied to relative abundances
tse <- transformAssay(tse, method = "relabundance")
# Perform PCoA
tse <- runMDS(tse, FUN = vegan::vegdist, method = "bray", name = "PCoA_BC", assay.type = "relabundance")


## ----plot-mds-bray-curtis, fig.cap="MDS plot based on the Bray-Curtis distances on the GlobalPattern dataset."----
# Create ggplot object
p <- plotReducedDim(tse, "PCoA_BC", colour_by = "Group")

# Add explained variance for each axis
e <- attr(reducedDim(tse, "PCoA_BC"), "eig");
rel_eig <- e/sum(e[e>0])		  
p <- p + labs(x = paste("PCoA 1 (", round(100 * rel_eig[[1]],1), "%", ")", sep = ""),
              y = paste("PCoA 2 (", round(100 * rel_eig[[2]],1), "%", ")", sep = ""))

print(p)


## ----plot-mds-nmds-comparison, fig.cap="Comparison of MDS and NMDS plots based on the Bray-Curtis or euclidean distances on the GlobalPattern dataset.", message=FALSE----
tse <- runMDS(tse, FUN = vegan::vegdist, name = "MDS_euclidean",
             method = "euclidean", assay.type = "counts")
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "NMDS_BC")
tse <- runNMDS(tse, FUN = vegan::vegdist, name = "NMDS_euclidean",
               method = "euclidean")
plots <- lapply(c("PCoA_BC", "MDS_euclidean", "NMDS_BC", "NMDS_euclidean"),
                plotReducedDim,
                object = tse,
                colour_by = "Group")

library(patchwork)
plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] +
  plot_layout(guides = "collect")


## -----------------------------------------------------------------------------
library(scater)
tse <- runMDS(tse, FUN = mia::calculateUnifrac, name = "Unifrac",
              tree = rowTree(tse),
              ntop = nrow(tse),
             assay.type = "counts")


## ----plot-unifrac, fig.cap="Unifrac distances scaled by MDS of the GlobalPattern dataset."----
plotReducedDim(tse, "Unifrac", colour_by = "Group")


## -----------------------------------------------------------------------------
tse <- runPCA(tse, name = "PCA", assay.type = "counts", ncomponents = 10)


## ----plot-pca, fig.cap="PCA plot on the GlobalPatterns data set containing sample from different sources."----
plotReducedDim(tse, "PCA", colour_by = "Group")


## -----------------------------------------------------------------------------
tse <- runTSNE(tse, name = "TSNE", assay.type = "counts", ncomponents = 3)


## ----plot-tsne, fig.cap="t-SNE plot on the GlobalPatterns data set containing sample from different sources."----
plotReducedDim(tse, "TSNE", colour_by = "Group", ncomponents = c(1:3))


## ----microbiome_RDA1----------------------------------------------------------
# Load required packages
if(!require("vegan")){
    install.packages("vegan")
    library("vegan")
}
if(!require("stringr")){
    install.packages("stringr")
    library("stringr")
}
if(!require("knitr")){
    install.packages("knitr")
    library("knitr")
}
# Load data
data(enterotype, package="mia")
# Covariates that are being analyzed
variable_names <- c("ClinicalStatus", "Gender", "Age")

# Apply relative transform
enterotype <- transformAssay(enterotype, method = "relabundance")

# Create a formula
formula <- as.formula(paste0("assay ~ ", str_c(variable_names, collapse = " + ")) )

# # Perform RDA
rda <- calculateRDA(enterotype, assay.type = "relabundance",
                    formula = formula, distance = "bray", na.action = na.exclude)
# Get the rda object
rda <- attr(rda, "rda")
# Calculate p-value and variance for whole model
# Recommendation: use 999 permutations instead of 99
set.seed(436)
permanova <- anova.cca(rda, permutations = 99)
# Create a data.frame for results
rda_info <- as.data.frame(permanova)["Model", ]

# Calculate p-value and variance for each variable
# by = "margin" --> the order or variables does not matter
set.seed(4585)
permanova <- anova.cca(rda, by = "margin",  permutations = 99)
# Add results to data.frame
rda_info <- rbind(rda_info, permanova)

# Add info about total variance
rda_info[ , "Total variance"] <- rda_info["Model", 2] +
    rda_info["Residual", 2]

# Add info about explained variance
rda_info[ , "Explained variance"] <- rda_info[ , 2] / 
    rda_info[ , "Total variance"]

# Loop through variables, calculate homogeneity
homogeneity <- list()
# Get colDtaa
coldata <- colData(enterotype)
# Get assay
assay <- t(assay(enterotype, "relabundance"))
for( variable_name in rownames(rda_info) ){
    # If data is continuous or discrete
    if( variable_name %in% c("Model", "Residual") ||
        length(unique(coldata[[variable_name]])) /
        length(coldata[[variable_name]]) > 0.2 ){
        # Do not calculate homogeneity for continuous data
        temp <- NA
    } else{
        # Calculate homogeneity for discrete data
        # Calculate homogeneity
        set.seed(413)
        temp <- anova(
            betadisper( 
                vegdist(assay, method = "bray"),
                group = coldata[[variable_name]] ),
            permutations = permutations )["Groups", "Pr(>F)"]
    }
    # Add info to the list
    homogeneity[[variable_name]] <- temp
}
# Add homogeneity to information
rda_info[["Homogeneity p-value (NULL hyp: distinct/homogeneous --> permanova suitable)"]] <-
    homogeneity

kable(rda_info)


## ----microbiome_RDA2----------------------------------------------------------
# Load ggord for plotting
if(!require("ggord")){
    if(!require("devtools")){
        install.packages("devtools")
        library("devtools")
    }
    install_github("https://github.com/fawda123/ggord/")
    library("ggord")
}
if(!require("ggplot2")){
    install.packages("ggplot2")
    library("ggplot2")
}
# Since na.exclude was used, if there were rows missing information, they were 
# dropped off. Subset coldata so that it matches with rda.
coldata <- coldata[ rownames(rda$CCA$wa), ]

# Adjust names
# Get labels of vectors
vec_lab_old <- rownames(rda$CCA$biplot)

# Loop through vector labels
vec_lab <- sapply(vec_lab_old, FUN = function(name){
    # Get the variable name
    variable_name <- variable_names[ str_detect(name, variable_names) ]
    # If the vector label includes also group name
    if( !any(name %in% variable_names) ){
        # Get the group names
        group_name <- unique( coldata[[variable_name]] )[ 
        which( paste0(variable_name, unique( coldata[[variable_name]] )) == name ) ]
        # Modify vector so that group is separated from variable name
        new_name <- paste0(variable_name, " \U2012 ", group_name)
    } else{
        new_name <- name
    }
    # Add percentage how much this variable explains, and p-value
    new_name <- expr(paste(!!new_name, " (", 
                           !!format(round( rda_info[variable_name, "Explained variance"]*100, 1), nsmall = 1), 
                           "%, ",italic("P"), " = ", 
                           !!gsub("0\\.","\\.", format(round( rda_info[variable_name, "Pr(>F)"], 3), 
                                                       nsmall = 3)), ")"))

    return(new_name)
})
# Add names
names(vec_lab) <- vec_lab_old

# Create labels for axis
xlab <- paste0("RDA1 (", format(round( rda$CCA$eig[[1]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")
ylab <- paste0("RDA2 (", format(round( rda$CCA$eig[[2]]/rda$CCA$tot.chi*100, 1), nsmall = 1 ), "%)")

# Create a plot        
plot <- ggord(rda, grp_in = coldata[["ClinicalStatus"]], vec_lab = vec_lab,
              alpha = 0.5,
              size = 4, addsize = -4,
              #ext= 0.7, 
              txt = 3.5, repel = TRUE, 
              #coord_fix = FALSE
          ) + 
    # Adjust titles and labels
    guides(colour = guide_legend("ClinicalStatus"),
           fill = guide_legend("ClinicalStatus"),
           group = guide_legend("ClinicalStatus"),
           shape = guide_legend("ClinicalStatus"),
           x = guide_axis(xlab),
           y = guide_axis(ylab)) +
    theme( axis.title = element_text(size = 10) )
plot


## -----------------------------------------------------------------------------
# Agglomerate to genus level
tse_Genus <- mergeFeaturesByRank(tse, rank="Genus")
# Convert to relative abundances
tse_Genus <- transformAssay(tse, method = "relabundance", assay.type="counts")
# Add info on dominant genus per sample
tse_Genus <- addPerSampleDominantFeatures(tse_Genus, assay.type="relabundance", name = "dominant_taxa")


## -----------------------------------------------------------------------------
tse_Genus <- runMDS(tse_Genus, FUN = vegan::vegdist,
              name = "PCoA_BC", assay.type = "relabundance")


## -----------------------------------------------------------------------------
# Getting the top taxa
top_taxa <- getTopFeatures(tse_Genus,top = 6, assay.type = "relabundance")

# Naming all the rest of non top-taxa as "Other"
most_abundant <- lapply(colData(tse_Genus)$dominant_taxa,
                   function(x){if (x %in% top_taxa) {x} else {"Other"}})

# Storing the previous results as a new column within colData
colData(tse_Genus)$most_abundant <- as.character(most_abundant)

# Calculating percentage of the most abundant
most_abundant_freq <- table(as.character(most_abundant))
most_abundant_percent <- round(most_abundant_freq/sum(most_abundant_freq)*100, 1)

# Retrieving the explained variance
e <- attr(reducedDim(tse_Genus, "PCoA_BC"), "eig");
var_explained <- e/sum(e[e>0])*100

# Visualization
plot <-plotReducedDim(tse_Genus,"PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(most_abundant_percent),"(",most_abundant_percent,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       color="")
plot


## -----------------------------------------------------------------------------
# Calculating the frequencies and percentages for both categories
freq_TRUE <- table(as.character(most_abundant[colData(tse_Genus)$Group==TRUE]))
freq_FALSE <- table(as.character(most_abundant[colData(tse_Genus)$Group==FALSE]))
percent_TRUE <- round(freq_TRUE/sum(freq_TRUE)*100, 1)
percent_FALSE <- round(freq_FALSE/sum(freq_FALSE)*100, 1)

# Visualization
plotReducedDim(tse_Genus[,colData(tse_Genus)$Group==TRUE],
                          "PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(percent_TRUE),"(",percent_TRUE,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       title = "Group = TRUE", color="")

plotReducedDim(tse_Genus[,colData(tse_Genus)$Group==FALSE],
                          "PCoA_BC", colour_by = "most_abundant") +
  scale_colour_manual(values = c("black", "blue", "lightblue", "darkgray", "magenta", "darkgreen", "red"),
                      labels=paste0(names(percent_FALSE),"(",percent_FALSE,"%)"))+
  labs(x=paste("PC 1 (",round(var_explained[1],1),"%)"),
       y=paste("PC 2 (",round(var_explained[2],1),"%)"),
       title = "Group = FALSE", color="")


## -----------------------------------------------------------------------------
if( !require(vegan) ){
    BiocManager::install("vegan")
    library("vegan")
}
# Agglomerate data to Species level
tse <- mergeFeaturesByRank(tse, rank = "Species")

# Set seed for reproducibility
set.seed(1576)
# We choose 99 random permutations. Consider applying more (999 or 9999) in your
# analysis. 
permanova <- adonis2(t(assay(tse,"relabundance")) ~ Group,
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
p_values <- c( permanova["Group", "Pr(>F)"], permanova2["Group", "Pr(>F)"] )
p_values <-as.data.frame(p_values)
rownames(p_values) <- c("adonis2", "dbRDA+anova.cca")
p_values


## -----------------------------------------------------------------------------
# Add taxa info
sppscores(dbrda) <- t(assay(tse,"relabundance"))
# Get coefficients
coef <- dbrda$CCA$v
# Get the taxa with biggest weights
top.coef <- head( coef[rev(order(abs(coef))), , drop = FALSE], 20)
# Sort weights in increasing order
top.coef <- top.coef[ order(top.coef), ]
# Get top names
top_names <- names(top.coef)[ order(abs(top.coef), decreasing = TRUE) ]


## ----plot-top-coef-anova, fig.cap=""------------------------------------------
ggplot(data.frame(x = top.coef,
                  y = factor(names(top.coef),
                                      unique(names(top.coef)))),
        aes(x = x, y = y)) +
    geom_bar(stat="identity") +
    labs(x="",y="",title="Top Taxa") +
    theme_bw()


## -----------------------------------------------------------------------------
anova( betadisper(vegdist(t(assay(tse, "counts"))), colData(tse)$Group) )


## ----sessionInfo, echo=FALSE, results='asis'----------------------------------
prettySessionInfo()

