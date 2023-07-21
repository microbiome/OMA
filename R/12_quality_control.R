## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ---- message=FALSE-----------------------------------------------------------
library(mia)


## ---- warning=FALSE, message=FALSE--------------------------------------------
# Load example data
library(miaTime)
library(miaViz)
data(hitchip1006)
tse <- hitchip1006

# Add relative abundances
tse <- transformCounts(tse, MARGIN = "samples", method = "relabundance")

# Use argument names
# assay.type / assay.type / assay.type
# depending on the mia package version
plotAbundanceDensity(tse, layout = "jitter", assay.type = "relabundance",
                     n = 40, point_size=1, point_shape=19, point_alpha=0.1) + 
                     scale_x_log10(label=scales::percent)


## ---- warning=FALSE, message=FALSE--------------------------------------------
plotAbundanceDensity(tse, layout = "density", assay.type = "relabundance",
                     n = 5, colour_by="nationality", point_alpha=1/10) +
    scale_x_log10()


## ----exploration-prevalence---------------------------------------------------
head(getPrevalence(tse, detection = 1/100, sort = TRUE, as_relative = TRUE))


## ----concepts_prevalence2-----------------------------------------------------
head(getPrevalence(tse, detection = 1, sort = TRUE, assay.type = "counts",
                   as_relative = FALSE))


## -----------------------------------------------------------------------------
# Agglomerate taxa abundances to Phylum level, and add the new table
# to the altExp slot
altExp(tse,"Phylum") <- agglomerateByRank(tse, "Phylum")
# Check prevalence for the Phylum abundance table from the altExp slot
head(getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = TRUE,
                   assay.type = "counts", as_relative = TRUE))


## -----------------------------------------------------------------------------
head(getPrevalence(tse, rank = "Phylum", detection = 1/100, sort = TRUE,
                   assay.type = "counts", as_relative = TRUE))


## ----core-members, message=FALSE, warning=FALSE, eval = FALSE-----------------
## getPrevalentFeatures(tse, detection = 0, prevalence = 50/100)
## prev <- getPrevalentFeatures(tse, detection = 0, prevalence = 50/100,
##                          rank = "Phylum", sort = TRUE)
## prev


## -----------------------------------------------------------------------------
rowData(altExp(tse,"Phylum"))$prevalence <- 
    getPrevalence(altExp(tse,"Phylum"), detection = 1/100, sort = FALSE,
                  assay.type = "counts", as_relative = TRUE)


## ---- message=FALSE, warning=FALSE--------------------------------------------
library(scater)
plotRowData(altExp(tse,"Phylum"), "prevalence", colour_by = "Phylum")


## -----------------------------------------------------------------------------
altExps(tse) <- splitByRanks(tse)
altExps(tse) <-
   lapply(altExps(tse),
          function(y){
              rowData(y)$prevalence <- 
                  getPrevalence(y, detection = 1/100, sort = FALSE,
                                assay.type = "counts", as_relative = TRUE)
              y
          })
top_phyla <- getTopFeatures(altExp(tse,"Phylum"),
                        method="prevalence",
                        top=5L,
                        assay.type="counts")
top_phyla_mean <- getTopFeatures(altExp(tse,"Phylum"),
                             method="mean",
                             top=5L,
                             assay.type="counts")
x <- unsplitByRanks(tse, ranks = taxonomyRanks(tse)[1:6])
x <- addTaxonomyTree(x)


## ----plot-prev-prev, message=FALSE, fig.cap="Prevalence of top phyla as judged by prevalence"----
library(miaViz)
plotRowTree(x[rowData(x)$Phylum %in% top_phyla,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")


## ----plot-prev-mean, message=FALSE, fig.cap="Prevalence of top phyla as judged by mean abundance"----
plotRowTree(x[rowData(x)$Phylum %in% top_phyla_mean,],
            edge_colour_by = "Phylum",
            tip_colour_by = "prevalence",
            node_colour_by = "prevalence")


## ---- message=FALSE-----------------------------------------------------------
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns 


## ----top-feature-taxo---------------------------------------------------------
# Pick the top taxa
top_features <- getTopFeatures(tse, method="median", top=10)

# Check the information for these
rowData(tse)[top_features, taxonomyRanks(tse)]


## ----lib-size-----------------------------------------------------------------
library(scater)
perCellQCMetrics(tse)
tse <- addPerCellQC(tse)
colData(tse)


## ----plot-viz-lib-size-1, fig.width=8, fig.height=4, fig.cap="Library size distribution."----
library(ggplot2)

p1 <- ggplot(as.data.frame(colData(tse))) +
        geom_histogram(aes(x = sum), color = "black", fill = "gray", bins = 30) +
        labs(x = "Library size", y = "Frequency (n)") + 
        # scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
        # labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), # Removes the grid
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) # Adds y-axis

library(dplyr)
df <- as.data.frame(colData(tse)) %>%
        arrange(sum) %>%
        mutate(index = 1:n())
p2 <- ggplot(df, aes(y = index, x = sum/1e6)) +
        geom_point() +	
        labs(x = "Library size (million reads)", y = "Sample index") +	
        theme_bw() +
        theme(panel.grid.major = element_blank(), # Removes the grid
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) # Adds y-axis

library(patchwork)
p1 + p2


## ----plot-viz-lib-size-2, fig.width=8, fig.height=4, fig.cap="Library sizes per sample."----
library(ggplot2)
# Sort samples by read count, order the factor levels, and store back to tse as DataFrame
# TODO: plotColData could include an option for sorting samples based on colData variables
colData(tse) <- as.data.frame(colData(tse)) %>%
                 arrange(X.SampleID) %>%
        	 mutate(X.SampleID = factor(X.SampleID, levels=X.SampleID)) %>%
		 DataFrame
plotColData(tse,"sum","X.SampleID", colour_by = "SampleType") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    labs(y = "Library size (N)", x = "Sample ID") 	    


## ----plot-viz-lib-size-3, fig.width=8, fig.height=4, fig.cap="Library sizes per sample type."----
plotColData(tse,"sum","SampleType", colour_by = "SampleType") + 
    theme(axis.text.x = element_text(angle = 45, hjust=1))


## ----sessionInfo, echo=FALSE, results='asis'----------------------------------
prettySessionInfo()

