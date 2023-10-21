## ----setup, echo=FALSE, results="asis"---------------------------------------------------------------------------------
library(rebook)
chapterPreamble()


## ----load-pkg-data-----------------------------------------------------------------------------------------------------
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns


## ----------------------------------------------------------------------------------------------------------------------
tse <- mia::estimateRichness(tse, 
                             assay.type = "counts", 
                             index = "observed", 
                             name="observed")

head(tse$observed)


## ----plot-div-shannon, message=FALSE, fig.cap="Shannon diversity estimates plotted grouped by sample type with colour-labeled barcode."----
library(scater)
plotColData(tse, 
            "observed", 
            "SampleType", 
            colour_by = "Final_Barcode") +
    theme(axis.text.x = element_text(angle=45,hjust=1)) + 
  ylab(expression(Richness[Observed]))



## ----estimate-shannon--------------------------------------------------------------------------------------------------
tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index = "shannon", 
                              name = "shannon")
head(tse$shannon)


## ----visualize-shannon-------------------------------------------------------------------------------------------------
library(ggsignif)
library(ggplot2)
library(patchwork)
library(ggsignif)

# Subsets the data. Takes only those samples that are from feces, skin, or tongue,
# and creates data frame from the collected data
df <- as.data.frame(colData(tse)[tse$SampleType %in% 
                 c("Feces", "Skin", "Tongue"), ])

# Changes old levels with new levels
df$SampleType <- factor(df$SampleType)

# For significance testing, all different combinations are determined
comb <- split(t(combn(levels(df$SampleType), 2)), 
           seq(nrow(t(combn(levels(df$SampleType), 2)))))

ggplot(df, aes(x = SampleType, y = shannon)) +
  # Outliers are removed, because otherwise each data point would be plotted twice; 
  # as an outlier of boxplot and as a point of dotplot.
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  geom_signif(comparisons = comb, map_signif_level = FALSE) +
  theme(text = element_text(size = 10))


## ----phylo-div-1-------------------------------------------------------------------------------------------------------
tse <- mia::estimateFaith(tse,
                          assay.type = "counts")
head(tse$faith)


## ----phylo-div-2-------------------------------------------------------------------------------------------------------
plots <- lapply(c("shannon", "faith"),
                plotColData,
                object = tse, colour_by = "SampleType")
plots[[1]] + plots[[2]] +
  plot_layout(guides = "collect")


## ----phylo-div-3-------------------------------------------------------------------------------------------------------
tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index = "faith", 
                              name = "faith")


## ----evenness-1--------------------------------------------------------------------------------------------------------
tse <- estimateEvenness(tse, 
                        assay.type = "counts", 
                        index="simpson")
head(tse$simpson)


## ----dominance-1-------------------------------------------------------------------------------------------------------
tse <- estimateDominance(tse, 
                         assay.type = "counts", 
                         index="relative")

head(tse$relative)


## ----rarity-1----------------------------------------------------------------------------------------------------------
tse <- mia::estimateDiversity(tse, 
                              assay.type = "counts",
                              index = "log_modulo_skewness")

head(tse$log_modulo_skewness)


## ----------------------------------------------------------------------------------------------------------------------
tse <- mia::estimateDivergence(tse,
                               assay.type = "counts",
                               reference = "median",
                               FUN = vegan::vegdist)


## ----plot-all-diversities, fig.width = 6.5-----------------------------------------------------------------------------
plots <- lapply(c("observed", "shannon", "simpson", "relative", "faith", "log_modulo_skewness"),
                plotColData,
                object = tse,
                x = "SampleType",
                colour_by = "SampleType")

plots <- lapply(plots, "+", 
                theme(axis.text.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.ticks.x = element_blank()))

((plots[[1]] | plots[[2]] | plots[[3]]) / 
(plots[[4]] | plots[[5]] | plots[[6]])) +
  plot_layout(guides = "collect")

