## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ---- message=FALSE-----------------------------------------------------------
library(mia)
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns 


## -----------------------------------------------------------------------------
checkTaxonomy(tse)


## -----------------------------------------------------------------------------
taxonomyRanks(tse)


## -----------------------------------------------------------------------------
rowData(tse)[,taxonomyRanks(tse)]


## -----------------------------------------------------------------------------
all(!taxonomyRankEmpty(tse, rank = "Kingdom"))
table(taxonomyRankEmpty(tse, rank = "Genus"))
table(taxonomyRankEmpty(tse, rank = "Species"))


## -----------------------------------------------------------------------------
head(getTaxonomyLabels(tse))


## -----------------------------------------------------------------------------
phylum <- !is.na(rowData(tse)$Phylum) & 
    vapply(data.frame(apply(rowData(tse)[,taxonomyRanks(tse)[3:7]],1L,is.na)),all,logical(1))
head(getTaxonomyLabels(tse[phylum,]))
head(getTaxonomyLabels(tse[phylum,], with_rank = TRUE))


## -----------------------------------------------------------------------------
head(getTaxonomyLabels(tse[phylum,], with_rank = TRUE, make_unique = FALSE))


## -----------------------------------------------------------------------------
head(getUniqueFeatures(tse, rank = "Phylum"))


## -----------------------------------------------------------------------------
taxonomyTree(tse)


## -----------------------------------------------------------------------------
tse <- addTaxonomyTree(tse)
tse


## -----------------------------------------------------------------------------
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
altExp(tse, "Family") <- mergeFeaturesByRank(tse, rank = "Family",
                                           agglomerateTree = TRUE)
altExp(tse, "Family")


## -----------------------------------------------------------------------------
assayNames(tse)
assayNames(altExp(tse, "Family"))


## -----------------------------------------------------------------------------
assay(altExp(tse, "Family"), "relabundance")[1:5,1:7]


## ----taxinfo_altexp_example---------------------------------------------------
assay(altExp(tse, "Family"), "counts")[1:5,1:7]


## -----------------------------------------------------------------------------
assay(tse, "pseudo") <- assay(tse, "counts") + 1
tse <- transformAssay(tse, assay.type = "pseudo", method = "relabundance")
tse <- transformAssay(x = tse, assay.type = "relabundance", method = "clr", 
                        pseudocount = 1, name = "clr_transformation")

head(assay(tse, "clr_transformation"))


## -----------------------------------------------------------------------------
tse <- transformAssay(tse, method = "pa")

head(assay(tse, "pa"))


## -----------------------------------------------------------------------------
# list of abundance tables that assays slot contains
assays(tse)


## -----------------------------------------------------------------------------
taxa.abund.cc1 <- getAbundanceSample(tse, 
                                     sample_id = "CC1",
                                     assay.type = "counts")
taxa.abund.cc1[1:10]


## -----------------------------------------------------------------------------
taxa.abundances <- getAbundanceFeature(tse, 
                                      feature_id = "Phylum:Bacteroidetes",
                                      assay.type = "counts")
taxa.abundances[1:10]


## ----sessionInfo, echo=FALSE, results='asis'----------------------------------
prettySessionInfo()

