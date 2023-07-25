## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----include=FALSE------------------------------------------------------------
# REMOVE THIS CHUNK ################################################################
# when bioc devel version has mergeSE 
if( !require(devtools) ){
    install.packages("devtools")
    library(devtools)
}
devtools::install_github("microbiome/mia", upgrade = "never")
library(mia)
####################################


## -----------------------------------------------------------------------------
library(mia)
data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
tse <- transformAssay(tse, MARGIN = "samples", method="relabundance")
molten_tse <- mia::meltAssay(tse,
                        add_row_data = TRUE,
                        add_col_data = TRUE,
                        assay.type = "relabundance")
molten_tse


## ----include = FALSE----------------------------------------------------------
# Load libraries and data
library(mia)
library(dplyr)
library(knitr)


## -----------------------------------------------------------------------------
# Store data into se and check dimensions
data("GlobalPatterns", package="mia")
tse <- GlobalPatterns
# Show dimensions (features x samples)
dim(tse) 


## -----------------------------------------------------------------------------
# Inspect possible values for SampleType
unique(tse$SampleType)

## ----eval = FALSE-------------------------------------------------------------
## # Show the frequency of each value
## tse$SampleType %>% table()

## ----echo = FALSE-------------------------------------------------------------
# Show the frequency of each value
tse$SampleType %>% table() %>% kable() %>%
    kableExtra::kable_styling("striped", latex_options="scale_down") %>% 
    kableExtra::scroll_box(width = "100%")


## -----------------------------------------------------------------------------
# Subset by sample
tse_subset_by_sample <- tse[ , tse$SampleType %in% c("Feces", "Skin", "Tongue")]

# Show dimensions
dim(tse_subset_by_sample)


## -----------------------------------------------------------------------------
# Inspect possible values for phylum
unique(rowData(tse)$Phylum)

## ----eval = FALSE-------------------------------------------------------------
## # Show the frequency of each value
## rowData(tse)$Phylum %>% table()

## ----echo = FALSE-------------------------------------------------------------
# Show te frequency of each value
rowData(tse)$Phylum %>% table() %>% kable() %>%
    kableExtra::kable_styling("striped", latex_options="scale_down") %>% 
    kableExtra::scroll_box(width = "100%")


## -----------------------------------------------------------------------------
# Subset by feature
tse_subset_by_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), ]

# Show dimensions
dim(tse_subset_by_feature)


## -----------------------------------------------------------------------------
# Agglomerate by phylum
tse_phylum <- tse %>% agglomerateByRank(rank = "Phylum")

# Subset by feature and remove NAs
tse_phylum_subset_by_feature <- tse_phylum[rowData(tse_phylum)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse_phylum)$Phylum), ]

# Show dimensions
dim(tse_phylum_subset_by_feature)


## -----------------------------------------------------------------------------
# Store features of interest into phyla
phyla <- c("Phylum:Actinobacteria", "Phylum:Chlamydiae")
# subset by feature
tse_phylum_subset_by_feature <- tse_phylum[phyla, ]
# Show dimensions
dim(tse_subset_by_feature)


## -----------------------------------------------------------------------------
# Subset by sample and feature and remove NAs
tse_subset_by_sample_feature <- tse[rowData(tse)$Phylum %in% c("Actinobacteria", "Chlamydiae") & !is.na(rowData(tse)$Phylum), tse$SampleType %in% c("Feces", "Skin", "Tongue")]

# Show dimensions
dim(tse_subset_by_sample_feature)


## -----------------------------------------------------------------------------
# Agglomerate data at Genus level 
tse_genus <- agglomerateByRank(tse, rank = "Genus")
# List bacteria that we want to include
genera <- c("Class:Thermoprotei", "Genus:Sulfolobus", "Genus:Sediminicola")
# Subset data
tse_genus_sub <- tse_genus[genera, ]

tse_genus_sub


## -----------------------------------------------------------------------------
# List total counts of each sample
colSums(assay(tse_genus_sub, "counts"))


## -----------------------------------------------------------------------------
# Remove samples that do not contain any bacteria
tse_genus_sub <- tse_genus_sub[ , colSums(assay(tse_genus_sub, "counts")) != 0 ]
tse_genus_sub


## -----------------------------------------------------------------------------
# Take only those samples that are collected from feces, skin, or tongue
tse_genus_sub <- tse_genus[ , tse_genus$SampleType %in% c("Feces", "Skin", "Tongue")]

tse_genus_sub


## -----------------------------------------------------------------------------
# What is the number of bacteria that are not present?
sum(rowSums(assay(tse_genus_sub, "counts")) == 0)


## -----------------------------------------------------------------------------
# Take only those bacteria that are present
tse_genus_sub <- tse_genus_sub[rowSums(assay(tse_genus_sub, "counts")) > 0, ]

tse_genus_sub


## ----splitbyRanks-------------------------------------------------------------
altExps(tse) <- splitByRanks(tse)
altExps(tse)


## ----splitOn------------------------------------------------------------------
splitOn(tse, "SampleType")


## ----modify-coldata-----------------------------------------------------------
# modify the Description entries
colData(tse)$Description <- paste(colData(tse)$Description, "modified description")

# view modified variable
head(tse$Description)


## ----add-coldata--------------------------------------------------------------
# simulate new data
new_data <- runif(ncol(tse))

# store new data as new variable in colData
colData(tse)$NewVariable <- new_data

# view new variable
head(tse$NewVariable)


## ----merge1-------------------------------------------------------------------
# Take subsets for demonstration purposes
tse1 <- tse[, 1]
tse2 <- tse[, 2]
tse3 <- tse[, 3]
tse4 <- tse[1:100, 4]


## ----merge2-------------------------------------------------------------------
# With inner join, we want to include all shared rows. When using mergeSEs function
# all samples are always preserved.
tse <- mergeSEs(list(tse1, tse2, tse3, tse4), join = "inner")
tse


## ----merge3-------------------------------------------------------------------
# Left join preserves all rows of the 1st object
tse <- mia::left_join(tse1, tse4, missing_values = 0)
tse

