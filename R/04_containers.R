## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----echo=FALSE---------------------------------------------------------------
knitr::include_graphics("general/figures/FigureOverviewV2_mod.png")


## -----------------------------------------------------------------------------
library(mia)
data(GlobalPatterns, package="mia")
tse <- GlobalPatterns
assays(tse)


## -----------------------------------------------------------------------------
assays(tse)


## -----------------------------------------------------------------------------
assay(tse, "counts")[1:5,1:7]


## -----------------------------------------------------------------------------
tse <- relAbundanceCounts(tse)
assays(tse)


## -----------------------------------------------------------------------------
assay(tse, "relabundance")[1:5,1:7]


## ----coldata------------------------------------------------------------------
colData(tse)


## ----rowdata------------------------------------------------------------------
rowData(tse)


## ----rowtree------------------------------------------------------------------
rowTree(tse)


## ----rowlinks-----------------------------------------------------------------
rowLinks(tse)


## -----------------------------------------------------------------------------
# Agglomerate the data to Phylym level
tse_phylum <- agglomerateByRank(tse, "Phylum")
# both have the same number of columns (samples)
dim(tse)
dim(tse_phylum)

# Add the new table as an alternative experiment
altExp(tse, "Phylum") <- tse_phylum
altExpNames(tse)

# Pick a sample subset: this acts on both altExp and assay data
tse[,1:10]
dim(altExp(tse[,1:10],"Phylum"))


## -----------------------------------------------------------------------------
#TODO: Find the right dataset to explain a non 1:1 sample relationship


## ----dada2_1, include=FALSE---------------------------------------------------
# Load objects
seqtab.nochim <- readRDS("data/dada2_seqtab.nochim")
taxa <- readRDS("data/dada2_taxa")


## ----dada2_2------------------------------------------------------------------
library(mia)
library(ggplot2)

if( !require("BiocManager") ){
    install.packages("BiocManager")
    library("BiocManager")
}

if( !require("Biostrings") ){
    BiocManager::install("Biostrings")
    library("Biostrings")
}
library(Biostrings)


## ----dada2_3------------------------------------------------------------------
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


## ----dada2_4------------------------------------------------------------------
# Create a list that contains assays
counts <- t(seqtab.nochim)
counts <- as.matrix(counts)
assays <- SimpleList(counts = counts)

# Convert colData and rowData into DataFrame
samdf <- DataFrame(samdf)
taxa <- DataFrame(taxa)

# Create TreeSE
tse <- TreeSummarizedExperiment(assays = assays,
                                colData = samdf,
                                rowData = taxa
                                )

# Remove mock sample like it is also done in DADA2 pipeline tutorial
tse <- tse[ , colnames(tse) != "mock"]


## ----dada2_5------------------------------------------------------------------
# Convert sequences into right format
dna <- Biostrings::DNAStringSet( rownames(tse) )
# Add sequences into referenceSeq slot
referenceSeq(tse) <- dna
# Convert rownames into ASV_number format
rownames(tse) <- paste0("ASV", seq( nrow(tse) ))
tse


## ----importingcsv1, message=FALSE---------------------------------------------
count_file <- "data/assay_taxa.csv"
tax_file <- "data/rowdata_taxa.csv"
sample_file <- "data/coldata.csv"

# Load files
counts  <- read.csv(count_file)   # Abundance table (e.g. ASV data; to assay data)
tax     <- read.csv(tax_file)     # Taxonomy table (to rowData)
samples <- read.csv(sample_file)  # Sample data (to colData)


## ----importingcsv2------------------------------------------------------------
# Add rownames and remove an additional column
rownames(counts) <- counts$X
counts$X <- NULL

# Add rownames and remove an additional column
rownames(samples) <- samples$X
samples$X <- NULL

# Add rownames and remove an additional column
rownames(tax) <- tax$X
tax$X <- NULL

# As an example:
# If e.g. samples do not match between colData and counts table, you must order 
# counts based on colData
if( any( colnames(counts) != rownames(samples) ) ){
    counts <- counts[ , rownames(samples) ]
}

# And same with rowData and counts...
if( any( rownames(counts) != rownames(tax) ) ){
    counts <- counts[ rownames(tax), ]
}


## ----importingcsv3------------------------------------------------------------
# Ensure that the data is in correct format

# counts should be in matrix format
counts <- as.matrix(counts)
# And it should be added to a SimpleList
assays <-  SimpleList(counts = counts)

# colData and rowData should be in DataFrame format
colData <- DataFrame(colData)
rowData <- DataFrame(rowData)

# Create a TreeSE
tse_taxa <- TreeSummarizedExperiment(assays = assays,
                                     colData = samples,
                                     rowData = tax)

tse_taxa


## ----importingcsv4, message=FALSE---------------------------------------------
count_file <- "data/assay_metabolites.csv"
sample_file <- "data/coldata.csv"

# Load files
counts  <- read.csv(count_file)  
samples <- read.csv(sample_file)

# Add rownames and remove an additional column
rownames(counts) <- counts$X
counts$X <- NULL
rownames(samples) <- samples$X
samples$X <- NULL

# Convert into right format
counts <- as.matrix(counts)
assays <-  SimpleList(concs = counts)
colData <- DataFrame(colData)

# Create a TreeSE
tse_metabolite <- TreeSummarizedExperiment(assays = assays,
                                           colData = samples)
tse_metabolite


## ----importingcsv5------------------------------------------------------------
# Create an ExperimentList that includes experiments
experiments <- ExperimentList(microbiome = tse_taxa, 
                              metabolite = tse_metabolite)

# Create a MAE
mae <- MultiAssayExperiment(experiments = experiments)

mae


## -----------------------------------------------------------------------------
biom_file_path <- "data/Aggregated_humanization2.biom"
sample_meta_file_path <- "data/Mapping_file_ADHD_aggregated.csv"
tree_file_path <- "data/Data_humanization_phylo_aggregation.tre"


## -----------------------------------------------------------------------------
library(mia)

# Imports the data
se <- loadFromBiom(biom_file_path)

# Check
se


## -----------------------------------------------------------------------------
assays(se)$counts[1:3, 1:3]


## -----------------------------------------------------------------------------
head(rowData(se))


## -----------------------------------------------------------------------------
names(rowData(se)) <- c("Kingdom", "Phylum", "Class", "Order", 
                        "Family", "Genus")

# Goes through the whole DataFrame. Removes '.*[kpcofg]__' from strings, where [kpcofg] 
# is any character from listed ones, and .* any character.
rowdata_modified <- BiocParallel::bplapply(rowData(se), 
                                           FUN = stringr::str_remove, 
                                           pattern = '.*[kpcofg]__')

# Genus level has additional '\"', so let's delete that also
rowdata_modified <- BiocParallel::bplapply(rowdata_modified, 
                                           FUN = stringr::str_remove, 
                                           pattern = '\"')

# rowdata_modified is a list, so it is converted back to DataFrame format. 
rowdata_modified <- DataFrame(rowdata_modified)

# And then assigned back to the SE object
rowData(se) <- rowdata_modified

# Now we have a nicer table
head(rowData(se))


## -----------------------------------------------------------------------------
head(colData(se))


## -----------------------------------------------------------------------------
# We use this to check what type of data it is
# read.table(sample_meta_file_path)

# It seems like a comma separated file and it does not include headers
# Let us read it and then convert from data.frame to DataFrame
# (required for our purposes)
sample_meta <- DataFrame(read.table(sample_meta_file_path, sep = ",", header = FALSE))

# Add sample names to rownames
rownames(sample_meta) <- sample_meta[,1]

# Delete column that included sample names
sample_meta[,1] <- NULL

# We can add headers
colnames(sample_meta) <- c("patient_status", "cohort", "patient_status_vs_cohort", "sample_name")

# Then it can be added to colData
colData(se) <- sample_meta


## -----------------------------------------------------------------------------
head(colData(se))


## -----------------------------------------------------------------------------
tse <- as(se, "TreeSummarizedExperiment")

# tse includes same data as se
tse


## -----------------------------------------------------------------------------
# Reads the tree file
tree <- ape::read.tree(tree_file_path)

# Add tree to rowTree
rowTree(tse) <- tree

# Check
tse


## ---- eval=FALSE--------------------------------------------------------------
## head(rowTree(tse))


## ---- message=FALSE-----------------------------------------------------------
library(mia)

# phyloseq example data
data(GlobalPatterns, package="phyloseq") 
GlobalPatterns_phyloseq <- GlobalPatterns
GlobalPatterns_phyloseq


## ---- message=FALSE-----------------------------------------------------------
# convert phyloseq to TSE
GlobalPatterns_TSE <- makeTreeSummarizedExperimentFromPhyloseq(GlobalPatterns_phyloseq) 
GlobalPatterns_TSE


## ---- message=FALSE-----------------------------------------------------------
# convert TSE to phyloseq
GlobalPatterns_phyloseq2 <- makePhyloseqFromTreeSummarizedExperiment(GlobalPatterns_TSE) 
GlobalPatterns_phyloseq2


## ---- message=FALSE, eval=FALSE-----------------------------------------------
## library(mia)
## data(package="mia")


## ---- message=FALSE-----------------------------------------------------------
data("GlobalPatterns", package="mia")
GlobalPatterns


## ---- message=FALSE, echo=FALSE-----------------------------------------------
help(GlobalPatterns)


## ---- message=FALSE-----------------------------------------------------------
library(microbiomeDataSets)
availableDataSets()


## ----eval=FALSE, message=FALSE------------------------------------------------
## # mae <- HintikkaXOData()
## # Since HintikkaXOData is now added to mia, we can load it directly from there
## # We suggest to check other datasets from microbiomeDataSets
## data(HintikkaXOData)
## mae <- HintikkaXOData


## ---- message=FALSE, eval=FALSE-----------------------------------------------
## library(curatedMetagenomicData)
## tse <- curatedMetagenomicData("Vatanen*", dryrun = FALSE, counts = TRUE)


## ----sessionInfo, echo=FALSE, results='asis'----------------------------------
prettySessionInfo()


