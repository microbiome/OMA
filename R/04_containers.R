## ----setup, echo=FALSE, results="asis"-----------------
library(rebook)
chapterPreamble()


## ----echo=FALSE----------------------------------------
knitr::include_graphics("general/figures/FigureOverviewV2_mod.png")


## ------------------------------------------------------
library(mia)
data("hitchip1006", package = "miaTime")
tse <- hitchip1006


## ------------------------------------------------------
assays(tse)


## ------------------------------------------------------
assay(tse, "counts")[1:5,1:7]


## ------------------------------------------------------
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
assays(tse)


## ------------------------------------------------------
assay(tse, "relabundance")[1:5,1:7]


## ----coldata-------------------------------------------
colData(tse)


## ----rowdata-------------------------------------------
rowData(tse)


## ----rowtree-------------------------------------------
rowTree(tse)


## ----rowlinks------------------------------------------
rowLinks(tse)


## ----altexp_agglomerate--------------------------------
tse_phylum <- mergeFeaturesByRank(tse, "Phylum", na.rm=TRUE)
# Both have the same number of columns (samples)
dim(tse)
dim(tse_phylum)


## ----altexp_agglomerate2-------------------------------
# Add the new data object to the original data object as an alternative experiment with the name "Phylum"
altExp(tse, "Phylum") <- tse_phylum

# Check the alternative experiment names available in the data
altExpNames(tse)


## ----altexp_agglomerate3-------------------------------
tse[,1:10]
dim(altExp(tse[,1:10],"Phylum"))


## ------------------------------------------------------
#TODO: Find the right dataset to explain a non 1:1 sample relationship


## ---- message=FALSE, eval=FALSE------------------------
## library(mia)
## data(package="mia")


## ---- message=FALSE------------------------------------
data("GlobalPatterns", package="mia")
GlobalPatterns


## ---- message=FALSE------------------------------------
library(microbiomeDataSets)
availableDataSets()


## ----eval=FALSE, message=FALSE-------------------------
## # mae <- HintikkaXOData()
## # Since HintikkaXOData is now added to mia, we can load it directly from there
## # We suggest to check other datasets from microbiomeDataSets
## data(HintikkaXOData, package = "mia")
## mae <- HintikkaXOData


## ---- message=FALSE, eval=FALSE------------------------
## library(curatedMetagenomicData)
## tse <- curatedMetagenomicData("Vatanen*", dryrun = FALSE, counts = TRUE)


## ---- message=FALSE, eval=FALSE------------------------
## library(MicroBioMap)
## cpd <- getCompendium()


## ----dada2_1, include=FALSE----------------------------
# Load objects
seqtab.nochim <- readRDS("data/dada2_seqtab.nochim")
taxa <- readRDS("data/dada2_taxa")


## ----dada2_2-------------------------------------------
library(mia)
library(ggplot2)
library(BiocManager)
library(Biostrings)


## ----dada2_3-------------------------------------------
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out


## ----dada2_4-------------------------------------------
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


## ----dada2_5-------------------------------------------
# Convert sequences into right format
dna <- Biostrings::DNAStringSet( rownames(tse) )
# Add sequences into referenceSeq slot
referenceSeq(tse) <- dna
# Convert rownames into ASV_number format
rownames(tse) <- paste0("ASV", seq( nrow(tse) ))
tse


## ----importingcsv1, message=FALSE----------------------
count_file  <- "data/assay_taxa.csv"
tax_file    <- "data/rowdata_taxa.csv"
sample_file <- "data/coldata.csv"

# Load files
counts  <- read.csv(count_file, row.names=1)   # Abundance table (e.g. ASV data; to assay data)
tax     <- read.csv(tax_file, row.names=1)     # Taxonomy table (to rowData)
samples <- read.csv(sample_file, row.names=1)  # Sample data (to colData)


## ----importingcsv2-------------------------------------
# Match rows and columns
counts <- counts[rownames(tax), rownames(samples)]

# Let us ensure that the data is in correct (numeric matrix) format:
counts <- as.matrix(counts)


## ----demodata, message=FALSE---------------------------
# coldata rownames match assay colnames
all(rownames(samples) == colnames(counts)) # our dataset
class(samples) # should be data.frame or DataFrame

# rowdata rownames match assay rownames
all(rownames(tax) == rownames(counts)) # our dataset
class(tax) # should be data.frame or DataFrame

# Counts 
class(counts) # should be a numeric matrix


## ----importingcsv3-------------------------------------
# Create a TreeSE
tse_taxa <- TreeSummarizedExperiment(assays =  SimpleList(counts = counts),
                                     colData = DataFrame(samples),
                                     rowData = DataFrame(tax))

tse_taxa


## ----importingcsv4, message=FALSE----------------------
count_file <- "data/assay_metabolites.csv"
sample_file <- "data/coldata.csv"

# Load files
counts  <- read.csv(count_file, row.names=1)  
samples <- read.csv(sample_file, row.names=1)

# Create a TreeSE for the metabolite data
tse_metabolite <- TreeSummarizedExperiment(assays = SimpleList(concs = as.matrix(counts)),
                                           colData = DataFrame(samples))

tse_metabolite


## ----importingcsv5-------------------------------------
# Create an ExperimentList that includes experiments
experiments <- ExperimentList(microbiome = tse_taxa, 
                              metabolite = tse_metabolite)

# Create a MAE
mae <- MultiAssayExperiment(experiments = experiments)

mae


## ------------------------------------------------------
biom_file_path <- "data/Aggregated_humanization2.biom"
sample_meta_file_path <- "data/Mapping_file_ADHD_aggregated.csv"
tree_file_path <- "data/Data_humanization_phylo_aggregation.tre"


## ------------------------------------------------------
library(mia)

# read biom and convert it to TreeSE
tse <- loadFromBiom(biom_file_path,
                    rankFromPrefix = TRUE,
                    removeTaxaPrefixes = TRUE)

# Check
tse


## ------------------------------------------------------
assay(tse, "counts")[1:3, 1:3]


## ------------------------------------------------------
head(rowData(tse))


## ------------------------------------------------------
# Genus level has additional '\"', so let's delete that also
rowdata_modified <- BiocParallel::bplapply(rowData(tse), 
                                           FUN = stringr::str_remove, 
                                           pattern = '\"')

# rowdata_modified is a list, so convert this back to DataFrame format. 
# and assign the cleaned data back to the TSE rowData
rowData(tse) <- DataFrame(rowdata_modified)

# Now we have a nicer table
head(rowData(tse))


## ------------------------------------------------------
head(colData(tse))


## ------------------------------------------------------
# CSV file with colnames in the first row and rownames in the first column
sample_meta <- read.csv(sample_meta_file_path,
                        sep = ",", row.names = 1)

# Add this sample data to colData of the taxonomic data object
# Note that the data must be given in a DataFrame format (required for our purposes)
colData(tse) <- DataFrame(sample_meta)


## ------------------------------------------------------
head(colData(tse))


## ------------------------------------------------------
# Reads the tree file
tree <- ape::read.tree(tree_file_path)

# Add tree to rowTree
rowTree(tse) <- tree

# Check
tse


## ---- eval=FALSE---------------------------------------
## head(rowTree(tse))


## ---- message=FALSE------------------------------------
library(mia)

# phyloseq example data
data(GlobalPatterns, package="phyloseq") 
GlobalPatterns_phyloseq <- GlobalPatterns
GlobalPatterns_phyloseq


## ---- message=FALSE------------------------------------
# convert phyloseq to TSE
GlobalPatterns_TSE <- makeTreeSummarizedExperimentFromPhyloseq(GlobalPatterns_phyloseq) 
GlobalPatterns_TSE


## ---- message=FALSE------------------------------------
# convert TSE to phyloseq
GlobalPatterns_phyloseq2 <- makePhyloseqFromTreeSummarizedExperiment(GlobalPatterns_TSE) 
GlobalPatterns_phyloseq2

