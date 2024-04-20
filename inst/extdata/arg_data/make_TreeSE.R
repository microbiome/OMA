library(mia)
library(dplyr)

## loading metadata
col_data <- read.table(
  "./sfile-microbiome-sample-detail.tsv",
  header = TRUE, sep = "\t", check.names = FALSE)
# Making a new ID variable to match the assembly data with the rest of the data
col_data <- col_data %>% mutate(tmp_ID = paste0(Study, "__", Sample_ID))
rownames(col_data) <- col_data$tmp_ID

## Loading the taxonomic abundance table
assay <- read.table(
  "./samples_with_1k_core_COG_ORFs.adult_stool.SGB_level.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# helper list for subsetting and renaming rownames
new_ids <- lapply(rownames(assay), function(old_id) {
  rownames(col_data[which(old_id==col_data$Sample_ID),]) }) %>%
  setNames(rownames(assay)) %>% unlist()
# subsetting
assay <- assay[names(new_ids),]
# renaming
rownames(assay) <- new_ids %>% unname()

## Loading rowData
row_data <- read.table(
  "./SGB_info.taxonomy.tsv",
  header = TRUE, sep = "\t",
  row.names = 1, # sgbid
  check.names = FALSE)
# subsetting to taxonomy ranks only
row_data <- row_data[,2:8]
colnames(row_data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                        "Species")

## common taxa
common_taxa <- intersect(rownames(row_data), colnames(assay))

## Loading Read-based normalized ARG abundance table
assay_read <- read.table(
  "./cpg_table.20220516_scg_valid_only.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
# helper list for subsetting and renaming rownames
new_ids <- lapply(rownames(assay_read), function(old_id) {
  rownames(col_data[which(old_id==col_data$Sample_ID),]) }) %>%
  setNames(rownames(assay_read)) %>% unlist()
# subsetting
assay_read <- assay_read[names(new_ids),]
# renaming
rownames(assay_read) <- new_ids %>% unname()

## Loading  assembly-based normalized ARG abundance table
assay_assembly <- read.table(
  "./DS3.SCG_normalized_ARG_abund.columns_CARD_ref.n_6104.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

## Loading rowData for arg read data
row_data_arg <- read.table(
  "ARG_family_info.tsv",
  header = TRUE, sep = "\t",
  row.names = 1, # ARG_family
  check.names = FALSE)
## common arg family
common_arg <- Reduce(intersect, list(rownames(row_data_arg),
                                     colnames(assay_read),
                                     colnames(assay_assembly)))

## Common ids
common_ids <- Reduce(intersect, list(rownames(assay_read),
                                     rownames(assay),
                                     rownames(assay_assembly),
                                     rownames(col_data)))

## Constructing the TreeSE object with the read assay and metadata, with an
# alternative experiment object for assembly data.

tse <- TreeSummarizedExperiment(
  assays=SimpleList(abundances=t(assay[common_ids, common_taxa])),
  rowData=row_data[common_taxa, ],
  colData=col_data[common_ids,]
)
altExp(tse, "read") <- SummarizedExperiment(
  assays = SimpleList(abundances=t(assay_read[common_ids, common_arg])),
  rowData = row_data_arg[common_arg, ])
altExp(tse, "assembly") <- SummarizedExperiment(
  assays = SimpleList(abundances=t(assay_assembly[common_ids, common_arg])),
  rowData = row_data_arg[common_arg, ])

## Saving the object
saveRDS(tse, "../arg_data.rds")
rm(common_ids, col_data, assay_assembly, assay_read, new_ids, assay, row_data,
   row_data_arg, common_arg, common_taxa)