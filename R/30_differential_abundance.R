## ----setup, echo=FALSE, results="asis"---------------------------------------------------------------------------------
library(rebook)
chapterPreamble()


## ----import-daa-data---------------------------------------------------------------------------------------------------
library(mia)
library(tidyverse)

# Import dataset
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

# Show patient status by cohort
table(tse$patient_status, tse$cohort) %>%
  knitr::kable()


## ----prep-daa-data-----------------------------------------------------------------------------------------------------
# Agglomerate by genus and subset by prevalence
tse <- subsetByPrevalentFeatures(tse,
                             rank = "Genus",
                             prevalence = 10 / 100)

# Transform count assay to relative abundances
tse <- transformAssay(tse,
                      assay.type = "counts",
                      method = "relabundance")


## ----run-aldex2--------------------------------------------------------------------------------------------------------
# Load package
library(ALDEx2)

# Generate Monte Carlo samples of the Dirichlet distribution for each sample.
# Convert each instance using the centered log-ratio transform.
# This is the input for all further analyses.
set.seed(123)
x <- aldex.clr(assay(tse), tse$patient_status)     


## ----aldex2-ttest------------------------------------------------------------------------------------------------------
# calculates expected values of the Welch's t-test and Wilcoxon rank
# test on the data returned by aldex.clr
x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)


## ----aldex2-effect-----------------------------------------------------------------------------------------------------
# Determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)

# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)


## ----plot-aldex2-------------------------------------------------------------------------------------------------------
par(mfrow = c(1, 2))

aldex.plot(aldex_out,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)


## ----aldex2-res--------------------------------------------------------------------------------------------------------
aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()


## ----run-ancombc, warning=FALSE----------------------------------------------------------------------------------------
# Load package
library(ANCOMBC)

# Run ANCOM-BC at the genus level and only including the prevalent genera
ancombc2_out <- ancombc2(data = tse,
                         assay.type = "counts",
                         fix_formula = "patient_status",
                         p_adj_method = "fdr",
                         prv_cut = 0,
                         group = "patient_status",
                         struc_zero = TRUE,
                         neg_lb = TRUE,
                         # multi group comparison is deactivated automatically
                         global = TRUE)


## ----ancombc-res-------------------------------------------------------------------------------------------------------
# store the FDR adjusted results 
ancombc2_out$res %>%
  dplyr::select(taxon, lfc_patient_statusControl, q_patient_statusControl) %>%
  filter(q_patient_statusControl < 0.05) %>%
  arrange(q_patient_statusControl) %>%
  head() %>%
  knitr::kable()


## ----run-maaslin2, warning=FALSE, results="hide"-----------------------------------------------------------------------
# Load package
library(Maaslin2)

# maaslin expects features as columns and samples as rows 
# for both the abundance table as well as metadata 

# We can specify different GLMs/normalizations/transforms.
# specifying a ref is especially important if you have more than 2 levels
maaslin2_out <- Maaslin2(input_data = as.data.frame(t(assay(tse))),
                         input_metadata = as.data.frame(colData(tse)),
                         output = "DAA example",
                         transform = "AST",
                         fixed_effects = "patient_status",
                         # you can also fit MLM by specifying random effects
                         # random_effects = c(...),
                         reference = "patient_status,Control",
                         normalization = "TSS",
                         standardize = FALSE,
                         # filtering was previously performed
                         min_prevalence = 0)


## ---- maaslin2-res-----------------------------------------------------------------------------------------------------
maaslin2_out$results %>%
  filter(qval < 0.05) %>%
  knitr::kable()


## ----run-linda---------------------------------------------------------------------------------------------------------
# Load package
library(MicrobiomeStat)

# Run LinDA
linda_out <- linda(feature.dat = as.data.frame(assay(tse)),
                   meta.dat = as.data.frame(colData(tse)),
                   formula = "~ patient_status",
                   alpha = 0.05,
                   prev.filter = 0,
                   mean.abund.filter = 0)


## ----linda-res---------------------------------------------------------------------------------------------------------
# List genera for which H0 could be rejected:
linda_out$output$patient_statusControl %>%
  filter(reject) %>%
  dplyr::select(stat, padj) %>%
  rownames_to_column(var = "feature") %>%
  knitr::kable()


## ----run-zicoseq-------------------------------------------------------------------------------------------------------
# Load package
library(GUniFrac)

set.seed(123)
zicoseq_out <- ZicoSeq(feature.dat = as.matrix(assay(tse)),
                       meta.dat = as.data.frame(colData(tse)),
                       grp.name = "patient_status",
                       feature.dat.type = "count",
                       return.feature.dat = TRUE,
                       prev.filter = 0,
                       mean.abund.filter = 0,
                       max.abund.filter = 0,
                       perm.no = 999)


## ----zicoseq-res-------------------------------------------------------------------------------------------------------
zicoseq_res <- cbind.data.frame(p.raw = zicoseq_out$p.raw,
                                p.adj.fdr = zicoseq_out$p.adj.fdr)

zicoseq_res %>%
  filter(p.adj.fdr < 0.05) %>%
  arrange(p.adj.fdr) %>%
  knitr::kable()


## ----plot-zicoseq------------------------------------------------------------------------------------------------------
## x-axis is the effect size: R2 * direction of coefficient
ZicoSeq.plot(ZicoSeq.obj = zicoseq_out,
             pvalue.type = 'p.adj.fdr')


## ----daa-data-libsize--------------------------------------------------------------------------------------------------
# Compute and store library size in colData
colData(tse)$library_size <- colSums(assay(tse, "counts"))


## ----run-adj-ancombc, warning=FALSE------------------------------------------------------------------------------------
# perform the analysis 
ancombc2_out <- ancombc2(tse,
                         assay.type = "counts",
                         fix_formula = "patient_status + cohort + library_size",
                         p_adj_method = "fdr",
                         lib_cut = 0,
                         group = "patient_status", 
                         struc_zero = TRUE, 
                         neg_lb = TRUE,
                         alpha = 0.05,
                         # multi-group comparison is deactivated automatically
                         global = TRUE)


## ----adj-ancombc-res---------------------------------------------------------------------------------------------------
ancombc2_out$res %>%
  dplyr::select(starts_with(c("taxon", "lfc", "q"))) %>%
  arrange(q_patient_statusControl) %>%
  head() %>%
  knitr::kable()


## ----run-adj-linda-----------------------------------------------------------------------------------------------------
linda_out <- linda(as.data.frame(assay(tse, "counts")),
                   as.data.frame(colData(tse)),
                   formula = "~ patient_status + cohort + library_size",
                   alpha = 0.05,
                   prev.filter = 0,
                   mean.abund.filter = 0)


## ----adj-linda-res-----------------------------------------------------------------------------------------------------
# Select results for the patient status
linda_res <- linda_out$output$patient_statusControl

linda_res %>%
  filter(reject) %>%
  dplyr::select(log2FoldChange, stat, padj) %>%
  rownames_to_column(var = "feature") %>%
  head() %>%
  knitr::kable()


## ----run-adj-zicoseq---------------------------------------------------------------------------------------------------
set.seed(123)
zicoseq_out <- ZicoSeq(feature.dat = as.matrix(assay(tse)),
                       meta.dat = as.data.frame(colData(tse)),
                       grp.name = "patient_status",
                       adj.name = c("cohort", "library_size"), 
                       feature.dat.type = "count",
                       return.feature.dat = TRUE,
                       prev.filter = 0,
                       mean.abund.filter = 0,
                       max.abund.filter = 0,
                       perm.no = 999)


## ----adj-zicoseq-res---------------------------------------------------------------------------------------------------
zicoseq_res <- cbind.data.frame(p.raw = zicoseq_out$p.raw,
                                p.adj.fdr = zicoseq_out$p.adj.fdr)

zicoseq_res %>%
  filter(p.adj.fdr < 0.05) %>%
  head() %>%
  knitr::kable()

