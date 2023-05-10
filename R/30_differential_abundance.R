## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----load-pkg-data------------------------------------------------------------
library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(knitr)
library(tidyverse)
library(ALDEx2)
library(Maaslin2)
library(MicrobiomeStat)
library(ANCOMBC)
library(GUniFrac)

# set random seed because some tools can randomly vary and then produce 
# different results:
set.seed(13253)

# we use a demo dataset and restrict it to two geo locations
# for easy illustration
data(peerj13075)
tse0 <- peerj13075
tse0 <- tse0[ ,tse0$Geographical_location %in% c("Pune", "Nashik")]
# Let us make this a factor
tse0$Geographical_location <- factor(tse0$Geographical_location)

# how many observations do we have per group?
count(as.data.frame(colData(tse0)), Geographical_location) %>% kable()


## -----------------------------------------------------------------------------
tse <- agglomerateByRank(tse0, rank = "genus") %>%
       transformCounts(assay.type = "counts", method = "relabundance", MARGIN = "samples") %>%
       subsetByPrevalentTaxa(detection = 0, prevalence = 10/100, assay.type = "relabundance")

# Add also clr abundances
tse <- transformCounts(tse, method="clr", pseudocount=1) # not bale to run


## ---- aldex2, eval=TRUE-------------------------------------------------------
# Generate Monte Carlo samples of the Dirichlet distribution for each sample.
# Convert each instance using the centered log-ratio transform.
# This is the input for all further analyses.
set.seed(254)
x <- aldex.clr(assay(tse), tse$Geographical_location)     


## ---- aldex2_ttest, eval=TRUE-------------------------------------------------
# calculates expected values of the Welch's t-test and Wilcoxon rank
# test on the data returned by aldex.clr
x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)


## ---- aldex2_efs, eval=TRUE---------------------------------------------------
# determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)

# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)


## ---- eval=TRUE---------------------------------------------------------------
par(mfrow = c(1, 2))
  aldex.plot(
    aldex_out, 
    type = "MA", 
    test = "welch", 
    xlab = "Log-ratio abundance",
    ylab = "Difference",
    cutoff = 0.05
  )
  aldex.plot(
    aldex_out, 
    type = "MW", 
    test = "welch",
    xlab = "Dispersion",
    ylab = "Difference",
    cutoff = 0.05
  )


## ---- eval=TRUE---------------------------------------------------------------
rownames_to_column(aldex_out, "genus") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than tt
  dplyr::select(genus, we.eBH, wi.eBH, effect, overlap) %>%
  kable()


## ----ancombc2, warning = FALSE, eval=TRUE-------------------------------------
# Run ANCOM-BC 
out <- ancombc2(
  data = tse,
  assay_name = "counts", 
  tax_level = "genus", 
  fix_formula = "Geographical_location", 
  p_adj_method = "fdr", 
  lib_cut = 0,
  prv_cut = 0,
  group = "Geographical_location", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)

# store the FDR adjusted results [test on v2.0.3] 
ancombc_result <- cbind.data.frame(taxid = out$res$taxon,
                                   ancombc = as.vector(out$res$q_Geographical_locationPune))

# store the FDR adjusted results [test on v1.2.2] 
# ancombc_result <- out$res$q_val %>% rownames_to_column('taxid') %>% dplyr::rename(ancombc = 2)


## ---- eval=TRUE---------------------------------------------------------------
kable(head(ancombc_result))


## ----maaslin2, results = "hide", eval=TRUE------------------------------------
# maaslin expects features as columns and samples as rows 
# for both the abundance table as well as metadata 

# We can specify different GLMs/normalizations/transforms.
# Let us use similar settings as in Nearing et al. (2021):
maaslin2_out <- Maaslin2(
  t(assay(tse)),
  data.frame(colData(tse)),
  output = "DAA example",
  transform = "AST",
  fixed_effects = "Geographical_location",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "Geographical_location,Pune",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)


## ---- maaslin2kable, eval=TRUE------------------------------------------------
kable(head(filter(maaslin2_out$results, qval <= 0.05)))


## ----linda, eval=TRUE---------------------------------------------------------
meta <- as.data.frame(colData(tse)) %>% dplyr::select(Geographical_location)
linda.res <- linda(
  as.data.frame(assay(tse)), 
  meta, 
  formula = '~Geographical_location', 
  alpha = 0.05, 
  prev.filter = 0, 
  mean.abund.filter = 0)

linda_out <- linda.res$output$Geographical_locationPune
# to scan the table for genera where H0 could be rejected:
kable(head(filter(as.data.frame(linda_out), reject)))


## ----ZicoSeq, eval=TRUE-------------------------------------------------------
set.seed(123)
meta <- as.data.frame(colData(tse))
zicoseq.obj <- GUniFrac::ZicoSeq(meta.dat = meta, 
                                 feature.dat = as.matrix(assay(tse)),
                                 grp.name = 'Geographical_location',
                                 adj.name = NULL, 
                                 feature.dat.type = 'count',
                                 prev.filter = 0,
                                 perm.no = 999,
                                 mean.abund.filter = 0,
                                 max.abund.filter = 0,
                                 return.feature.dat = T)
zicoseq_out <- cbind.data.frame(p.raw=zicoseq.obj$p.raw, p.adj.fdr=zicoseq.obj$p.adj.fdr) 
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))


## ----ZicoSeqplot, eval=TRUE---------------------------------------------------
## x-axis is the effect size: R2 * direction of coefficient
ZicoSeq.plot(ZicoSeq.obj = zicoseq.obj, meta.dat = meta, pvalue.type ='p.adj.fdr')


## ----comparison, eval=TRUE----------------------------------------------------
aldex_result <- rownames_to_column(aldex_out, "taxid") %>% dplyr::select(taxid, aldex2 = wi.eBH)
summ <- full_join(aldex_result, ancombc_result, by = "taxid") %>%
  full_join(
    dplyr::select(maaslin2_out$results, taxid = feature, maaslin2 = qval), 
    by = "taxid") %>%
    full_join(linda_out %>% rename(LinDA=padj) %>% dplyr::select(LinDA)%>% rownames_to_column('taxid') ) %>%
  full_join(zicoseq_out %>% dplyr::select(p.adj.fdr) %>% rename(ZicoSeq = p.adj.fdr) %>% rownames_to_column('taxid')) %>%
  mutate(
    across(c(aldex2: ZicoSeq), ~ .x <= 0.05),
    # the following line would be necessary without prevalence filtering 
    # as some methods output NA
    #across(-taxid, function(x) ifelse(is.na(x), FALSE, x)),
    score = rowSums(across(c(aldex2, ancombc, maaslin2, LinDA, ZicoSeq)))
  )

# Mark all NAs as FALSE
summ[is.na(summ)] <- FALSE

# This is how it looks like:
kable(head(summ))


## ---- eval=TRUE---------------------------------------------------------------
# how many genera were identified by each method?
summarise(summ, across(where(is.logical), sum)) %>%
  kable()
# which genera are identified by all methods?
filter(summ, score == 5) %>% kable()


## ---- daplotting, eval=TRUE, fig.width=20,fig.height=4------------------------
# Create a jittered boxplot for each genus 
assay.type <- "relabundance"
plot_data <- data.frame(t(assay(tse, assay.type)))
plot_data$Geographical_location <- tse$Geographical_location
plots <- pmap(dplyr::select(summ, taxid, score), function(taxid, score) {
  ggplot(plot_data, aes_string(x="Geographical_location", y=taxid)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    # scale_y_log10() + # log trans will cause 0 values missing 
    scale_y_sqrt() + 
    labs(title=glue::glue("{taxid}"), x="", y=glue::glue("Abundance ({assay.type})")) +    
    theme_bw() +
    theme(legend.position = "none")
})

# now we can show only those genera that have at least score 5 (or 4 or 3 or 2 or 1)
robust_plots <- plots[summ$score == 4 & !is.na(summ$score)] 

# to display this nicely in the book we use patchwork here:
# (we show first ones)
robust_plots[[1]] + 
  robust_plots[[2]] + 
  robust_plots[[3]] + 
  robust_plots[[4]] +
  robust_plots[[5]] +
  plot_layout(nrow = 1)
  
# or if we have most trust in any specific method we can show genera that 
# are differentially abundant according to that method and then look in the
# title how many methods also identified it (we only show first 6 here):
ancombc_plots <- plots[summ$ancombc & !is.na(summ$score)] 
ancombc_plots[[1]] + 
  ancombc_plots[[2]] + 
  ancombc_plots[[3]] + 
  ancombc_plots[[4]] +
  ancombc_plots[[5]] +
  ancombc_plots[[6]] +
  plot_layout(nrow = 1)


## ----ancombc_adj, eval=TRUE---------------------------------------------------
# perform the analysis 
ancombc_cov <- ancombc2(
  data = tse,
  assay_name = "counts",
  tax_level = "genus",
  fix_formula = "Geographical_location + Age", 
  p_adj_method = "fdr", 
  lib_cut = 0, 
  group = "Geographical_location", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
# now the model answers the question: holding Age constant, are 
# bacterial taxa differentially abundant? Or, if that is of interest,
# holding phenotype constant, is Age associated with bacterial abundance?
# Again we only show the first 6 entries.
kable(head(ancombc_cov$res$q_val))


## ----linda_adj, eval=TRUE-----------------------------------------------------
otu.tab <- as.data.frame(assay(tse))
meta <- as.data.frame(colData(tse))
linda_cov <- linda(
  otu.tab, 
  meta, 
  formula = '~ Geographical_location + Age', 
  alpha = 0.05, 
  prev.filter = 0, 
  mean.abund.filter = 0)
linda.res <- linda_cov$output$Geographical_locationPune
kable(head(filter(linda.res, reject==T)))



## ----ZicoSeq_adj, eval=TRUE---------------------------------------------------
set.seed(123)
zicoseq.obj <- GUniFrac::ZicoSeq(meta.dat = as.data.frame(colData(tse)) , 
                                 feature.dat = as.matrix(assay(tse)),
                                 grp.name = 'Geographical_location',
                                 adj.name = 'Gender', 
                                 feature.dat.type = 'count',
                                 prev.filter = 0,
                                 perm.no = 999,
                                 mean.abund.filter = 0,
                                 max.abund.filter = 0,
                                 return.feature.dat = T)
zicoseq_out <- cbind.data.frame(p.raw=zicoseq.obj$p.raw, p.adj.fdr=zicoseq.obj$p.adj.fdr) 
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))


## ----sessionInfo, echo=FALSE, results='asis'----------------------------------
prettySessionInfo()

