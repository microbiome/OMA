## ----setup, echo=FALSE, results="asis"----------------------------------------
library(rebook)
chapterPreamble()


## ----load-pkg-data------------------------------------------------------------
library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(ALDEx2)
library(Maaslin2)
library(MicrobiomeStat)
library(knitr)
library(tidyverse)
library(ANCOMBC)

# set random seed because some tools can randomly vary and then produce 
# different results:
set.seed(13253)

# we use a demo dataset and restrict it to two geo locations
# for easy illustration
data(peerj13075)
tse <- peerj13075
tse <- tse[ ,tse$Geographical_location %in% c("Pune", "Nashik")]
# Let us make this a factor
tse$Geographical_location <- factor(tse$Geographical_location)

# how many observations do we have per group?
count(as.data.frame(colData(tse)), Geographical_location) %>% kable()


## -----------------------------------------------------------------------------
tse <- subsetByPrevalentTaxa(tse, detection = 0, prevalence = 0.1)


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
  select(genus, we.eBH, wi.eBH, effect, overlap) %>%
  kable()


## ----ancombc2, warning = FALSE, eval=TRUE-------------------------------------
# perform the analysis 
out <- ancombc2(
  data = tse,
  tax_level="genus",
  fix_formula = "Geographical_location", 
  p_adj_method = "fdr", 
  prv_cut = 0, # prev filtering has been done above already
  lib_cut = 0, 
  group = "Geographical_location", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 20), # use max_iter >= 100 on real data 
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
# store the results in res 
res <- out$res


## ---- eval=TRUE---------------------------------------------------------------
kable(head(res))


## ----maaslin2, results = "hide", eval=TRUE------------------------------------
# maaslin expects features as columns and samples as rows 
# for both the asv/otu table as well as meta data 
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
# you can specifiy different GLMs/normalizations/transforms. We used similar
# settings as in Nearing et al. (2021) here:
fit_data <- Maaslin2(
  asv,
  meta_data,
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
kable(head(filter(fit_data$results, qval <= 0.05)))


## ----linda, eval=TRUE---------------------------------------------------------
otu.tab <- as.data.frame(assay(tse))
meta <- as.data.frame(colData(tse)) %>% select(Geographical_location)
res <- linda(
  otu.tab, 
  meta, 
  formula = '~Geographical_location', 
  alpha = 0.05, 
  prev.filter = 0, 
  mean.abund.filter = 0)

# to scan the table for genera where H0 could be rejected:
kable(head(filter(as.data.frame(res$output$Geographical_locationPune), reject)))


## ----linda2, eval=TRUE--------------------------------------------------------
# change genus names to otu ids for ancombc results to make it joinable with others
id_switch <- as.data.frame(rowData(tse)) %>%
  rownames_to_column("taxid") %>%
  select(taxid, genus)
abc_res <- select(out$res, genus = taxon, ancombc = diff_Geographical_locationPune) %>%
  left_join(id_switch, by = "genus") %>%
  select(-genus)

# join all results together
summ <- full_join(
    rownames_to_column(aldex_out, "taxid") %>%
      select(taxid, aldex2 = wi.eBH),
    abc_res,
    by = "taxid") %>%
  full_join(
    select(fit_data$results, taxid = feature, maaslin2 = qval), 
    by = "taxid") %>%
    full_join(
      rownames_to_column(as.data.frame(res$output$Geographical_locationPune), "taxid") %>%
        select(taxid, LinDA = reject), 
      by = "taxid") %>%
  mutate(
    across(c(aldex2, maaslin2), ~ .x <= 0.05),
    # the following line would be necessary without prevalence filtering 
    # as some methods output NA
    #across(-genus, function(x) ifelse(is.na(x), FALSE, x)),
    ancombc = ifelse(is.na(ancombc), FALSE, ancombc),
    score = rowSums(across(c(aldex2, ancombc, maaslin2, LinDA))),
  )

# This is how it looks like:
kable(head(summ))


## ---- eval=TRUE---------------------------------------------------------------
# how many genera were identified by each method?
summarise(summ, across(where(is.logical), sum)) %>%
  kable()
# which genera are identified by all methods?
filter(summ, score == 4) %>% kable()


## ---- daplotting, eval=TRUE, fig.width=20,fig.height=4------------------------
# Data
data(peerj13075)
tse <- peerj13075

# Add relative abundances and clr abundances
tse <- transformCounts(tse, method="relabundance")
tse <- transformCounts(tse, method="clr", pseudocount=1)

# Subset to prevalent taxa (exclude rare taxa at 10 percent prevalence using 0 detection threshold):
# do the subsetting based on the relative abundance assay
tse <- subsetByPrevalentTaxa(tse, detection = 0, prevalence = 10/100, assay.type="relabundance")

# Subset to certain geolocations
tse <- tse[ ,tse$Geographical_location %in% c("Pune", "Nashik")]

# Let us make the geo location a factor
tse$Geographical_location <- factor(tse$Geographical_location)

# Create a jittered boxplot for each genus 
assay.type <- "relabundance"
plot_data <- data.frame(t(assay(tse, assay.type)))
plot_data$Geographical_location <- tse$Geographical_location
plots <- pmap(select(summ, taxid, score), function(taxid, score) {
  ggplot(plot_data, aes_string(x="Geographical_location", y=taxid)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    scale_y_log10() + 
    labs(title=glue::glue("{taxid}"), x="", y=glue::glue("Abundance ({assay.type})")) +    
    theme_bw() +
    theme(legend.position = "none")
})

# now we can show only those genera that have at least score 3 (or 2 or 1)
robust_plots <- plots[summ$score == 4 & !is.na(summ$score)] 

# to display this nicely in the book we use patchwork here:
# (we show first ones)
robust_plots[[1]] + 
  robust_plots[[2]] + 
  robust_plots[[3]] + 
  robust_plots[[4]] +
  robust_plots[[5]] +
  robust_plots[[6]] +
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
  ancombc_plots[[6]] 


## ---- eval=TRUE---------------------------------------------------------------
# FIXME: switch to a faster example / method
out_cov = ancombc2(
  data = tse, 
  fix_formula = "Geographical_location + Age", # here we add Age to the model
  p_adj_method = "fdr", 
  prv_cut = 0,  # we did that already
  lib_cut = 0, 
  group = "Geographical_location",
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 20),
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)

# now the model answers the question: holding Age constant, are 
# bacterial taxa differentially abundant? Or, if that is of interest,
# holding phenotype constant, is Age associated with bacterial abundance?
# Again we only show the first 6 entries.
kable(head(out_cov$res))


## ----sessionInfo, echo=FALSE, results='asis'----------------------------------
prettySessionInfo()

