# Differential abundance {#differential-abundance}

<script>
document.addEventListener("click", function (event) {
    if (event.target.classList.contains("rebook-collapse")) {
        event.target.classList.toggle("active");
        var content = event.target.nextElementSibling;
        if (content.style.display === "block") {
            content.style.display = "none";
        } else {
            content.style.display = "block";
        }
    }
})
</script>

<style>
.rebook-collapse {
  background-color: #eee;
  color: #444;
  cursor: pointer;
  padding: 18px;
  width: 100%;
  border: none;
  text-align: left;
  outline: none;
  font-size: 15px;
}

.rebook-content {
  padding: 0 18px;
  display: none;
  overflow: hidden;
  background-color: #f1f1f1;
}
</style>



## Differential abundance analysis

This section provides an overview and examples of *differential
abundance analysis (DAA)* based on one of the [openly available
datasets](https://microbiome.github.io/mia/reference/mia-datasets.html)
in mia to illustrate how to perform differential abundance analysis
(DAA). DAA identifies differences in the abundances of individual
taxonomic groups between two or more groups (e.g. treatment vs
control). This can be performed at any phylogenetic level.

We perform DAA to identify biomarkers and/or gain understanding of a
complex system by looking at its isolated components. For example,
identifying that a bacterial taxon is different between a patient
group with disease *X* vs a healthy control group might lead to
important insights into the pathophysiology. Changes in the microbiota
might be cause or a consequence of a disease. Either way, it can
help to understand the system as a whole. Be aware that this approach
has also been criticized recently [@Quinn2021].


### Examples and tools

There are many tools to perform DAA. The most popular tools, without going into
evaluating whether or not they perform well for this task, are:  

- ALDEx2 [@Gloor2016] 
- ANCOMBC [@ancombc2020]
- corncob [@Martin2021]
- DESeq2 [@Love2014] 
- edgeR [@Chen2016]
- lefser [@Khlebrodova2021]
- Maaslin2 [@Mallick2020]
- metagenomeSeq [@Paulson2017]
- limma [@Ritchie2015]
- LinDA [@Zhou2022]
- [t-test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test)  
- [Wilcoxon test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test)  


We recommend to have a look at @Nearing2022
who compared all these listed methods across 38
different datasets. Because different methods use different approaches
(parametric vs non-parametric, different normalization techiniques, assumptions
etc.), results can differ between methods. 
Unfortunately, as @Nearing2022 point out, they
can differ substantially. Therefore, it is highly recommended to pick several
methods to get an idea about how robust and potentially reproducible your
findings are depending on the method. In this section we demonstrate 4 methods
that can be recommended based on recent literature (ANCOM-BC, ALDEx2, Maaslin2
and LinDA) and we will compare the results between them.
Note that the purpose of this section is to show how to perform DAA in R, not
how to correctly do causal inference. Depending on your experimental setup
and your theory, you must determine how to specify any model exactly. 
E.g., there might be confounding factors that might drive (the absence of)
differences between the shown groups that we ignore here for simplicity.
However, we will show how you could include covariates in those models.
Furthermore, we picked a dataset that merely has microbial abundances in a TSE
object as well as a grouping variable in the sample data. We simplify the
analysis by only including 2 of the 3 groups. 





```r
library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(ALDEx2)
library(Maaslin2)
library(MicrobiomeStat)
library(knitr)
library(tidyverse)

# set random seed because some tools can randomly vary and then produce 
# different results:
set.seed(13253)

# For ANCOMBC we need development version
# Can be installed with:
# remotes::install_github("FrederickHuangLin/ANCOMBC")
library(ANCOMBC)

# we use a demo dataset and restrict it to two geo locations
# for easy illustration
data(peerj13075)
tse <- peerj13075
tse <- tse[ ,tse$Geographical_location %in% c("Pune", "Nashik")]
# Let us make this a factor
tse$Geographical_location <- factor(tse$Geographical_location)

# how many observations do we have per group?
count(as.data.frame(colData(tse)), Geographical_location) %>% kable()
```


\begin{tabular}{l|r}
\hline
Geographical\_location & n\\
\hline
Nashik & 11\\
\hline
Pune & 36\\
\hline
\end{tabular}

### Prevalence Filtering 

Before we jump to our analyses, we may want to perform prevalence filtering.
@Nearing2022 found that applying a 10% threshold
for the prevalence of the taxa generally resulted in more robust results. 
Some tools have builtin arguments for that. By applying the threshold to our
input data, we can make sure it is applied for all tools. Below we show how to
do this in `mia`:


```r
tse <- subsetByPrevalentTaxa(tse, detection = 0, prevalence = 0.1)
```


### ALDEx2

In this section, we will show how to perform a simple ALDEx2 analysis. 
If you wanted to pick a single method, this method could be recommended to use.
According to the developers experience, it tends to identify the common
features identified by other methods. This statement is in line with a recent
independent evaluation by @Nearing2022.  
Please also have a look at the more extensive 
[vignette](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html) 
that covers this flexible tool in more depth. ALDEx2 estimates technical
variation within each sample per taxon by utilizing the Dirichlet distribution.
It furthermore applies the centered-log-ratio transformation (or closely
related log-ratio transforms). Depending on the experimental setup, it will
perform a two sample Welch's T-test and Wilcoxon-test or a one-way ANOVA and
Kruskal-Wallis-test. For more complex study designs, there is a possibility to 
utilize the `glm` functionality within ALDEx2. The Benjamini-Hochberg procedure
is applied in any case to correct for multiple testing. Below we show a simple
example that illustrates the workflow.



```r
# Generate Monte Carlo samples of the Dirichlet distribution for each sample.
# Convert each instance using the centered log-ratio transform.
# This is the input for all further analyses.
set.seed(254)
x <- aldex.clr(assay(tse), tse$Geographical_location)     
```


The t-test:


```r
# calculates expected values of the Welch's t-test and Wilcoxon rank
# test on the data returned by aldex.clr
x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)
```


Effect sizes:


```r
# determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)

# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)
```

Now, we can create a so called Bland-Altman or MA plot (left). It shows the
association between the relative abundance and the magnitude of the difference
per sample. Next to that, we can also create a plot that shows the dispersion
on the x-axis instead of log-ratio abundance. Red dots represent genera that are
differentially abundant ($q \leq 0.1$) between the 2 groups. Black points are
rare taxa and grey ones are abundant taxa. The dashed line represent an effect
size of 1. See @Gloor2016 to learn more about these plots.


```r
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
```

The evaluation as differential abundant in above plots is based on the
corrected pvalue. According to the ALDEx2 developers, the safest
approach is to identify those features where the 95% CI of the effect
size does not cross 0. As we can see in below table, this is not the
case for any of the identified genera (see overlap column, which
indicates the proportion of overlap). Also, the authors recommend to
focus on effect sizes and CIs rather than interpreting the pvalue. To
keep the comparison simple, we will here use the pvalue as decision
criterion. But please be aware that the effect size together with the
CI is a better answer to the question we are typically interested in
(see also [this
article](https://www.nature.com/articles/d41586-019-00857-9)).



```r
rownames_to_column(aldex_out, "genus") %>%
  filter(wi.eBH <= 0.05)  %>% # here we chose the wilcoxon output rather than tt
  select(genus, we.eBH, wi.eBH, effect, overlap) %>%
  kable()
```

### ANCOM-BC

The analysis of composition of microbiomes with bias correction
(ANCOM-BC) [@Das2020] is a recently developed method for differential
abundance testing. It is based on an earlier published approach
[@Mandal2015].  The previous version of ANCOM was among the methods
that produced the most consistent results and is probably a
conservative approach [@Nearing2022].  However, the new ANCOM-BC
method operates quite differently compared to the former ANCOM method.

As the only method, ANCOM-BC incorporates the so called *sampling
fraction* into the model. The latter term could be empirically
estimated by the ratio of the library size to the microbial
load. According to the authors, variations in this sampling fraction
would bias differential abundance analyses if ignored.  Furthermore,
this method provides p-values and confidence intervals for each
taxon. It also controls the FDR and it is computationally simple to
implement.

As we will see below, to obtain results, all that is needed is to pass
a tse object to the `ancombc()` function. Below, we first specify a
formula.  In this formula, other covariates could potentially be
included to adjust for confounding. We show this further below.
Please check the [function documentation](https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc.html)
to learn about the additional arguments that we specify below.

We recommend to use  or higher in real case studies.


```r
# perform the analysis 
out <- ancombc2(
  data = tse,
  tax_level="family",
  fix_formula = "Geographical_location", 
  p_adj_method = "fdr", 
  prv_cut = 0, # prev filtering has been done above already
  lib_cut = 0, 
  group = "Geographical_location", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 20),  
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
# store the results in res 
res <- out$res
```

The object `out` contains all model output. Again, see the 
[documentation of the function](https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc.html) 
under **Value** for an explanation of all the output objects. Our question
whether taxa are differentially abundant can be answered by looking at the
`res` object, which now contains dataframes with the coefficients, 
standard errors, p-values and q-values. Conveniently, there is a dataframe
`diff_abn`. Here, for each taxon it is indicated whether it is differentially
abundant between the groups (again, keep in mind that the answer is not 
black-white). Below we show the first 6 entries of this dataframe:  


```r
kable(head(res))
```



### MaAsLin2 

Next, we will illustrate how to use MaAsLin2, which is the next generation of
MaAsLin. As it is based on generalized linear models, it is flexible for different study designs and covariate
structures. The official package tutorial can be found [here](https://github.com/biobakery/biobakery/wiki/maaslin2). 


```r
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
```

Which genera are identified as differentially abundant? (leave out "head" to see all).


```r
kable(head(filter(fit_data$results, qval <= 0.05)))
```

A folder will be created that is called like the above specified
output.  It contains also figures to visualize the difference between
genera for the significant ones.


### LinDA 

Lastly, we cover linear models for differential abundance analysis of
microbiome compositional data (@Zhou2022). This tool is very similar
to ANCOMBC with few differences: 1) LinDA correct for the
compositional bias differently using the mode of all regression
coefficients. 2) The authors claim that it runs 100-1000x faster than
ANCOMBC and 3) it support hierarchical models.  The latter could be
ignored as ANCOMBC will be supporting hierarchical models with the
next release. Nevertheless, LinDA seems a promising tool that achieves
the best power/fdr trade-off together with ANCOMBC according to the
authors. The speed might make it the choice for bigger datasets or
datasets with a very high number of features.



```r
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
kable(head(filter(as.data.frame(res$output), Geographical_locationPune.reject)))
```


### Comparison of the methods

When we compare the methods in the context of a research question, we could
look at e.g. at whether they agree based on the applied decision criterion
(e.g. adjusted p value < 0.05). That is what we illustrate here. First we will 
look at how many taxa were identified by each method to begin with. In the next
step we will look at the intersection of identified taxa. To achieve that, we
first create a dataframe that summarises the decision criterion for each method
and shows a score from 0 to 3 indicating how many methods agreed on a particular
taxon.




```r
# Rename "taxon" column from ancombc results  so that it match with others
colnames(out$res)[colnames(out$res) == "taxon"] <- "genus"

summ <- full_join(
    rownames_to_column(aldex_out, "genus") %>%
      select(genus, aldex2 = wi.eBH),
    out$res %>%
      select(genus, ancombc = diff_Geographical_locationPune),
    by = "genus") %>%
  full_join(
    select(fit_data$results, genus = feature, maaslin2 = qval), 
    by = "genus") %>%
    full_join(
      rownames_to_column(as.data.frame(res$output), "genus") %>%
        select(genus, LinDA = Geographical_locationPune.reject), 
      by = "genus") %>%
  mutate(
    across(c(aldex2, maaslin2), ~ .x <= 0.05),
    # the following line would be necessary without prevalence filtering 
    # as some methods output NA
    #across(-genus, function(x) ifelse(is.na(x), FALSE, x)),
    score = rowSums(across(c(aldex2, ancombc, maaslin2, LinDA)))
  )

# This is how it looks like:
kable(head(summ))
```

Now we can answer our questions:


```r
# how many genera were identified by each method?
summarise(summ, across(where(is.logical), sum)) %>%
  kable()
# which genera are identified by all methods?
filter(summ, score == 4) %>% kable()
```

We see that each method identified at least 9 genera as differentially
abundant. Eight of those that were identified by ALDEx2,
were also identified by the other methods. We could plot the data for
any method or for those taxa that were identified by all methods:




```r
plot_data <- data.frame(t(assay(tse)))
plot_data$Geographical_location <- tse$Geographical_location
# create a plot for each genus where the score is indicated in the title
plots <- pmap(select(summ, genus, score), function(genus, score) {
  ggplot(plot_data, aes("Geographical_location", genus)) +
    geom_boxplot(aes(fill = Geographical_location), outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ggtitle(glue::glue("Score {score}")) +
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
  plot_layout(nrow = 2)
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
```



### Confounding variables

To perform causal inference, it is crucial that the method is able to include
covariates in the model. This is not possible with e.g. the Wilcoxon test.
Other methods such as both ANCOM methods, ALDEx2, LinDA, MaAsLin2 and others
allow this. Below we show how to include a covariate in ANCOM-BC.
It is very similar for all the methods that allow this. Since in this dataset
there are no covariates, I first simulate a new variable and add it to the TSE
object.



```r
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
```

In the next section of this book chapter we cover methods that can also take
into account the phylogenetic information of bacterial taxa to perform 
group-wise associations.


## Tree-based methods

### Group-wise associations testing based on balances

For testing associations based on balances, check the philr R/Bioconductor package.


## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.2.1 (2022-06-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ANCOMBC_2.1.2                  lubridate_1.9.2               
 [3] forcats_1.0.0                  stringr_1.5.0                 
 [5] dplyr_1.1.0                    purrr_1.0.1                   
 [7] readr_2.1.4                    tidyr_1.3.0                   
 [9] tibble_3.1.8                   ggplot2_3.4.1                 
[11] tidyverse_2.0.0                knitr_1.42                    
[13] MicrobiomeStat_1.1             Maaslin2_1.10.0               
[15] ALDEx2_1.28.1                  zCompositions_1.4.0-1         
[17] truncnorm_1.0-8                NADA_1.6-1.1                  
[19] survival_3.5-3                 MASS_7.3-58.2                 
[21] tidySummarizedExperiment_1.6.1 patchwork_1.1.2               
[23] mia_1.7.7                      MultiAssayExperiment_1.24.0   
[25] TreeSummarizedExperiment_2.1.4 Biostrings_2.66.0             
[27] XVector_0.38.0                 SingleCellExperiment_1.20.0   
[29] SummarizedExperiment_1.28.0    Biobase_2.58.0                
[31] GenomicRanges_1.50.2           GenomeInfoDb_1.34.9           
[33] IRanges_2.32.0                 S4Vectors_0.36.1              
[35] BiocGenerics_0.44.0            MatrixGenerics_1.10.0         
[37] matrixStats_0.63.0-9003        BiocStyle_2.24.0              
[39] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] estimability_1.4.1          coda_0.19-4                
  [3] bit64_4.0.5                 multcomp_1.4-22            
  [5] irlba_2.3.5.1               DelayedArray_0.24.0        
  [7] data.table_1.14.8           rpart_4.1.19               
  [9] doParallel_1.0.17           RCurl_1.98-1.10            
 [11] generics_0.1.3              ScaledMatrix_1.6.0         
 [13] TH.data_1.1-1               timeSeries_4021.105        
 [15] RSQLite_2.3.0               proxy_0.4-27               
 [17] bit_4.0.5                   tzdb_0.3.0                 
 [19] DirichletMultinomial_1.40.0 viridis_0.6.2              
 [21] xfun_0.37                   fBasics_4021.93            
 [23] hms_1.1.2                   evaluate_0.20              
 [25] DEoptimR_1.0-11             fansi_1.0.4                
 [27] readxl_1.4.2                igraph_1.4.1               
 [29] DBI_1.1.3                   htmlwidgets_1.6.1          
 [31] Rmpfr_0.9-1                 CVXR_1.0-11                
 [33] ellipsis_0.3.2              energy_1.7-11              
 [35] backports_1.4.1             bookdown_0.32              
 [37] permute_0.9-7               deldir_1.0-6               
 [39] sparseMatrixStats_1.10.0    vctrs_0.5.2                
 [41] cachem_1.0.7                withr_2.5.0                
 [43] robustbase_0.95-0           emmeans_1.8.4-1            
 [45] checkmate_2.1.0             vegan_2.6-4                
 [47] treeio_1.22.0               getopt_1.20.3              
 [49] cluster_2.1.4               gsl_2.1-8                  
 [51] ape_5.7                     dir.expiry_1.4.0           
 [53] lazyeval_0.2.2              crayon_1.5.2               
 [55] pkgconfig_2.0.3             nlme_3.1-162               
 [57] vipor_0.4.5                 nnet_7.3-18                
 [59] rlang_1.0.6                 spatial_7.3-16             
 [61] lifecycle_1.0.3             sandwich_3.0-2             
 [63] filelock_1.0.2              phyloseq_1.40.0            
 [65] rsvd_1.0.5                  cellranger_1.1.0           
 [67] rngtools_1.5.2              graph_1.74.0               
 [69] Matrix_1.5-3                lpsymphony_1.24.0          
 [71] zoo_1.8-11                  Rhdf5lib_1.18.2            
 [73] boot_1.3-28.1               base64enc_0.1-3            
 [75] beeswarm_0.4.0              png_0.1-8                  
 [77] viridisLite_0.4.1           stabledist_0.7-1           
 [79] rootSolve_1.8.2.3           bitops_1.0-7               
 [81] rhdf5filters_1.8.0          blob_1.2.3                 
 [83] DelayedMatrixStats_1.20.0   doRNG_1.8.6                
 [85] decontam_1.18.0             jpeg_0.1-10                
 [87] DECIPHER_2.26.0             beachmat_2.14.0            
 [89] scales_1.2.1                memoise_2.0.1              
 [91] magrittr_2.0.3              plyr_1.8.8                 
 [93] zlibbioc_1.44.0             compiler_4.2.1             
 [95] RColorBrewer_1.1-3          clue_0.3-64                
 [97] lme4_1.1-31                 cli_3.6.0                  
 [99] ade4_1.7-22                 lmerTest_3.1-3             
[101] htmlTable_2.4.1             Formula_1.2-5              
[103] mgcv_1.8-41                 tidyselect_1.2.0           
[105] stringi_1.7.12              yaml_2.3.7                 
[107] BiocSingular_1.14.0         latticeExtra_0.6-30        
[109] ggrepel_0.9.3               grid_4.2.1                 
[111] tools_4.2.1                 lmom_2.9                   
[113] timechange_0.2.0            parallel_4.2.1             
[115] rstudioapi_0.14             foreign_0.8-84             
[117] foreach_1.5.2               statip_0.2.3               
[119] optparse_1.7.3              gridExtra_2.3              
[121] gld_2.6.6                   stable_1.1.6               
[123] RcppZiggurat_0.1.6          digest_0.6.31              
[125] BiocManager_1.30.20         Rcpp_1.0.10                
[127] scuttle_1.8.4               httr_1.4.5                 
[129] Rdpack_2.4                  colorspace_2.1-0           
[131] XML_3.99-0.13               modeest_2.4.0              
[133] splines_4.2.1               yulab.utils_0.0.6          
[135] rmutil_1.1.10               statmod_1.5.0              
[137] tidytree_0.4.2              expm_0.999-7               
[139] scater_1.26.1               multtest_2.52.0            
[141] Exact_3.2                   plotly_4.10.1              
[143] xtable_1.8-4                gmp_0.7-1                  
[145] jsonlite_1.8.4              nloptr_2.0.3               
[147] CodeDepends_0.6.5           timeDate_4022.108          
[149] Rfast_2.0.7                 R6_2.5.1                   
[151] Hmisc_4.8-0                 pillar_1.8.1               
[153] htmltools_0.5.4             glue_1.6.2                 
[155] fastmap_1.1.1               minqa_1.2.5                
[157] BiocParallel_1.32.5         BiocNeighbors_1.16.0       
[159] class_7.3-21                codetools_0.2-19           
[161] pcaPP_2.0-3                 mvtnorm_1.1-3              
[163] utf8_1.2.3                  lattice_0.20-45            
[165] numDeriv_2016.8-1.1         ggbeeswarm_0.7.1           
[167] DescTools_0.99.48           interp_1.1-3               
[169] biglm_0.9-2.1               rmarkdown_2.20             
[171] biomformat_1.24.0           munsell_0.5.0              
[173] e1071_1.7-13                rhdf5_2.40.0               
[175] GenomeInfoDbData_1.2.9      iterators_1.0.14           
[177] reshape2_1.4.4              gtable_0.3.1               
[179] rbibutils_2.2.13           
```
</div>

