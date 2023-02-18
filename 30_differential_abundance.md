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
x <- aldex.clr(
  reads = assay(tse),
  conds = tse$Geographical_location
)
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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-2-1.pdf)<!-- --> 

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


\begin{tabular}{l|r|r|r|r}
\hline
genus & we.eBH & wi.eBH & effect & overlap\\
\hline
OTU194 & 0.0570 & 0.0156 & 0.8731 & 0.1634\\
\hline
OTU562 & 0.0959 & 0.0255 & -0.7378 & 0.1660\\
\hline
OTU773 & 0.0317 & 0.0032 & 1.1905 & 0.1092\\
\hline
OTU860 & 0.0839 & 0.0353 & -0.8652 & 0.1761\\
\hline
OTU1075 & 0.0410 & 0.0083 & -1.1272 & 0.1321\\
\hline
OTU1235 & 0.0331 & 0.0311 & -0.8827 & 0.1548\\
\hline
OTU1680 & 0.0883 & 0.0381 & -0.9757 & 0.1787\\
\hline
OTU2529 & 0.0994 & 0.0458 & -0.8732 & 0.1816\\
\hline
\end{tabular}

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


\begin{tabular}{l|r|r|r|r|r|r|r|r|r|r|l|l}
\hline
taxon & lfc\_(Intercept) & lfc\_Geographical\_locationPune & se\_(Intercept) & se\_Geographical\_locationPune & W\_(Intercept) & W\_Geographical\_locationPune & p\_(Intercept) & p\_Geographical\_locationPune & q\_(Intercept) & q\_Geographical\_locationPune & diff\_(Intercept) & diff\_Geographical\_locationPune\\
\hline
OTU2 & 0.0399 & -0.0570 & 0.1675 & 0.1915 & 0.2383 & -0.2975 & 0.8117 & 0.7661 & 0.8718 & 0.8463 & FALSE & FALSE\\
\hline
OTU15 & 0.6874 & -0.9024 & 0.1947 & 0.2225 & 3.5304 & -4.0547 & 0.0004 & 0.0001 & 0.0032 & 0.0004 & TRUE & TRUE\\
\hline
OTU53 & 0.1243 & -0.1672 & 0.7823 & 0.8940 & 0.1589 & -0.1870 & 0.8737 & 0.8516 & 0.8969 & 0.8701 & FALSE & FALSE\\
\hline
OTU87 & 0.1347 & -0.1807 & 0.1938 & 0.2215 & 0.6952 & -0.8161 & 0.4869 & 0.4145 & 0.6596 & 0.5616 & FALSE & FALSE\\
\hline
OTU99 & 0.2716 & -0.3594 & 0.1635 & 0.1869 & 1.6608 & -1.9231 & 0.0967 & 0.0545 & 0.2793 & 0.1504 & FALSE & FALSE\\
\hline
OTU111 & 0.0237 & -0.0358 & 0.1677 & 0.1917 & 0.1413 & -0.1868 & 0.8876 & 0.8519 & 0.8994 & 0.8701 & FALSE & FALSE\\
\hline
\end{tabular}



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


```r
# which genera are identified as differentially abundant? (leave out "head" to
# see all)
kable(head(filter(fit_data$results, qval <= 0.05)))
```


\begin{tabular}{l|l|l|r|r|r|l|r|r|r}
\hline
feature & metadata & value & coef & stderr & pval & name & qval & N & N.not.zero\\
\hline
OTU1053 & Geographical\_location & Pune & -0.0080 & 0.0011 & 0 & Geographical\_locationPune & 0 & 47 & 9\\
\hline
OTU860 & Geographical\_location & Pune & -0.0373 & 0.0059 & 0 & Geographical\_locationPune & 0 & 47 & 13\\
\hline
OTU1075 & Geographical\_location & Pune & -0.1295 & 0.0207 & 0 & Geographical\_locationPune & 0 & 47 & 27\\
\hline
OTU1980 & Geographical\_location & Pune & -0.0395 & 0.0062 & 0 & Geographical\_locationPune & 0 & 47 & 9\\
\hline
OTU611 & Geographical\_location & Pune & -0.0274 & 0.0045 & 0 & Geographical\_locationPune & 0 & 47 & 10\\
\hline
OTU2335 & Geographical\_location & Pune & -0.0089 & 0.0015 & 0 & Geographical\_locationPune & 0 & 47 & 10\\
\hline
\end{tabular}

```r
# A folder will be created that is called like the above specified output.
# It contains also figures to visualize the difference between genera 
# for the significant ones.
```

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
```

```
## 0  features are filtered!
## The filtered data has  47  samples and  262  features will be tested!
## Pseudo-count approach is used.
## Fit linear models ...
## Completed.
```

```r
# to scan the table for genera where H0 could be rejected:
kable(head(filter(as.data.frame(res$output), Geographical_locationPune.reject)))
```


\begin{tabular}{l|r|r|r|r|r|r|l|r}
\hline
  & Geographical\_locationPune.baseMean & Geographical\_locationPune.log2FoldChange & Geographical\_locationPune.lfcSE & Geographical\_locationPune.stat & Geographical\_locationPune.pvalue & Geographical\_locationPune.padj & Geographical\_locationPune.reject & Geographical\_locationPune.df\\
\hline
OTU15 & 1194.9 & -1.9113 & 0.3579 & -5.340 & 0.0000 & 0.0000 & TRUE & 45\\
\hline
OTU22 & 393.5 & -0.6683 & 0.2184 & -3.060 & 0.0037 & 0.0160 & TRUE & 45\\
\hline
OTU76 & 837.0 & -1.8013 & 0.3598 & -5.006 & 0.0000 & 0.0001 & TRUE & 45\\
\hline
OTU127 & 938.2 & -1.7558 & 0.3912 & -4.488 & 0.0000 & 0.0004 & TRUE & 45\\
\hline
OTU170 & 416.9 & -0.7518 & 0.2351 & -3.197 & 0.0025 & 0.0123 & TRUE & 45\\
\hline
OTU194 & 869.3 & 4.3671 & 1.2254 & 3.564 & 0.0009 & 0.0054 & TRUE & 45\\
\hline
\end{tabular}


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


\begin{tabular}{l|l|l|l|l|r}
\hline
genus & aldex2 & ancombc & maaslin2 & LinDA & score\\
\hline
OTU2 & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
OTU15 & FALSE & TRUE & TRUE & TRUE & 3\\
\hline
OTU22 & FALSE & NA & TRUE & TRUE & NA\\
\hline
OTU53 & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
OTU69 & FALSE & NA & FALSE & FALSE & NA\\
\hline
OTU76 & FALSE & NA & TRUE & TRUE & NA\\
\hline
\end{tabular}

Now we can answer our questions:


```r
# how many genera were identified by each method?
summarise(summ, across(where(is.logical), sum)) %>%
  kable()
```


\begin{tabular}{r|r|r|r}
\hline
aldex2 & ancombc & maaslin2 & LinDA\\
\hline
8 & NA & 67 & 75\\
\hline
\end{tabular}

```r
# which genera are identified by all methods?
filter(summ, score == 4) %>% kable()
```


\begin{tabular}{l|l|l|l|l|r}
\hline
genus & aldex2 & ancombc & maaslin2 & LinDA & score\\
\hline
OTU773 & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
OTU860 & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
OTU1075 & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
OTU1235 & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
OTU1680 & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
OTU2529 & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
\end{tabular}

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
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> 

```r
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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-11-2.pdf)<!-- --> 



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


\begin{tabular}{l|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|r|l|l|l|l}
\hline
taxon & lfc\_(Intercept) & lfc\_Geographical\_locationPune & lfc\_AgeElderly & lfc\_AgeMiddle\_age & se\_(Intercept) & se\_Geographical\_locationPune & se\_AgeElderly & se\_AgeMiddle\_age & W\_(Intercept) & W\_Geographical\_locationPune & W\_AgeElderly & W\_AgeMiddle\_age & p\_(Intercept) & p\_Geographical\_locationPune & p\_AgeElderly & p\_AgeMiddle\_age & q\_(Intercept) & q\_Geographical\_locationPune & q\_AgeElderly & q\_AgeMiddle\_age & diff\_(Intercept) & diff\_Geographical\_locationPune & diff\_AgeElderly & diff\_AgeMiddle\_age\\
\hline
OTU2 & 0.0397 & -0.0948 & 0.0893 & 0.0126 & 0.1725 & 0.2388 & 0.2299 & 0.2344 & 0.2304 & -0.3971 & 0.3884 & 0.0537 & 0.8178 & 0.6913 & 0.6977 & 0.9572 & 0.8853 & 0.9266 & 0.9125 & 0.9952 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU15 & 0.7028 & -0.7558 & -0.2434 & -0.1577 & 0.1987 & 0.2751 & 0.2648 & 0.2700 & 3.5364 & -2.7476 & -0.9190 & -0.5839 & 0.0004 & 0.0060 & 0.3581 & 0.5593 & 0.0031 & 0.0387 & 0.8827 & 0.9860 & TRUE & TRUE & FALSE & FALSE\\
\hline
OTU53 & 0.0205 & -1.0766 & 1.4992 & 1.1534 & 0.7871 & 1.0895 & 1.0494 & 1.0701 & 0.0260 & -0.9881 & 1.4287 & 1.0778 & 0.9793 & 0.3231 & 0.1531 & 0.2811 & 0.9878 & 0.6462 & 0.7721 & 0.7870 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU87 & 0.1765 & 0.1238 & -0.4584 & -0.4487 & 0.1906 & 0.2638 & 0.2540 & 0.2590 & 0.9259 & 0.4691 & -1.8046 & -1.7322 & 0.3545 & 0.6390 & 0.0711 & 0.0832 & 0.5502 & 0.8903 & 0.5713 & 0.6117 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU99 & 0.3029 & -0.1602 & -0.2768 & -0.3341 & 0.1636 & 0.2264 & 0.2179 & 0.2222 & 1.8521 & -0.7074 & -1.2703 & -1.5035 & 0.0640 & 0.4793 & 0.2040 & 0.1327 & 0.2171 & 0.7723 & 0.8133 & 0.6117 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU111 & 0.0023 & -0.1442 & 0.1352 & 0.2461 & 0.1710 & 0.2367 & 0.2279 & 0.2324 & 0.0135 & -0.6093 & 0.5935 & 1.0592 & 0.9892 & 0.5424 & 0.5528 & 0.2895 & 0.9892 & 0.8310 & 0.8900 & 0.7870 & FALSE & FALSE & FALSE & FALSE\\
\hline
\end{tabular}

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
 [1] ANCOMBC_2.1.2                  forcats_1.0.0                 
 [3] stringr_1.5.0                  dplyr_1.1.0                   
 [5] purrr_1.0.1                    readr_2.1.4                   
 [7] tidyr_1.3.0                    tibble_3.1.8                  
 [9] ggplot2_3.4.1                  tidyverse_1.3.2               
[11] knitr_1.42                     MicrobiomeStat_1.1            
[13] Maaslin2_1.10.0                ALDEx2_1.28.1                 
[15] zCompositions_1.4.0-1          truncnorm_1.0-8               
[17] NADA_1.6-1.1                   survival_3.5-3                
[19] MASS_7.3-58.2                  tidySummarizedExperiment_1.6.1
[21] patchwork_1.1.2                mia_1.7.5                     
[23] MultiAssayExperiment_1.24.0    TreeSummarizedExperiment_2.1.4
[25] Biostrings_2.66.0              XVector_0.38.0                
[27] SingleCellExperiment_1.20.0    SummarizedExperiment_1.28.0   
[29] Biobase_2.58.0                 GenomicRanges_1.50.2          
[31] GenomeInfoDb_1.34.9            IRanges_2.32.0                
[33] S4Vectors_0.36.1               BiocGenerics_0.44.0           
[35] MatrixGenerics_1.10.0          matrixStats_0.63.0-9003       
[37] BiocStyle_2.24.0               rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] estimability_1.4.1          coda_0.19-4                
  [3] bit64_4.0.5                 irlba_2.3.5.1              
  [5] DelayedArray_0.24.0         data.table_1.14.6          
  [7] rpart_4.1.19                doParallel_1.0.17          
  [9] RCurl_1.98-1.10             generics_0.1.3             
 [11] ScaledMatrix_1.6.0          timeSeries_4021.105        
 [13] RSQLite_2.2.20              proxy_0.4-27               
 [15] bit_4.0.5                   tzdb_0.3.0                 
 [17] xml2_1.3.3                  lubridate_1.9.2            
 [19] assertthat_0.2.1            DirichletMultinomial_1.40.0
 [21] viridis_0.6.2               gargle_1.3.0               
 [23] xfun_0.37                   fBasics_4021.93            
 [25] hms_1.1.2                   evaluate_0.20              
 [27] DEoptimR_1.0-11             fansi_1.0.4                
 [29] dbplyr_2.3.0                readxl_1.4.2               
 [31] igraph_1.4.0                DBI_1.1.3                  
 [33] htmlwidgets_1.6.1           googledrive_2.0.0          
 [35] hash_2.2.6.2                Rmpfr_0.9-1                
 [37] CVXR_1.0-11                 ellipsis_0.3.2             
 [39] energy_1.7-11               backports_1.4.1            
 [41] bookdown_0.32               permute_0.9-7              
 [43] deldir_1.0-6                sparseMatrixStats_1.10.0   
 [45] vctrs_0.5.2                 cachem_1.0.6               
 [47] withr_2.5.0                 robustbase_0.95-0          
 [49] emmeans_1.8.4-1             checkmate_2.1.0            
 [51] vegan_2.6-4                 treeio_1.22.0              
 [53] getopt_1.20.3               cluster_2.1.4              
 [55] gsl_2.1-8                   ape_5.7                    
 [57] dir.expiry_1.4.0            lazyeval_0.2.2             
 [59] crayon_1.5.2                labeling_0.4.2             
 [61] pkgconfig_2.0.3             nlme_3.1-162               
 [63] vipor_0.4.5                 nnet_7.3-18                
 [65] rlang_1.0.6                 spatial_7.3-16             
 [67] lifecycle_1.0.3             filelock_1.0.2             
 [69] phyloseq_1.40.0             modelr_0.1.10              
 [71] rsvd_1.0.5                  cellranger_1.1.0           
 [73] rngtools_1.5.2              graph_1.74.0               
 [75] Matrix_1.5-3                lpsymphony_1.24.0          
 [77] Rhdf5lib_1.18.2             boot_1.3-28.1              
 [79] base64enc_0.1-3             reprex_2.0.2               
 [81] beeswarm_0.4.0              googlesheets4_1.0.1        
 [83] png_0.1-8                   viridisLite_0.4.1          
 [85] stabledist_0.7-1            rootSolve_1.8.2.3          
 [87] bitops_1.0-7                rhdf5filters_1.8.0         
 [89] blob_1.2.3                  DelayedMatrixStats_1.20.0  
 [91] doRNG_1.8.6                 decontam_1.18.0            
 [93] jpeg_0.1-10                 DECIPHER_2.26.0            
 [95] beachmat_2.14.0             scales_1.2.1               
 [97] memoise_2.0.1               magrittr_2.0.3             
 [99] plyr_1.8.8                  zlibbioc_1.44.0            
[101] compiler_4.2.1              RColorBrewer_1.1-3         
[103] clue_0.3-64                 lme4_1.1-31                
[105] cli_3.6.0                   ade4_1.7-22                
[107] lmerTest_3.1-3              pbapply_1.7-0              
[109] htmlTable_2.4.1             Formula_1.2-4              
[111] mgcv_1.8-41                 tidyselect_1.2.0           
[113] stringi_1.7.12              highr_0.10                 
[115] yaml_2.3.7                  BiocSingular_1.14.0        
[117] latticeExtra_0.6-30         ggrepel_0.9.3              
[119] grid_4.2.1                  lmom_2.9                   
[121] tools_4.2.1                 timechange_0.2.0           
[123] parallel_4.2.1              rstudioapi_0.14            
[125] logging_0.10-108            foreign_0.8-84             
[127] foreach_1.5.2               statip_0.2.3               
[129] optparse_1.7.3              gridExtra_2.3              
[131] gld_2.6.6                   farver_2.1.1               
[133] stable_1.1.6                RcppZiggurat_0.1.6         
[135] digest_0.6.31               BiocManager_1.30.19        
[137] Rcpp_1.0.10                 broom_1.0.3                
[139] scuttle_1.8.4               httr_1.4.4                 
[141] Rdpack_2.4                  colorspace_2.1-0           
[143] rvest_1.0.3                 XML_3.99-0.13              
[145] fs_1.6.1                    modeest_2.4.0              
[147] splines_4.2.1               yulab.utils_0.0.6          
[149] rmutil_1.1.10               statmod_1.5.0              
[151] expm_0.999-7                tidytree_0.4.2             
[153] scater_1.26.1               Exact_3.2                  
[155] multtest_2.52.0             plotly_4.10.1              
[157] xtable_1.8-4                gmp_0.7-1                  
[159] jsonlite_1.8.4              nloptr_2.0.3               
[161] CodeDepends_0.6.5           timeDate_4022.108          
[163] Rfast_2.0.7                 R6_2.5.1                   
[165] Hmisc_4.8-0                 pillar_1.8.1               
[167] htmltools_0.5.4             glue_1.6.2                 
[169] fastmap_1.1.0               minqa_1.2.5                
[171] BiocParallel_1.32.5         BiocNeighbors_1.16.0       
[173] class_7.3-21                codetools_0.2-19           
[175] pcaPP_2.0-3                 mvtnorm_1.1-3              
[177] utf8_1.2.3                  lattice_0.20-45            
[179] numDeriv_2016.8-1.1         ggbeeswarm_0.7.1           
[181] DescTools_0.99.47           interp_1.1-3               
[183] biglm_0.9-2.1               rmarkdown_2.20             
[185] biomformat_1.24.0           munsell_0.5.0              
[187] e1071_1.7-13                rhdf5_2.40.0               
[189] GenomeInfoDbData_1.2.9      iterators_1.0.14           
[191] haven_2.5.1                 reshape2_1.4.4             
[193] gtable_0.3.1                rbibutils_2.2.13           
```
</div>

