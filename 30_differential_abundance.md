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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-2-1.pdf)<!-- --> 

The evaluation as differential abundant in above plots is based on the
corrected p-value. According to the ALDEx2 developers, the safest
approach is to identify those features where the 95% CI of the effect
size does not cross 0. As we can see in below table, this is not the
case for any of the identified genera (see overlap column, which
indicates the proportion of overlap). Also, the authors recommend to
focus on effect sizes and CIs rather than interpreting the p-value. To
keep the comparison simple, we will here use the p-value as decision
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
OTU194 & 0.0553 & 0.0151 & 0.9552 & 0.1563\\
\hline
OTU562 & 0.0714 & 0.0266 & -0.7147 & 0.1648\\
\hline
OTU611 & 0.1204 & 0.0468 & -0.7630 & 0.1745\\
\hline
OTU773 & 0.0281 & 0.0036 & 1.1930 & 0.1037\\
\hline
OTU860 & 0.0866 & 0.0380 & -0.8501 & 0.1733\\
\hline
OTU1075 & 0.0374 & 0.0083 & -1.1059 & 0.1278\\
\hline
OTU1235 & 0.0280 & 0.0253 & -0.9094 & 0.1702\\
\hline
OTU1680 & 0.0834 & 0.0341 & -0.9270 & 0.1915\\
\hline
OTU2529 & 0.1113 & 0.0449 & -0.8634 & 0.1929\\
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

Note that the original method was implemented in the `ancombc()` function (see 
[extended tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html)).
The method has since then been updated and new features have been added to enable
multi-group comparisons and repeated measurements among other improvements. 
We do not cover the more advanced features of ANCOMBC in this tutorial 
as these features are documented in detail in this 
[tutorial](https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html).

We now proceed with a simple example.  First, we specify a formula. In this 
formula, other covariates could potentially be included to adjust for 
confounding. We show this further below. Again, please make sure to check the 
[function documentation](https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc.html)
as well as the linked tutorials to learn about the additional arguments 
that we specify.



```r
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
Abyssicoccus & 0.0399 & -0.0570 & 0.1675 & 0.1915 & 0.2383 & -0.2975 & 0.8117 & 0.7661 & 0.8718 & 0.8463 & FALSE & FALSE\\
\hline
Acidaminococcus & 0.6874 & -0.9024 & 0.1947 & 0.2225 & 3.5304 & -4.0547 & 0.0004 & 0.0001 & 0.0032 & 0.0004 & TRUE & TRUE\\
\hline
Acinetobacter & 0.1243 & -0.1672 & 0.7823 & 0.8940 & 0.1589 & -0.1870 & 0.8737 & 0.8516 & 0.8969 & 0.8701 & FALSE & FALSE\\
\hline
Actinomyces & 0.1347 & -0.1807 & 0.1938 & 0.2215 & 0.6952 & -0.8161 & 0.4869 & 0.4145 & 0.6596 & 0.5616 & FALSE & FALSE\\
\hline
Actinoplanes & 0.2716 & -0.3594 & 0.1635 & 0.1869 & 1.6608 & -1.9231 & 0.0967 & 0.0545 & 0.2793 & 0.1504 & FALSE & FALSE\\
\hline
Aerococcus & 0.0237 & -0.0358 & 0.1677 & 0.1917 & 0.1413 & -0.1868 & 0.8876 & 0.8519 & 0.8994 & 0.8701 & FALSE & FALSE\\
\hline
\end{tabular}



### MaAsLin2 

Next, we will illustrate the use of MaAsLin2. The method is based on
generalized linear models. It is flexible for different study designs
and covariate structures. For more details, you can check their [official
tutorial](https://github.com/biobakery/biobakery/wiki/maaslin2).



```r
# maaslin expects features as columns and samples as rows 
# for both the asv/otu table as well as meta data 
asv <- t(assay(tse))
meta_data <- data.frame(colData(tse))
# We can specify different GLMs/normalizations/transforms.
# Let us use similar settings as in Nearing et al. (2021):
fit_data <- Maaslin2(
  asv,
  meta_data,
  output = "DAA example",
  transform = "AST",
  fixed_effects = "Geographical_location",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  # reference = "Geographical_location,Pune",  
  normalization = "TSS",
  standardize = FALSE,
  min_prevalence = 0 # prev filterin already done
)
```

Which genera are identified as differentially abundant? (leave out "head" to see all).


```r
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
kable(head(filter(as.data.frame(res$output$Geographical_locationPune), reject)))
```


\begin{tabular}{l|r|r|r|r|r|r|l|r}
\hline
  & baseMean & log2FoldChange & lfcSE & stat & pvalue & padj & reject & df\\
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
```


\begin{tabular}{l|l|l|l|l|r}
\hline
taxid & aldex2 & ancombc & maaslin2 & LinDA & score\\
\hline
OTU2 & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
OTU15 & FALSE & TRUE & TRUE & TRUE & 3\\
\hline
OTU22 & FALSE & FALSE & TRUE & TRUE & 2\\
\hline
OTU53 & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
OTU69 & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
OTU76 & FALSE & FALSE & TRUE & TRUE & 2\\
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
9 & 31 & 67 & 75\\
\hline
\end{tabular}

```r
# which genera are identified by all methods?
filter(summ, score == 4) %>% kable()
```


\begin{tabular}{l|l|l|l|l|r}
\hline
taxid & aldex2 & ancombc & maaslin2 & LinDA & score\\
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

We see that each method identified at least some genera as differentially
abundant. Many of those that were identified by ALDEx2,
were also identified by the other methods. Let us plot the data for
any method or for those taxa that were identified by all methods:



```r
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
```

![](30_differential_abundance_files/figure-latex/daplotting-1.pdf)<!-- --> 

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

![](30_differential_abundance_files/figure-latex/daplotting-2.pdf)<!-- --> 



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
OTU2 & 0.0441 & -0.1005 & 0.0892 & 0.0118 & 0.1727 & 0.2390 & 0.2301 & 0.2346 & 0.2553 & -0.4203 & 0.3876 & 0.0504 & 0.7985 & 0.6742 & 0.6983 & 0.9598 & 0.8739 & 0.9143 & 0.9128 & 0.9925 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU15 & 0.7071 & -0.7615 & -0.2435 & -0.1584 & 0.1988 & 0.2752 & 0.2649 & 0.2702 & 3.5567 & -2.7669 & -0.9190 & -0.5865 & 0.0004 & 0.0057 & 0.3581 & 0.5576 & 0.0029 & 0.0365 & 0.8827 & 0.9891 & TRUE & TRUE & FALSE & FALSE\\
\hline
OTU53 & 0.0248 & -1.0822 & 1.4991 & 1.1526 & 0.7869 & 1.0893 & 1.0491 & 1.0698 & 0.0315 & -0.9935 & 1.4289 & 1.0774 & 0.9749 & 0.3205 & 0.1530 & 0.2813 & 0.9749 & 0.6300 & 0.7718 & 0.7889 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU87 & 0.1808 & 0.1182 & -0.4585 & -0.4494 & 0.1907 & 0.2640 & 0.2541 & 0.2592 & 0.9481 & 0.4476 & -1.8039 & -1.7342 & 0.3431 & 0.6544 & 0.0712 & 0.0829 & 0.5236 & 0.9143 & 0.5707 & 0.6097 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU99 & 0.3073 & -0.1658 & -0.2769 & -0.3349 & 0.1637 & 0.2266 & 0.2182 & 0.2225 & 1.8767 & -0.7314 & -1.2693 & -1.5053 & 0.0606 & 0.4645 & 0.2043 & 0.1323 & 0.2096 & 0.7479 & 0.8240 & 0.6097 & FALSE & FALSE & FALSE & FALSE\\
\hline
OTU111 & 0.0066 & -0.1498 & 0.1351 & 0.2453 & 0.1712 & 0.2369 & 0.2281 & 0.2326 & 0.0388 & -0.6325 & 0.5925 & 1.0548 & 0.9691 & 0.5271 & 0.5535 & 0.2915 & 0.9749 & 0.8026 & 0.8901 & 0.7889 & FALSE & FALSE & FALSE & FALSE\\
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
R version 4.3.0 (2023-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] doRNG_1.8.6                     rngtools_1.5.2                 
 [3] foreach_1.5.2                   ANCOMBC_2.2.0                  
 [5] lubridate_1.9.2                 forcats_1.0.0                  
 [7] stringr_1.5.0                   dplyr_1.1.2                    
 [9] purrr_1.0.1                     readr_2.1.4                    
[11] tidyr_1.3.0                     tibble_3.2.1                   
[13] ggplot2_3.4.2                   tidyverse_2.0.0                
[15] knitr_1.42                      MicrobiomeStat_1.1             
[17] Maaslin2_1.7.3                  ALDEx2_1.32.0                  
[19] zCompositions_1.4.0-1           truncnorm_1.0-9                
[21] NADA_1.6-1.1                    survival_3.5-5                 
[23] MASS_7.3-59                     tidySummarizedExperiment_1.10.0
[25] patchwork_1.1.2                 mia_1.7.11                     
[27] MultiAssayExperiment_1.26.0     TreeSummarizedExperiment_2.1.4 
[29] Biostrings_2.68.0               XVector_0.40.0                 
[31] SingleCellExperiment_1.22.0     SummarizedExperiment_1.30.0    
[33] Biobase_2.60.0                  GenomicRanges_1.52.0           
[35] GenomeInfoDb_1.36.0             IRanges_2.34.0                 
[37] S4Vectors_0.38.0                BiocGenerics_0.46.0            
[39] MatrixGenerics_1.12.0           matrixStats_0.63.0-9003        
[41] BiocStyle_2.28.0                rebook_1.9.0                   

loaded via a namespace (and not attached):
  [1] bitops_1.0-7                DirichletMultinomial_1.42.0
  [3] RColorBrewer_1.1-3          doParallel_1.0.17          
  [5] httr_1.4.5                  numDeriv_2016.8-1.1        
  [7] backports_1.4.1             tools_4.3.0                
  [9] utf8_1.2.3                  R6_2.5.1                   
 [11] vegan_2.6-4                 lazyeval_0.2.2             
 [13] mgcv_1.8-42                 rhdf5filters_1.12.0        
 [15] permute_0.9-7               withr_2.5.0                
 [17] gridExtra_2.3               cli_3.6.1                  
 [19] logging_0.10-108            biglm_0.9-2.1              
 [21] sandwich_3.0-2              labeling_0.4.2             
 [23] mvtnorm_1.1-3               robustbase_0.95-1          
 [25] pbapply_1.7-0               proxy_0.4-27               
 [27] yulab.utils_0.0.6           foreign_0.8-84             
 [29] scater_1.28.0               decontam_1.20.0            
 [31] readxl_1.4.2                rstudioapi_0.14            
 [33] RSQLite_2.3.1               generics_0.1.3             
 [35] Matrix_1.5-4                biomformat_1.28.0          
 [37] ggbeeswarm_0.7.1            fansi_1.0.4                
 [39] DescTools_0.99.48           DECIPHER_2.28.0            
 [41] lifecycle_1.0.3             multcomp_1.4-23            
 [43] yaml_2.3.7                  rhdf5_2.44.0               
 [45] grid_4.3.0                  blob_1.2.4                 
 [47] crayon_1.5.2                dir.expiry_1.8.0           
 [49] lattice_0.21-8              beachmat_2.16.0            
 [51] CodeDepends_0.6.5           pillar_1.9.0               
 [53] optparse_1.7.3              statip_0.2.3               
 [55] boot_1.3-28.1               gld_2.6.6                  
 [57] estimability_1.4.1          codetools_0.2-19           
 [59] glue_1.6.2                  data.table_1.14.8          
 [61] Rdpack_2.4                  vctrs_0.6.2                
 [63] treeio_1.24.0               cellranger_1.1.0           
 [65] gtable_0.3.3                cachem_1.0.7               
 [67] xfun_0.39                   rbibutils_2.2.13           
 [69] Rfast_2.0.7                 coda_0.19-4                
 [71] pcaPP_2.0-3                 modeest_2.4.0              
 [73] timeDate_4022.108           iterators_1.0.14           
 [75] statmod_1.5.0               gmp_0.7-1                  
 [77] TH.data_1.1-2               ellipsis_0.3.2             
 [79] nlme_3.1-162                phyloseq_1.44.0            
 [81] bit64_4.0.5                 filelock_1.0.2             
 [83] fBasics_4022.94             irlba_2.3.5.1              
 [85] vipor_0.4.5                 rpart_4.1.19               
 [87] Hmisc_5.0-1                 colorspace_2.1-0           
 [89] DBI_1.1.3                   nnet_7.3-18                
 [91] ade4_1.7-22                 Exact_3.2                  
 [93] tidyselect_1.2.0            emmeans_1.8.5              
 [95] timeSeries_4021.105         bit_4.0.5                  
 [97] compiler_4.3.0              graph_1.78.0               
 [99] htmlTable_2.4.1             BiocNeighbors_1.18.0       
[101] expm_0.999-7                DelayedArray_0.25.0        
[103] plotly_4.10.1               bookdown_0.33              
[105] checkmate_2.2.0             scales_1.2.1               
[107] DEoptimR_1.0-12             spatial_7.3-16             
[109] digest_0.6.31               minqa_1.2.5                
[111] rmarkdown_2.21              base64enc_0.1-3            
[113] htmltools_0.5.5             pkgconfig_2.0.3            
[115] lme4_1.1-33                 sparseMatrixStats_1.12.0   
[117] lpsymphony_1.28.0           highr_0.10                 
[119] stabledist_0.7-1            fastmap_1.1.1              
[121] rlang_1.1.0                 htmlwidgets_1.6.2          
[123] DelayedMatrixStats_1.22.0   farver_2.1.1               
[125] energy_1.7-11               zoo_1.8-12                 
[127] jsonlite_1.8.4              BiocParallel_1.34.0        
[129] BiocSingular_1.16.0         RCurl_1.98-1.12            
[131] magrittr_2.0.3              Formula_1.2-5              
[133] scuttle_1.10.0              GenomeInfoDbData_1.2.10    
[135] Rhdf5lib_1.22.0             munsell_0.5.0              
[137] Rcpp_1.0.10                 ape_5.7-1                  
[139] viridis_0.6.2               RcppZiggurat_0.1.6         
[141] CVXR_1.0-11                 stringi_1.7.12             
[143] rootSolve_1.8.2.3           stable_1.1.6               
[145] zlibbioc_1.46.0             plyr_1.8.8                 
[147] parallel_4.3.0              ggrepel_0.9.3              
[149] lmom_2.9                    splines_4.3.0              
[151] hash_2.2.6.2                multtest_2.56.0            
[153] hms_1.1.3                   igraph_1.4.2               
[155] reshape2_1.4.4              ScaledMatrix_1.7.1         
[157] rmutil_1.1.10               XML_3.99-0.14              
[159] evaluate_0.20               BiocManager_1.30.20        
[161] nloptr_2.0.3                tzdb_0.3.0                 
[163] getopt_1.20.3               clue_0.3-64                
[165] rsvd_1.0.5                  xtable_1.8-4               
[167] Rmpfr_0.9-2                 e1071_1.7-13               
[169] tidytree_0.4.2              viridisLite_0.4.1          
[171] class_7.3-21                gsl_2.1-8                  
[173] lmerTest_3.1-3              memoise_2.0.1              
[175] beeswarm_0.4.0              cluster_2.1.4              
[177] timechange_0.2.0           
```
</div>


