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
- [ALDEx2](https://bioconductor.org/packages/release/bioc/html/ALDEx2.html)   
- [ANCOM-BC](https://bioconductor.org/packages/release/bioc/html/ANCOMBC.html)  
- [corncob](https://cran.r-project.org/web/packages/corncob/index.html)  
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)  
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)  
- [LEFse](https://bioconductor.org/packages/release/bioc/html/lefser.html)  
- [limma voom](https://bioconductor.org/packages/release/bioc/html/limma.html)  
- [LinDA](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5)  
- [MaAsLin2](https://www.bioconductor.org/packages/release/bioc/html/Maaslin2.html)  
- [metagenomeSeq](https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html)  
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

# For ANCOMBC we need development version
# Can be installed with:
# remotes::install_github("FrederickHuangLin/ANCOMBC")
library(ANCOMBC)

# we use the dmn_se dataset and restrict it to 
# obese vs lean for easy illustration
data(dmn_se)
se <- dmn_se
# To enable all features and advantages of TreeSE, we convert the object from SE to TreeSE
tse <- as(se, "TreeSummarizedExperiment")
tse <- tse[ ,colData(tse)$pheno != "Overwt"]
colData(tse)$pheno <- fct_drop(colData(tse)$pheno, "Overwt")
# how many observations do we have per group?
count(as.data.frame(colData(tse)), pheno) %>% kable()
```


\begin{tabular}{l|r}
\hline
pheno & n\\
\hline
Lean & 61\\
\hline
Obese & 193\\
\hline
\end{tabular}

```r
# set a seed because some tools can randomly vary and then produce 
# different results:
set.seed(1)
```

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
# Convert each instance using the centred log-ratio transform.
# This is the input for all further analyses.
x <- aldex.clr(
  reads = assay(tse),
  conds = colData(tse)$pheno, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)
# calculates expected values of the Welch's t-test and Wilcoxon rank test on
# the data returned by aldex.clr
x_tt <- aldex.ttest(
  x, 
  paired.test = FALSE, 
  verbose = FALSE)
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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-3-1.pdf)<!-- --> 

The evaluation as differential abundant in above plots is based on the
corrected pvalue. According to the ALDEx2 developers, the safest approach is to
identify those features where the 95% CI of the
effect size does not cross 0. As we can see in below table, this is not the
case for any of the identified 
genera (see overlap column, which indicates the proportion of overlap). Also,
the authors recommend to focus on effect sizes and CIs rather than interpreting
the pvalue. To keep the comparison simple, we will here use the pvalue as 
decision criterion. But please be aware that the effect size together with the
CI is a better answer to the question we are typically interested in
(see also [this article](https://www.nature.com/articles/d41586-019-00857-9)).



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
Alistipes & 0.0009 & 0.0001 & -0.3823 & 0.2979\\
\hline
Barnesiella & 0.0442 & 0.0066 & -0.3229 & 0.3489\\
\hline
Catenibacterium & 0.0266 & 0.0330 & 0.2713 & 0.3718\\
\hline
Lactobacillus & 0.0282 & 0.0183 & 0.2983 & 0.3537\\
\hline
Megasphaera & 0.0000 & 0.0001 & 0.5249 & 0.2758\\
\hline
Oscillibacter & 0.0004 & 0.0014 & -0.3681 & 0.3291\\
\hline
Parabacteroides & 0.0541 & 0.0133 & -0.2832 & 0.3509\\
\hline
Phascolarctobacterium & 0.0238 & 0.0077 & -0.3491 & 0.3404\\
\hline
Uknown & 0.0786 & 0.0439 & -0.2474 & 0.3852\\
\hline
\end{tabular}



### ANCOM-BC

The analysis of composition of microbiomes with bias correction 
(ANCOM-BC) [@Huang2020] 
is a recently developed method for differential abundance testing. It is based
on an 
earlier published approach [@Mandal2015]. 
The previous version of ANCOM was among the methods that produced the 
most consistent results and is probably a conservative approach
[@Nearing2022]. 
However, the new ANCOM-BC method operates quite differently compared to the
former ANCOM method.


As the only method, ANCOM-BC incorporates the so called *sampling fraction*
into the model. The latter term could be empirically estimated by the ratio of
the library size to the microbial load. According to the authors, variations in
this sampling fraction would bias differential abundance analyses if ignored.
Furthermore, this method provides p-values and confidence intervals for each
taxon. It also controls the FDR and it is computationally simple to implement. 


As we will see below, to obtain results, all that is needed is to pass 
a tse object to the `ancombc()` function. Below, we first specify a formula.
In this formula, other covariates could potentially be included to adjust for
confounding. We show this further below. 
Please check the
[function documentation](https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc.html) 
to learn about the additional arguments that we specify below.


```r
# perform the analysis 
out <- ancombc(
  x = tse, 
  formula = "pheno", 
  p_adj_method = "fdr", 
  prv_cut = 0, # no prev filtering necessary anymore 
  lib_cut = 0, 
  group = "pheno", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
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
kable(head(res$diff_abn))
```


\begin{tabular}{l|l}
\hline
  & phenoObese\\
\hline
Acetanaerobacterium & TRUE\\
\hline
Acetivibrio & FALSE\\
\hline
Acidaminococcus & TRUE\\
\hline
Akkermansia & FALSE\\
\hline
Alistipes & TRUE\\
\hline
Allisonella & FALSE\\
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
  fixed_effects = "pheno",
  # random_effects = c(...), # you can also fit MLM by specifying random effects
  # specifying a ref is especially important if you have more than 2 levels
  reference = "pheno,Lean",  
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
Megasphaera & pheno & Obese & 0.0489 & 0.0093 & 0 & phenoObese & 0e+00 & 254 & 78\\
\hline
Barnesiella & pheno & Obese & -0.0297 & 0.0068 & 0 & phenoObese & 2e-04 & 254 & 111\\
\hline
Parabacteroides & pheno & Obese & -0.0219 & 0.0050 & 0 & phenoObese & 2e-04 & 254 & 163\\
\hline
Phascolarctobacterium & pheno & Obese & -0.0325 & 0.0072 & 0 & phenoObese & 2e-04 & 254 & 99\\
\hline
Alistipes & pheno & Obese & -0.0523 & 0.0123 & 0 & phenoObese & 3e-04 & 254 & 227\\
\hline
Desulfovibrio & pheno & Obese & -0.0134 & 0.0032 & 0 & phenoObese & 3e-04 & 254 & 72\\
\hline
\end{tabular}

```r
# A folder will be created that is called like the above specified output.
# It contains also figures to visualize the difference between genera 
# for the significant ones.
```

### LinDA 

Lastly, we cover linear models for differential abundance analysis of
microbiome compositional data (@Zhou2022). This tool is very
similar to ANCOMBC with few differences: 1) LinDA correct for the compositional
bias differently 
using the mode of all regression coefficients. 2) The authors claim that it 
runs 100-1000x faster than ANCOMBC and 3) it support hierarchical models.
The latter could be ignored as ANCOMBC will be supporting hierarchical models
with the next release. Nevertheless, LinDA seems a promising tool that achieves 
the best power/fdr trade-off together with ANCOMBC according to the authors. The
speed might make it the choice for bigger datasets or datasets with a very high
number of features.



```r
otu.tab <- as.data.frame(assay(tse))
meta <- as.data.frame(colData(tse)) %>% select(pheno)
res <- linda(
  otu.tab, 
  meta, 
  formula = '~pheno', 
  alpha = 0.05, 
  prev.filter = 0, 
  mean.abund.filter = 0)
```

```
## 0  features are filtered!
## The filtered data has  254  samples and  55  features will be tested!
## Imputation approach is used.
## Fit linear models ...
## Completed.
```

```r
# to scan the table for genera where H0 could be rejected:
kable(head(filter(as.data.frame(res$output), phenoObese.reject)))
```


\begin{tabular}{l|r|r|r|r|r|r|l|r}
\hline
  & phenoObese.baseMean & phenoObese.log2FoldChange & phenoObese.lfcSE & phenoObese.stat & phenoObese.pvalue & phenoObese.padj & phenoObese.reject & phenoObese.df\\
\hline
Acetanaerobacterium & 918.5 & -0.7591 & 0.1619 & -4.687 & 0.0000 & 0.0000 & TRUE & 252\\
\hline
Acidaminococcus & 372.9 & 0.7394 & 0.2327 & 3.178 & 0.0017 & 0.0051 & TRUE & 252\\
\hline
Akkermansia & 827.6 & -0.5110 & 0.2168 & -2.357 & 0.0192 & 0.0422 & TRUE & 252\\
\hline
Alistipes & 23944.9 & -1.7646 & 0.3127 & -5.644 & 0.0000 & 0.0000 & TRUE & 252\\
\hline
Asaccharobacter & 621.0 & -0.4913 & 0.1467 & -3.350 & 0.0009 & 0.0039 & TRUE & 252\\
\hline
Bacteroides & 272640.7 & -0.8874 & 0.2686 & -3.304 & 0.0011 & 0.0040 & TRUE & 252\\
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
summ <- full_join(
    rownames_to_column(aldex_out, "genus") %>%
      select(genus, aldex2 = wi.eBH),
    rownames_to_column(out$res$diff_abn, "genus") %>%
      select(genus, ancombc = phenoObese),
    by = "genus") %>%
  full_join(
    select(fit_data$results, genus = feature, maaslin2 = qval), 
    by = "genus") %>%
    full_join(
      rownames_to_column(as.data.frame(res$output), "genus") %>%
        select(genus, LinDA = phenoObese.reject), 
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
Acetanaerobacterium & FALSE & TRUE & TRUE & TRUE & 3\\
\hline
Acetivibrio & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
Acidaminococcus & FALSE & TRUE & TRUE & TRUE & 3\\
\hline
Akkermansia & FALSE & FALSE & FALSE & TRUE & 1\\
\hline
Alistipes & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Allisonella & FALSE & FALSE & FALSE & FALSE & 0\\
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
9 & 22 & 16 & 25\\
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
Alistipes & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Barnesiella & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Catenibacterium & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Lactobacillus & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Megasphaera & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Oscillibacter & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Parabacteroides & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
Phascolarctobacterium & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
\end{tabular}

We see that each method identified at least 9 genera as differentially
abundant. Eight of those that were identified by ALDEx2,
were also identified by the other methods. We could plot the data for
any method or for those taxa that were identified by all methods:




```r
plot_data <- data.frame(t(assay(tse)))
plot_data$pheno <- colData(tse)$pheno
# create a plot for each genus where the score is indicated in the title
plots <- pmap(select(summ, genus, score), function(genus, score) {
  ggplot(plot_data, aes_string("pheno", genus)) +
    geom_boxplot(aes(fill = pheno), outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ggtitle(glue::glue("Score {score}")) +
    theme_bw() +
    theme(legend.position = "none")
})

# now we can show only those genera that have at least score 3 (or 2 or 1)
robust_plots <- plots[summ$score == 4] 

# to display this nicely in the book we use patchwork here:
# (we show first 8)
robust_plots[[1]] + 
  robust_plots[[2]] + 
  robust_plots[[3]] + 
  robust_plots[[4]] +
  robust_plots[[5]] +
  robust_plots[[6]] +
  robust_plots[[7]] +
  robust_plots[[8]] +
  plot_layout(nrow = 2)
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-12-1.pdf)<!-- --> 

```r
# or if we have most trust in any specific method we can show genera that 
# are differentially abundant according to that method and then look in the
# title how many methods also identified it (we only show first 6 here):
ancombc_plots <- plots[summ$ancombc] 
ancombc_plots[[1]] + 
  ancombc_plots[[2]] + 
  ancombc_plots[[3]] + 
  ancombc_plots[[4]] +
  ancombc_plots[[5]] +
  ancombc_plots[[6]] 
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-12-2.pdf)<!-- --> 



### Confounding variables

To perform causal inference, it is crucial that the method is able to include
covariates in the model. This is not possible with e.g. the Wilcoxon test.
Other methods such as both ANCOM methods, ALDEx2, LinDA, MaAsLin2 and others
allow this. Below we show how to include a covariate in ANCOM-BC.
It is very similar for all the methods that allow this. Since in this dataset
there are no covariates, I first simulate a new variable and add it to the TSE
object.



```r
# to join new data to existing colData we need to put rownames as a column 
colData(tse)$sample_id <- rownames(colData(tse))
# simulate a covariate that I will add to the colData.
df_sim <- tibble(
  sample_id = colData(tse)$sample_id,
  age = rnorm(n = length(colData(tse)$sample_id))
)
# an easy way to join data is to use dplyr functions. The package 
# tidySummarizedExperiment enables this functionality
tse <- full_join(tse, df_sim, by = "sample_id")
out_cov = ancombc(
  x = tse, 
  formula = "pheno + age", # here we add age to the model
  p_adj_method = "fdr", 
  prv_cut = 0,  # we did that already
  lib_cut = 0, 
  group = "pheno", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
# now the model answers the question: holding age constant, are 
# bacterial taxa differentially abundant? Or, if that is of interest,
# holding phenotype constant, is age associated with bacterial abundance?
# Again we only show the first 6 entries.
kable(head(out_cov$res$diff_abn))
```


\begin{tabular}{l|l|l}
\hline
  & phenoObese & age\\
\hline
Acetanaerobacterium & TRUE & FALSE\\
\hline
Acetivibrio & FALSE & FALSE\\
\hline
Acidaminococcus & TRUE & FALSE\\
\hline
Akkermansia & FALSE & FALSE\\
\hline
Alistipes & TRUE & FALSE\\
\hline
Allisonella & FALSE & FALSE\\
\hline
\end{tabular}

In the next section of this book chapter we cover methods that can also take
into account the phylogenetic information of bacterial taxa to perform 
group-wise associations.









## Tree-based methods

### Group-wise associations testing based on balances with fido

[TreeSummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/TreeSummarizedExperiment.html) 
frequently includes a Phylogenetic tree along with associated data about the
experiment (at `colData`), that holds covariates which can be used for
analyzing group-wise associations. 

Such an analysis could be performed with the function `pibble` from the `fido`
package, that offers a Multinomial Logistic-Normal Linear Regression model; see
[vignette](https://jsilve24.github.io/fido/articles/introduction-to-fido.html) of package.


The following presents such an exemplary analysis based on the 
data of @Sprockett2020 available
through `microbiomeDataSets` package.



```r
if (!require(fido)){
  # installing the fido package
  devtools::install_github("jsilve24/fido")
}
```

Loading the libraries and importing data:


```r
library(fido)
library(mia)
```


```r
library(microbiomeDataSets)
tse <- SprockettTHData()
```




We pick three covariates ("Sex","Age_Years","Delivery_Mode") during this
analysis as an example, and beforehand we check for missing data:



```r
cov_names <- c("Sex","Age_Years","Delivery_Mode")
na_counts <- apply(is.na(colData(tse)[,cov_names]), 2, sum)
na_summary<-as.data.frame(na_counts,row.names=cov_names)
```

We drop samples with na values at the covariates (features) under analysis: 



```r
tse <- tse[ , !is.na(colData(tse)$Delivery_Mode) ]
tse <- tse[ , !is.na(colData(tse)$Age_Years) ]
```

We agglomerate the data at a Phylum rank.
Note: Large assay data (along with the covariates/features data) could prevent the analysis later,
since the computation will construct matrices that would not always fit memory.


```r
tse_phylum <- agglomerateByRank(tse, "Phylum")
```

We extract the counts assay and feature data to build the model matrix having
an extra row of ones presenting the intercept for the regression task later: 


```r
Y <- assays(tse_phylum)$counts
# design matrix
# taking 3 covariates
sample_data<-as.data.frame(colData(tse_phylum)[,cov_names])
X <- t(model.matrix(~Sex+Age_Years+Delivery_Mode,data=sample_data))
```

Building the parameters for the `pibble` call to build the model; see more at [vignette](https://jsilve24.github.io/fido/articles/introduction-to-fido.html):


```r
n_taxa<-nrow(Y)
upsilon <- n_taxa+3
Omega <- diag(n_taxa)
G <- cbind(diag(n_taxa-1), -1)
Xi <- (upsilon-n_taxa)*G%*%Omega%*%t(G)
Theta <- matrix(0, n_taxa-1, nrow(X))
Gamma <- diag(nrow(X))
```

Automatically initializing the priors and visualizing their distributions:


```r
priors <- pibble(NULL, X, upsilon, Theta, Gamma, Xi)
names_covariates(priors) <- rownames(X)
plot(priors, pars="Lambda") + ggplot2::xlim(c(-5, 5))
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-23-1.pdf)<!-- --> 

Estimating the posterior by including the data at `Y`.
Note: Some computational failures could occur (see [discussion](https://github-wiki-see.page/m/jsilve24/fido/wiki/Frequently-Asked-Questions))
the arguments `multDirichletBoot` `calcGradHess` could be passed in such case.


```r
priors$Y <- Y 
posterior <- refit(priors, optim_method="adam", multDirichletBoot=0.5) # ,, calcGradHess=FALSE
```

Printing a summary about the posterior predictive distribution:


```r
ppc_summary(posterior)
```

```
## Proportions of Observations within 95% Credible Interval: 0.998
```
Plotting the summary of the posterior distributions of the regression parameters:


```r
names_categories(posterior) <- rownames(Y)
plot(posterior,par="Lambda",focus.cov=rownames(X)[2:4])
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-26-1.pdf)<!-- --> 

Seemingly the covariate "Age_Years" does not have effect on the model as "Delivery_Mode" would,
and "Sex" to some extent. Let's take a closer look at the two latter ones:


```r
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-27-1.pdf)<!-- --> 

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
 [1] fido_1.0.3                     ANCOMBC_1.99.1                
 [3] forcats_0.5.2                  stringr_1.4.1                 
 [5] dplyr_1.0.10                   purrr_0.3.4                   
 [7] readr_2.1.3                    tidyr_1.2.1                   
 [9] tibble_3.1.8                   ggplot2_3.3.6                 
[11] tidyverse_1.3.2                knitr_1.40                    
[13] MicrobiomeStat_1.1             Maaslin2_1.10.0               
[15] ALDEx2_1.28.1                  zCompositions_1.4.0-1         
[17] truncnorm_1.0-8                NADA_1.6-1.1                  
[19] survival_3.4-0                 MASS_7.3-58.1                 
[21] tidySummarizedExperiment_1.6.1 patchwork_1.1.2               
[23] mia_1.5.16                     MultiAssayExperiment_1.22.0   
[25] TreeSummarizedExperiment_2.1.4 Biostrings_2.64.1             
[27] XVector_0.36.0                 SingleCellExperiment_1.18.1   
[29] SummarizedExperiment_1.26.1    Biobase_2.56.0                
[31] GenomicRanges_1.48.0           GenomeInfoDb_1.32.4           
[33] IRanges_2.30.1                 S4Vectors_0.34.0              
[35] BiocGenerics_0.42.0            MatrixGenerics_1.8.1          
[37] matrixStats_0.62.0-9003        BiocStyle_2.24.0              
[39] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] coda_0.19-4                 bit64_4.0.5                
  [3] irlba_2.3.5.1               DelayedArray_0.22.0        
  [5] data.table_1.14.2           rpart_4.1.16               
  [7] doParallel_1.0.17           RCurl_1.98-1.9             
  [9] generics_0.1.3              ScaledMatrix_1.4.1         
 [11] microbiome_1.18.0           timeSeries_4021.104        
 [13] RSQLite_2.2.17              proxy_0.4-27               
 [15] bit_4.0.4                   tzdb_0.3.0                 
 [17] xml2_1.3.3                  lubridate_1.8.0            
 [19] assertthat_0.2.1            DirichletMultinomial_1.38.0
 [21] viridis_0.6.2               gargle_1.2.1               
 [23] xfun_0.33                   fBasics_4021.92            
 [25] ggdist_3.2.0                hms_1.1.2                  
 [27] evaluate_0.16               DEoptimR_1.0-11            
 [29] fansi_1.0.3                 dbplyr_2.2.1               
 [31] readxl_1.4.1                igraph_1.3.5               
 [33] DBI_1.1.3                   htmlwidgets_1.5.4          
 [35] tensorA_0.36.2              googledrive_2.0.0          
 [37] hash_2.2.6.2                ellipsis_0.3.2             
 [39] energy_1.7-10               backports_1.4.1            
 [41] bookdown_0.29               permute_0.9-7              
 [43] deldir_1.0-6                sparseMatrixStats_1.8.0    
 [45] vctrs_0.4.2                 abind_1.4-5                
 [47] tidybayes_3.0.2             cachem_1.0.6               
 [49] withr_2.5.0                 robustbase_0.95-0          
 [51] checkmate_2.1.0             vegan_2.6-2                
 [53] treeio_1.20.2               getopt_1.20.3              
 [55] cluster_2.1.4               gsl_2.1-7.1                
 [57] ape_5.6-2                   dir.expiry_1.4.0           
 [59] lazyeval_0.2.2              crayon_1.5.2               
 [61] labeling_0.4.2              pkgconfig_2.0.3            
 [63] nlme_3.1-159                vipor_0.4.5                
 [65] nnet_7.3-18                 rlang_1.0.6                
 [67] spatial_7.3-15              lifecycle_1.0.2            
 [69] filelock_1.0.2              phyloseq_1.40.0            
 [71] modelr_0.1.9                rsvd_1.0.5                 
 [73] distributional_0.3.1        cellranger_1.1.0           
 [75] rngtools_1.5.2              graph_1.74.0               
 [77] Matrix_1.5-1                lpsymphony_1.24.0          
 [79] Rhdf5lib_1.18.2             boot_1.3-28                
 [81] base64enc_0.1-3             reprex_2.0.2               
 [83] beeswarm_0.4.0              googlesheets4_1.0.1        
 [85] png_0.1-7                   viridisLite_0.4.1          
 [87] stabledist_0.7-1            rootSolve_1.8.2.3          
 [89] bitops_1.0-7                rhdf5filters_1.8.0         
 [91] blob_1.2.3                  DelayedMatrixStats_1.18.1  
 [93] doRNG_1.8.2                 decontam_1.16.0            
 [95] jpeg_0.1-9                  DECIPHER_2.24.0            
 [97] beachmat_2.12.0             scales_1.2.1               
 [99] memoise_2.0.1               magrittr_2.0.3             
[101] plyr_1.8.7                  zlibbioc_1.42.0            
[103] compiler_4.2.1              RColorBrewer_1.1-3         
[105] clue_0.3-61                 lme4_1.1-30                
[107] cli_3.4.1                   ade4_1.7-19                
[109] lmerTest_3.1-3              pbapply_1.5-0              
[111] htmlTable_2.4.1             Formula_1.2-4              
[113] mgcv_1.8-40                 tidyselect_1.1.2           
[115] stringi_1.7.8               highr_0.9                  
[117] yaml_2.3.5                  svUnit_1.0.6               
[119] BiocSingular_1.12.0         latticeExtra_0.6-30        
[121] ggrepel_0.9.1               grid_4.2.1                 
[123] tools_4.2.1                 lmom_2.9                   
[125] parallel_4.2.1              rstudioapi_0.14            
[127] logging_0.10-108            foreign_0.8-83             
[129] foreach_1.5.2               statip_0.2.3               
[131] optparse_1.7.3              gridExtra_2.3              
[133] gld_2.6.5                   posterior_1.3.1            
[135] farver_2.1.1                Rtsne_0.16                 
[137] stable_1.1.6                RcppZiggurat_0.1.6         
[139] digest_0.6.29               BiocManager_1.30.18        
[141] Rcpp_1.0.9                  broom_1.0.1                
[143] scuttle_1.6.3               httr_1.4.4                 
[145] Rdpack_2.4                  colorspace_2.0-3           
[147] rvest_1.0.3                 XML_3.99-0.11              
[149] fs_1.5.2                    modeest_2.4.0              
[151] splines_4.2.1               yulab.utils_0.0.5          
[153] rmutil_1.1.9                statmod_1.4.37             
[155] expm_0.999-6                tidytree_0.4.1             
[157] scater_1.24.0               Exact_3.2                  
[159] multtest_2.52.0             plotly_4.10.0              
[161] jsonlite_1.8.2              nloptr_2.0.3               
[163] CodeDepends_0.6.5           timeDate_4021.106          
[165] Rfast_2.0.6                 R6_2.5.1                   
[167] Hmisc_4.7-1                 pillar_1.8.1               
[169] htmltools_0.5.3             glue_1.6.2                 
[171] fastmap_1.1.0               minqa_1.2.4                
[173] BiocParallel_1.30.3         BiocNeighbors_1.14.0       
[175] class_7.3-20                codetools_0.2-18           
[177] pcaPP_2.0-2                 mvtnorm_1.1-3              
[179] utf8_1.2.2                  lattice_0.20-45            
[181] arrayhelpers_1.1-0          numDeriv_2016.8-1.1        
[183] ggbeeswarm_0.6.0            DescTools_0.99.46          
[185] interp_1.1-3                biglm_0.9-2.1              
[187] rmarkdown_2.16              biomformat_1.24.0          
[189] munsell_0.5.0               e1071_1.7-11               
[191] rhdf5_2.40.0                GenomeInfoDbData_1.2.8     
[193] iterators_1.0.14            haven_2.5.1                
[195] reshape2_1.4.4              gtable_0.3.1               
[197] rbibutils_2.2.9            
```
</div>
