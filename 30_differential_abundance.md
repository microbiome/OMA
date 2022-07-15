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
# now the data from df_sim is in the tse object and we can again repeat
# the steps as above:
assay(tse)
```

```
##                       TS1.2 TS10.2 TS100.2 TS100 TS101.2 TS103.2 TS103 TS104
## Acetanaerobacterium       0      0       0     1       0       0     0     0
## Acetivibrio               0      0       0     0       0       0     0     0
## Acidaminococcus           0      0       1     0       0       0     0     0
## Akkermansia               1      0       1     1       1       0     0     0
## Alistipes                41      0       3    23       0      30     1    51
## Allisonella               0      0       2     0       1       0     0     0
## Anaerostipes              0      0       0     0       0       0     0     0
## Anaerotruncus            37      9       0     0       9       1     0     1
## Asaccharobacter           0      0       0     0       0       0     0     1
## Bacteroides             194    227      56   555      34     552   124   415
## Barnesiella              21      0       0     5       0       5     1     0
## Bifidobacterium           0      0       0     0       0       0     0     0
## Butyricicoccus            1      4       0     2       0       4     1     3
## Butyricimonas             1      0       0     2       0       0     0     0
## Catenibacterium           0      0     139   126       0       0     0     0
## Clostridium               1      0       1     3       1       0     1     0
## Collinsella               0      7      10     2       0       3     1     0
## Coprobacillus             7     32       2     3       0       0    11     0
## Coprococcus               3      2      42    32       2      11    20    17
## Desulfovibrio             3      0       0     0       0       0     0     0
## Dialister                 9      8       0     0       0       0     0     0
## Dorea                     3     17      87    86      81      42    39    23
## Eggerthella               0      2       1     1       3       0     0     1
## Eubacterium               2      0      63   110     148       0     0     0
## Faecalibacterium         83     93      89   126      34     261    33   131
## Gordonibacter             0      0       0     0       1       0     0     0
## Hespellia                 0      0      34     0       1      16     3     3
## Holdemania                0      0       0     0       0       0     0     1
## Lachnobacterium           0      0       0     0       0       0     0     0
## Lactobacillus             0      0       0     4       2       0     0     0
## Lactococcus               0      0       0     0       0       0     0     0
## Lactonifactor             0      0      13    18       0      13     2     1
## Marvinbryantia            0      0       1     1       0       0     0     0
## Megasphaera               0      0      11    10       0       0     0     0
## Odoribacter               0      0       0     5       0       7     0     1
## Oribacterium              0      0       0     0       0       0     0     0
## Oscillibacter           113      0       3    12       1      29     2    37
## Parabacteroides           5      1       1     0       0       2     2     0
## Paraprevotella            0      0       0     9       0       0     0     0
## Parasutterella            1      0       0     0       0       0     3     0
## Phascolarctobacterium     0      0       0     0       2       9     6    10
## Prevotella                0      0       0     0       1       0     0     0
## Pseudobutyrivibrio        6      5       0     0       7       0     1     0
## Roseburia                 8     42      19    58      19      17    87    19
## Ruminococcus             12      0      24    25       3       0     0    10
## Slackia                   0      0       0     0       0       0     0     0
## Sporacetigenium           0      0      54    25       0       0     0     0
## Sporobacter               0      0       0     0       0       0     0     0
## Streptococcus             0      1       9     4       7       3     2     1
## Subdoligranulum          52      8      13    33      34      25     5    10
## Sutterella                0      0       1     6       0      54    15     0
## Tepidibacter              0      4       8     0       0       0     0     0
## Turicibacter              0      0       0     0       0       1     0     0
## Uknown                  206    196     429   365     564     408   248   287
## Veillonella               0      1       0     0       1       0     2     0
##                       TS105.2 TS105 TS106.2 TS106 TS107.2 TS107 TS109.2 TS109
## Acetanaerobacterium         1     1       0     0       0     0       1     1
## Acetivibrio                 0     0       0     0       0     0       0     0
## Acidaminococcus             0     0       0     0       0    11       0     0
## Akkermansia                 0     0       0     0       0     1       2     0
## Alistipes                   9     2       2     4      11    26      17    21
## Allisonella                 0     0       0     0       0     0       0     0
## Anaerostipes                0     0       0     0       0     0       0     0
## Anaerotruncus               0     0       2     0       0     0       9     0
## Asaccharobacter             0     0       1     1       1     1       0     0
## Bacteroides               792    60     453   275     112   292     375   709
## Barnesiella                 0     0       0     0       0     0       9     6
## Bifidobacterium             0     0       1     1       0     0       0     0
## Butyricicoccus              0     0       1     4       0     1       4     3
## Butyricimonas               0     0       0     0       0     0       4     8
## Catenibacterium             0     0       0     0     184     6       0     0
## Clostridium                 0     0       0     3       1     5       0     1
## Collinsella                 0     0      23     5     105    37       0     1
## Coprobacillus               0     0       1    21       5     0       9    25
## Coprococcus                 7     9      17    19       5     2       7    10
## Desulfovibrio               0     0       0     0       0     0      10     4
## Dialister                  10     7       1     0      63    31       0     0
## Dorea                       0     0      18    41     118    58      19     8
## Eggerthella                 0     1       1     0       0     2       0     0
## Eubacterium                 2     0       0     0     124    51       4     0
## Faecalibacterium          138    18     149   115     327    76     145   143
## Gordonibacter               0     0       0     0       0     0       1     0
## Hespellia                   0    31       8     1       0     2       0     0
## Holdemania                  0     0       0     0       0     0       0     1
## Lachnobacterium             0     0       0     0       0     0       5     0
## Lactobacillus               0     0       0     0       0     0       0     0
## Lactococcus                 0     2       0     0       0     0       0     0
## Lactonifactor               0     0       1     0       0     8       0     0
## Marvinbryantia              0     0       0     0       0     0       0     0
## Megasphaera                 0     0       0     0      56    44       0     0
## Odoribacter                 0     0       0     0       0     0       0     2
## Oribacterium                0     0       0     1       0     2       0     0
## Oscillibacter               6     0       2     1       3     2      11    12
## Parabacteroides             0     0       1     2       1     5       1    16
## Paraprevotella              0     0       0     0       0     0       0     1
## Parasutterella              0     0       0     0       0     0      12     2
## Phascolarctobacterium       0     0       0     0       0     0      22    31
## Prevotella                  0     0       0     0      11    72       5    15
## Pseudobutyrivibrio          0     0       1     0       0     0       0     0
## Roseburia                  50    29      22    78      50    34      51    35
## Ruminococcus                3     9       0    12       0     0      55    25
## Slackia                     0     0       0     0       2     4       0     0
## Sporacetigenium             0     0       3    33      24     2       0     0
## Sporobacter                 0     0       0     1       0     0       0     0
## Streptococcus               0     1      20    25      11    11       1     1
## Subdoligranulum             4     7       4     0     104    23       1    15
## Sutterella                  0     0       0     0       2    18       0     0
## Tepidibacter                0     0       0     0       0     0       0     0
## Turicibacter                0     0       0     1      20     5       0     0
## Uknown                    157   599     372   475     507   333     340   316
## Veillonella                 0     0       3     7       0     0       0     0
##                       TS10 TS11.2 TS110.2 TS110 TS111.2 TS111 TS115.2 TS115
## Acetanaerobacterium      0      0       4     0       0     2       0     3
## Acetivibrio              0      0       1     1       0     0       0     1
## Acidaminococcus          0      0       0     0       0     0       0     2
## Akkermansia              0      0       0     0       2     2       0     0
## Alistipes                0     16     268    41      39    65       0     7
## Allisonella              0      0       0     0       0     0       0     1
## Anaerostipes             0      1       0     1       0     1       0     0
## Anaerotruncus           27     10       8    10      14     9      24     1
## Asaccharobacter          0      0       2     1       1     5       0     0
## Bacteroides            434    613     367    99     338   342       7    32
## Barnesiella              0      0      13    15       5    26       0     0
## Bifidobacterium          0      0       0     0       0     0       0     0
## Butyricicoccus           0      1       2     1       1     1       1     0
## Butyricimonas            0      0       3     0       1     0       2     2
## Catenibacterium          0      0       0     0       0     0      52    59
## Clostridium              1      0       0     3       0     0       2     2
## Collinsella             25      6       0     2       0     5      33    34
## Coprobacillus          113      2     112    70       5    31       0     0
## Coprococcus              5      8      13    16      22    24      22    41
## Desulfovibrio            0      0       0     0       0     1       0     0
## Dialister                6     20       0     0       0     0      19    27
## Dorea                   46     15       9    15      24    31      13    20
## Eggerthella              0      0       1     2       0     0       0     0
## Eubacterium              1      3       6     0       0     0       8    13
## Faecalibacterium       361    248     484    40     188   123     115   238
## Gordonibacter            0      0       0     0       0     2       0     0
## Hespellia                0      7      18     0       1    23       0     0
## Holdemania               0      1       0     1       2     2       0     0
## Lachnobacterium          6      0       7     0       0     0       0     0
## Lactobacillus            0      0       0     0       0     0       7     2
## Lactococcus              3      0       0     0       0     0       0     0
## Lactonifactor            0      0       0     0       0    11       2     2
## Marvinbryantia           0      3       1     1       0     0       1     0
## Megasphaera              0      0       0     0       0     0       5     4
## Odoribacter              0      0       9     1       2     5       2     2
## Oribacterium            13      8       0     0       0    21       1     0
## Oscillibacter            3      0     151    11      45    68      51    57
## Parabacteroides          1      2      91     1      51    39       0     0
## Paraprevotella           0      0      23     2       0     0       0     2
## Parasutterella           0      0       0     0       8     3       0     0
## Phascolarctobacterium    0      0       0     0      21    45       0     1
## Prevotella               0      0       0     0       6    24      60   107
## Pseudobutyrivibrio       0     39       0     0       0     0       0     2
## Roseburia               82    100      50    10      87   118      16    13
## Ruminococcus             0      0      31    54      85    26       8    23
## Slackia                  0      0       0     0       0     1       0     0
## Sporacetigenium          4      3       0     0       0     0       0     0
## Sporobacter              5      0       1     0       3     5       0     0
## Streptococcus            7      5       1     8      14    39       1     0
## Subdoligranulum         11      1      38   160      29    39       8    16
## Sutterella               0      0       0     0       0     0       1     5
## Tepidibacter             0      0       0     0       0     0       3     0
## Turicibacter             1      0       0     1       2    12       1     4
## Uknown                 447    580    1000   535     645   927     242   254
## Veillonella              3      4       0     0       0     0       0     0
##                       TS116.2 TS116 TS117.2 TS117 TS118.2 TS118 TS119.2 TS119
## Acetanaerobacterium         1     4       0     0       0     0       0     2
## Acetivibrio                 0     0       0     0       0     0       0     0
## Acidaminococcus             8    53       0     0       0     0       0     0
## Akkermansia                 0     0       0     0       0     0       0     0
## Alistipes                   6    14       7    36      14    13       9     2
## Allisonella                 1     1       3     1       0     0       0     0
## Anaerostipes                0     0       0     0       0     0       0     0
## Anaerotruncus               0    63       0     0      36    11       2     1
## Asaccharobacter             0     0       0     0       0     0       1     1
## Bacteroides                21   224     493   359     185   754     509   517
## Barnesiella                 0     0       0     0       1     0       0     0
## Bifidobacterium             0     1       1     2       1     0       1     1
## Butyricicoccus              0     2       2     4       3     2       3     6
## Butyricimonas               6     1       4     3       6     7       8    11
## Catenibacterium            71   119       0     0       0     0       0     0
## Clostridium                 0     1      14    11       1     1       0     0
## Collinsella                 8    51      19    35      53     4      32    32
## Coprobacillus               0     0      52    75       8     6      38    25
## Coprococcus               131    63       5     0       6    10      37    18
## Desulfovibrio               6    12       0     0       0     2       0     0
## Dialister                  29    22       0     0       0     0      11    19
## Dorea                      12    20     105    55      41    26      17    26
## Eggerthella                 0     0       0     1       0     0       0     0
## Eubacterium                29    50       1     6     120    54       1     2
## Faecalibacterium          137   136     247   332     123   295     153    13
## Gordonibacter               0     0       0     0       0     0       0     0
## Hespellia                   1     0       0     3       0     1       5     4
## Holdemania                  0     0       1     0       0     4       1     0
## Lachnobacterium             2    10       0     1       0     0       0     0
## Lactobacillus               2     0       0     0       0     0       0     1
## Lactococcus                 0     0       0     0       0     0       0     0
## Lactonifactor               0     5       0     0       0     0       6     1
## Marvinbryantia              1     1       0     0       0     0       0     0
## Megasphaera                28    35       0     0       8     5      42   204
## Odoribacter                 5    16       0     0       2     4       0     0
## Oribacterium                0     1       0     0       0     4       0     0
## Oscillibacter              59   487       0     5       8    11       7     5
## Parabacteroides             1     2       0     0       0     0       4     5
## Paraprevotella             13     4       0     0       0     0       0     0
## Parasutterella              0     0       0     0       0     6       0     0
## Phascolarctobacterium       0     0       0     0       0     0       0     0
## Prevotella                 46    17       0     0       4     0       0     0
## Pseudobutyrivibrio          0     0       1    25       0     0       0     0
## Roseburia                  27    14      83    76      25    26      68    39
## Ruminococcus                7    26      29    27       0     1      58    67
## Slackia                     1     1       1     0       0     0       0     0
## Sporacetigenium             2     1       4     0       0     0       9     0
## Sporobacter                 0     0       0     0       0     1       0     0
## Streptococcus               1     2       0     3       2     0       1     2
## Subdoligranulum             6    22      28    30      56    14     146   133
## Sutterella                  1     5      87    59       0     0       0     0
## Tepidibacter                0     0       0     0       0     0       0     0
## Turicibacter                0     0       0     3       1     0       0     0
## Uknown                    387  1114     290   311     477   238     657   562
## Veillonella                 0     0       0     0       0     0       0     0
##                       TS11 TS12.2 TS120.2 TS120 TS124.2 TS124 TS126.2 TS126
## Acetanaerobacterium      0      0       3     0       1     6       2     4
## Acetivibrio              0      0       0     0       0     1       0     0
## Acidaminococcus          0      0       6     0       0     0       0     0
## Akkermansia              0      0       0     0       0     0       0     0
## Alistipes               21      0       1     0      73    29      24    12
## Allisonella              0      0       1     0       0     0       0     0
## Anaerostipes             0      0       0     0       0     0       0     0
## Anaerotruncus            2      0      34     2       1     2       3     3
## Asaccharobacter          0      0       0     0       0     0       0     0
## Bacteroides            707    254     130     3     759   507     192    33
## Barnesiella              0      0       2     1       1     1      61     5
## Bifidobacterium          1      0       1     0       0     0       0     0
## Butyricicoccus           2      0       2     1       0     1       1     0
## Butyricimonas            0      0       8     0       0     0       8     4
## Catenibacterium          0      1       0     0       0     0      22    27
## Clostridium              1      0       1     0       0     0       0     1
## Collinsella              9      6      22     4      12     5      17    11
## Coprobacillus            0      0      15     1       7     3       1     0
## Coprococcus             13      6      41     1      20     6      65    30
## Desulfovibrio            2      0       2     0       0     1       1     0
## Dialister                3     19       0     0       0     0       0    14
## Dorea                    0     31      37     1      15    68       6     5
## Eggerthella              0      0       0     0       1     1       0     0
## Eubacterium              0      0     107     1       8    28      18     7
## Faecalibacterium       281    295     178     4     124    29      88    77
## Gordonibacter            0      1       0     0       0     0       0     0
## Hespellia                6      0       0     0       1    34      12     3
## Holdemania               0      0       0     0       1     2       1     0
## Lachnobacterium          0      0       0     0       0     0       3     0
## Lactobacillus            0      0      14     0       1     0     175   100
## Lactococcus              0      0       0     0       0     0       5     1
## Lactonifactor            0      0      22     0       0     2       2     0
## Marvinbryantia           0      0       0     0       0     0       0     0
## Megasphaera              0      0      26     1       0     0      18    10
## Odoribacter              0      0       9     0       0     0       8     1
## Oribacterium             6      0       0     0       0     0       0     0
## Oscillibacter            8      0      53     2      58    42     200    11
## Parabacteroides          0      0       1     0       1     0       1     0
## Paraprevotella           0      8       2     1       0     0       2     0
## Parasutterella           0      0       0     0       0     0       0     1
## Phascolarctobacterium    0      0       0     0      36    35       3     8
## Prevotella               0      0     114     4       0     1       0     1
## Pseudobutyrivibrio       0      0       0     0       1     0       0     0
## Roseburia              105     30      21     0      20    48       9     9
## Ruminococcus             0      0       3     2      36    45      19     9
## Slackia                  0      0       0     0       0     0       0     0
## Sporacetigenium          0      1       0     0       0     0      14     0
## Sporobacter              0      0       0     0       0     0       0     0
## Streptococcus            8     11      13     0       3     1     393   332
## Subdoligranulum          0     44      38     6       6    15      13    39
## Sutterella               1      0      31     2       0     4       0     0
## Tepidibacter             0      2       1     0       0     0       2     0
## Turicibacter             0      0       0     0       0     2       1     2
## Uknown                 460    476     303    16     310   679     238   299
## Veillonella              5      3       0     0       0     0       0     0
##                       TS127.2 TS127 TS128.2 TS128 TS129.2 TS129 TS12 TS13.2
## Acetanaerobacterium         1     2       0     0       0     2    0      0
## Acetivibrio                 0     3       0     4       0     0    2      0
## Acidaminococcus             0     0       0     0       0     0    0      0
## Akkermansia                 0    20       0     0       7     4    1      0
## Alistipes                  25   180      16    46      35    32    0     48
## Allisonella                 0     0       0     0       0     0    0      0
## Anaerostipes                0     1       0     0       0     1    0      0
## Anaerotruncus               1   305       3     6       1     1    0      0
## Asaccharobacter             0     0       0     2       0     1    0      1
## Bacteroides               204   356     308  1121     764   526  524     81
## Barnesiella                30    60       4     8       0     3    0     23
## Bifidobacterium             0     0       0     0       0     0    0      0
## Butyricicoccus              7     4       1     1       2     5    1      1
## Butyricimonas               1     0       0     0       3     6    0      2
## Catenibacterium             0     0       3     0       0     0    0      0
## Clostridium                 4     1       7     3       1     1    2      0
## Collinsella                 4     8      25    63       0    10   19      9
## Coprobacillus               0     1       0     5       9    59    0    104
## Coprococcus                26    66       9    17       6     5   14     54
## Desulfovibrio               0     1       0     0       0     0    0      9
## Dialister                   8    12       2     0       0     0    2      0
## Dorea                       5    46      15    16       3     7   73     24
## Eggerthella                 0     5       0     0       0     0    0      1
## Eubacterium                 0     0       1     3       0     0    2      0
## Faecalibacterium          157    37     338   500     258   119  132     90
## Gordonibacter               0     0       1     0       0     0    0      0
## Hespellia                   0     0       1     0       0     0    4      4
## Holdemania                  1     1       3     2       2     0    0      0
## Lachnobacterium             0     0      18     0      41     0    0      0
## Lactobacillus               0     0       1     0       0     0    0      0
## Lactococcus                 0     0       0     0       1     1    0      0
## Lactonifactor               0     1       0    13       1     0    1      3
## Marvinbryantia              0     0       0     0       0     0    0      0
## Megasphaera                 0     0       0     0       0     0    0      0
## Odoribacter                 4    16       0     4       0     0    0      1
## Oribacterium                3    13       0    13       0     0    1      0
## Oscillibacter              44   134      17    53      46    38    0     31
## Parabacteroides             0     4       0     9       0     2    0      1
## Paraprevotella              0     0       0     0       0     0   14      0
## Parasutterella              3     5       2    15       1     0    0      0
## Phascolarctobacterium       2     0       9    36      17    19    0      5
## Prevotella                  0     0       1     0       0     1    0      3
## Pseudobutyrivibrio          1     0       0     0       0     0    0      0
## Roseburia                  73    34     158    95      49    21   62     15
## Ruminococcus               18    57      97    86      15    11    0     28
## Slackia                     0     0       0     0       0     0    0      1
## Sporacetigenium             1    19      21    28       3     2    0      1
## Sporobacter                 0     0       0     0       0     0    0      0
## Streptococcus               0     1      13    17       1     8   23      5
## Subdoligranulum            32    36     127    61      23     5   82     12
## Sutterella                  0     0       0     0       0     0    0     11
## Tepidibacter                0     1       3     0       0     0    1      0
## Turicibacter                2    15       0    17       0     0    0      0
## Uknown                    568   681     672   476     272   407  413    291
## Veillonella                 0     0       6     0       0     0    3      0
##                       TS130.2 TS130 TS131.2 TS131 TS132.2 TS132 TS133.2 TS133
## Acetanaerobacterium         0     0       1     0       0     0       0     0
## Acetivibrio                 0     0       0     0       0     0       0     0
## Acidaminococcus             0     0       0     0       0     0       0     0
## Akkermansia                 0     0       0     0       0     0       0     0
## Alistipes                   0     0       3     4      31    58       7    27
## Allisonella                 0     0       0     1       0     0       0     0
## Anaerostipes                0     0       1     2       0     0       0     0
## Anaerotruncus               2     1       0     3      43     0       0     6
## Asaccharobacter             0     0       0     0       0     0       0     0
## Bacteroides              1110   900     606   632     913   540     895   336
## Barnesiella                 0     0       0     0       0     0       0     0
## Bifidobacterium             0     0       0     1       1     0       3     0
## Butyricicoccus              3     0       3     4       3     0       6     0
## Butyricimonas               0     0       0     0       0     0       0     0
## Catenibacterium             0     3       0     0       0     0     100    61
## Clostridium                 2     2       1     9       1     0       0     6
## Collinsella                15     0       0     0      51     8      25    20
## Coprobacillus               4     0       0     1      31    14       4     3
## Coprococcus                 6     3       0     0      14     4      28     2
## Desulfovibrio               0     0       0     0       0     0       0     0
## Dialister                   9     0      27    30      20     0      16     5
## Dorea                      17     2      92   134      53    24      45    32
## Eggerthella                 0     0       5     5       1     2       0     0
## Eubacterium                 2     8       0     4       0     1       0     0
## Faecalibacterium          203   200      12    39      83    84     315   229
## Gordonibacter               0     0       0     0       0     0       0     0
## Hespellia                  17    20      27    15       2    38       2     2
## Holdemania                  0     1       0     0       1     1       2     0
## Lachnobacterium             0     0       0     0       0     0       0     0
## Lactobacillus               0     0      17    40       1     0       0     0
## Lactococcus                 0     0       0     0       0     0       0     0
## Lactonifactor               0    13       0     7       0     0       2     3
## Marvinbryantia              0     0       0     0       0     0       0     0
## Megasphaera                 0     0       0     0      39     0      28    39
## Odoribacter                 0     0       1     0       0    11       0     0
## Oribacterium                0     0       0     0       0     0       0     1
## Oscillibacter              11     2       4     1     127    49      79    26
## Parabacteroides             0     0       4     1       2     7       0     1
## Paraprevotella              0     0       0     0       0     0       0     0
## Parasutterella              0     0       0     0      10     5      17    13
## Phascolarctobacterium       0     0       0     0       1    23       0     0
## Prevotella                  0     2       0     0       0     0       0     0
## Pseudobutyrivibrio          0     0       0     1       0     0       0   121
## Roseburia                  56    62      16    96      63    30      33    42
## Ruminococcus                0    12       0     0      10     9      13    17
## Slackia                     0     0       0     0       4     0       0     0
## Sporacetigenium             0     0       4     8      40     0       1     0
## Sporobacter                 8     0       0     0       0     0       1     0
## Streptococcus               1     0       2     2      18     2       3    12
## Subdoligranulum            22     3       0     1      44    23      58    78
## Sutterella                 12     0       0     0       0     0       0     0
## Tepidibacter                1     0       0     0       7     0       1     0
## Turicibacter                0     0       0     0       2     0       0     0
## Uknown                    273   300     895   614     539   409     353   120
## Veillonella                 1     0       0     0       0     2       0     0
##                       TS134.2 TS134 TS135.2 TS135 TS136.2 TS137.2 TS137 TS138.2
## Acetanaerobacterium         0     1       1     0       2       1     3       0
## Acetivibrio                 0     0       0     0       0       0     0       0
## Acidaminococcus             0     0       5    14       0       0     0       4
## Akkermansia                 0     0       1     1       4       0     0       0
## Alistipes                  63    17      52    43      40      11    26       5
## Allisonella                 0     0       0     0       0       0     4       0
## Anaerostipes                0     1       0     0       0       0     0       0
## Anaerotruncus               0     0       1     4       3       1     4       4
## Asaccharobacter             0     1       0     1       0       0     0       0
## Bacteroides               314   128     655   609     154      86   476     280
## Barnesiella                 0     0       0     0       9       1     6       0
## Bifidobacterium             0     0       0     0       0       0     0       0
## Butyricicoccus              2     4       1     2       1       1     0       1
## Butyricimonas               0     0       2     2       6       0     2       0
## Catenibacterium             0     2       0     0       0       0     0       0
## Clostridium                 0     0       0     0       8       4     2       1
## Collinsella                 8     8       0     0      14       6     6       5
## Coprobacillus               9     1      19    29       1      18    11       1
## Coprococcus                27     8       9     9      54       7    20       1
## Desulfovibrio               0     0       0     0       0       0     3       0
## Dialister                  14     2       0     0       0       0     0       0
## Dorea                       1    10       8    15      15      43    42      10
## Eggerthella                 1     0       1     3       1       0     0       0
## Eubacterium                 4     0       2     1      21      88    23       0
## Faecalibacterium          301   120      64   167     215      37   391      72
## Gordonibacter               0     0       0     0       1       0     0       0
## Hespellia                  20    21       2     3       0       0     0       6
## Holdemania                  1     1       0     0       0       0     0       0
## Lachnobacterium             0     0       0     0       0       0     2       2
## Lactobacillus               0     0       0     0       3      24    10       0
## Lactococcus                 0     0       0     0       0       6     0       0
## Lactonifactor               0     0       3    15       0       0     4       1
## Marvinbryantia              0     0       1     0       0       1     3       0
## Megasphaera                 0     0       0     5       0       0     0       4
## Odoribacter                 0     0       8     3       9       6     9       1
## Oribacterium                6     3       0     0       0       0     7       0
## Oscillibacter              45    14      47    33      59      27    35       1
## Parabacteroides             0     0       2     4       2       0     0       0
## Paraprevotella              0     0       0     0       6       2    20       0
## Parasutterella              0     0      11    20       0       0     1       0
## Phascolarctobacterium       0     0      38    43       0       0     0       6
## Prevotella                  0     0       0     0     202       2   147       0
## Pseudobutyrivibrio          0     0       0     0       4       0     0       2
## Roseburia                  90    28      40   159     119      43    27      13
## Ruminococcus               64     0      23    21     191      26    95       0
## Slackia                     0     0       0     0       0       0     0       0
## Sporacetigenium             0     5       0     0       6       1     2       0
## Sporobacter                 0     0       0     0       0       0     0       0
## Streptococcus               5     0       3     1       1       1     0       1
## Subdoligranulum            75    39      25    64      71     284    22       5
## Sutterella                  4     0       0     0       0       0    69       0
## Tepidibacter                0     0       0     0       0       0     0       0
## Turicibacter                1     2       0     0       0       0     5       0
## Uknown                    276   237     230   291     696     461   234     215
## Veillonella                 0     0       0     0       0       0     0       0
##                       TS138 TS139.2 TS139 TS13 TS140.2 TS140 TS141.2 TS141
## Acetanaerobacterium       0       0     0    6       1     0       0     0
## Acetivibrio               1       0     0    0       0     0       0     0
## Acidaminococcus          16       1     0    0       0     0       8     2
## Akkermansia               0       0     0    0       0     0       0     0
## Alistipes                 2      16    10   78       1     5      25    13
## Allisonella               0       0     0    0       1     0       2     3
## Anaerostipes              0       0     0    1       0     0       0     0
## Anaerotruncus            18       3    18    3       5     0       0     1
## Asaccharobacter           0       0     0    0       0     0       0     0
## Bacteroides             137     338    20 1065      54    97     347   327
## Barnesiella               0       4     3   58       5     0       3     0
## Bifidobacterium           0       1     1    0       0     0       0     0
## Butyricicoccus            5       1     1    0       2     3       2     4
## Butyricimonas             0       0     1    2       1     1       8     1
## Catenibacterium           0     126   104    0       5     5      17     5
## Clostridium               3       6     0    3       1     0       1     0
## Collinsella              15       3    34    2       1     0       2     0
## Coprobacillus            13       0     1   70       0     5       0     4
## Coprococcus              16      47    26  171       8    47      42    22
## Desulfovibrio             0       7     1   14       0     5       1     1
## Dialister                 0       0     0    0       0     1       0     0
## Dorea                    94      10   108   43       2    31      16     1
## Eggerthella               0       0     0    2       0     0       0     0
## Eubacterium               0     104    85    0       6    50      42    35
## Faecalibacterium        119     146   251  648     325   686     322   488
## Gordonibacter             0       0     0    0       0     0       0     0
## Hespellia                 8       1    23   11       0     0       9     0
## Holdemania                0       1     0    3       0     0       0     0
## Lachnobacterium           0       0     0   73       2    11       0     0
## Lactobacillus             0       0     0    0       5     3       1     0
## Lactococcus               0       0     1    0       0     0       0     0
## Lactonifactor             8       4     0   18       3     0       1     8
## Marvinbryantia            0       0     1    0       1     0       0     0
## Megasphaera              25       0     0    0       0     0       0     0
## Odoribacter               1       0     2    7       1     5       2     2
## Oribacterium              1       0     0    0       0     1       5     4
## Oscillibacter             8      16     5  122      37    49      35     6
## Parabacteroides           0       1     0    2       0     0       1     0
## Paraprevotella            0       4     1    0       5     5       9     7
## Parasutterella            0       0     0    0       0     0       0     0
## Phascolarctobacterium    21       0     0   35       0     0       0     0
## Prevotella                0       0     0    0     243    49       0     0
## Pseudobutyrivibrio        2       7     2    0       0     0       1     1
## Roseburia                40     141    41   55      13     6     123    58
## Ruminococcus             56      33   119   88       5    29      14     9
## Slackia                   0       0     0    1       0     0       0     0
## Sporacetigenium           0       0     0    0       5     0       0     4
## Sporobacter               0       2     2    0       1     3       0     0
## Streptococcus             1       2     9    7       1     0      16    16
## Subdoligranulum          42      25   163  327      15    14      21     8
## Sutterella                0      29     4   11      11    19      71    77
## Tepidibacter              0       0    11    0       0     0       0     0
## Turicibacter              0       1     2    1       0     0       0     0
## Uknown                  952     387   690 1308     160   337     278   177
## Veillonella               0       0     0    0       0     0       0     4
##                       TS142.2 TS142 TS143.2 TS143 TS144.2 TS144 TS145.2 TS145
## Acetanaerobacterium         1     1       0     0       1     1       6     1
## Acetivibrio                 3     0       0     1       0     0       2     0
## Acidaminococcus             0     0      21    30      24     2       0     0
## Akkermansia                 0     0       0     0       0     0      24    15
## Alistipes                 150    15      44    18      70    16     312    53
## Allisonella                 0     0       1     1       2     0       0     0
## Anaerostipes                0     0       0     1       0     0       1     0
## Anaerotruncus               2     1       2    52      76    21      14    28
## Asaccharobacter             5     3       0     0       0     0       0     0
## Bacteroides               936   137     525   459     461   168     651   540
## Barnesiella                 0     0       2     6       0     0       1     0
## Bifidobacterium             0     1       0     1       0     1       1     0
## Butyricicoccus              0     1       2     7       0     2       0     0
## Butyricimonas              22     7       3     2       9     3       0     0
## Catenibacterium            22    31       6     0      21    14       0     0
## Clostridium                 1     1       0     2       1     0       2     0
## Collinsella                10    55       6    41      13     5       0     0
## Coprobacillus               0     0       0     0       1     1      23    53
## Coprococcus                 6     6      14    53      18    10      29    29
## Desulfovibrio               2     0       0     0       1     1       0     0
## Dialister                   4     6      16    23       0     0       0     0
## Dorea                      26    90      12    67      30    16       7     0
## Eggerthella                 0     0       1     0       0     0       4     1
## Eubacterium                51   105       1     1       9    13       0     0
## Faecalibacterium          165   191     348   229     124   185      67   116
## Gordonibacter               0     0       0     0       0     0       0     0
## Hespellia                  29     1      12     1       0     0       1     0
## Holdemania                  0     0       0     1       0     0       2     1
## Lachnobacterium             0     0       0     0       0     0       0     0
## Lactobacillus               7     6       3     2      27    17       0     0
## Lactococcus                 0     0       0     0       0     0       0     0
## Lactonifactor               0     0       0     0       3     3       0     0
## Marvinbryantia              0     0       0     0       0     1       0     1
## Megasphaera                 0     0      27    51      92    26       0     0
## Odoribacter                 0     0       3     0       6     2       6     3
## Oribacterium                0     0       0     0       0     1       0     0
## Oscillibacter              77     2      18    36      68    21     431   112
## Parabacteroides             1     0       3     2       1     0      14     3
## Paraprevotella              0     0       2     0       3     2       0     0
## Parasutterella             11     0       0     2       1     0       0     0
## Phascolarctobacterium       0     0       0     0       0     0      39    40
## Prevotella                  2     0      41     1     146   168       0     0
## Pseudobutyrivibrio          0     0       0     0       0     0       0     0
## Roseburia                   6    13      35   108      25    38       8    30
## Ruminococcus                0     0      15     1       5    30      23     0
## Slackia                     0     0       0     0       0     0       0     0
## Sporacetigenium             0     0       2    16       0     0       0     0
## Sporobacter                 9     0       1     0       0     0       0     4
## Streptococcus               3     2      12     2      20    69       0     2
## Subdoligranulum            53    78     115   134      87    33      15    43
## Sutterella                  0     0       0     0      11     3       4    33
## Tepidibacter                0     0       0     0       0     1       0     0
## Turicibacter                1     1       0     2       0     0       1     0
## Uknown                    198   412     364   485     269   153     304   240
## Veillonella                 0     0       1     0       2     0       0     0
##                       TS146.2 TS146 TS147 TS148 TS149 TS14 TS150 TS151.2 TS151
## Acetanaerobacterium         5     4     0     2     0    2     2       2     0
## Acetivibrio                 2     0     0     0     0    0     0       0     0
## Acidaminococcus             0     0     0     0     0    0     0       0     0
## Akkermansia                 8     0     0     4     0    1    28       8     0
## Alistipes                 205   264    70    51    33   28     3     156   112
## Allisonella                 0     0     1     0     0    0     0       0     0
## Anaerostipes                0     0     0     0     1    2     3       0     0
## Anaerotruncus             120   105     0    44    56    3     2       7     2
## Asaccharobacter             1     4     0     0     0    9     0       0     0
## Bacteroides               490   342   549   719  1775   97    29     450   387
## Barnesiella                44    47     0     9     9    1     0       0     0
## Bifidobacterium             0     0     0     0     0    7     0       0     1
## Butyricicoccus              1     1     2     4    12    2    22       0     2
## Butyricimonas               7     7     0     4     8    1     0       0     0
## Catenibacterium            29    41     0     0     0    0     0      30   122
## Clostridium                 1     6     0     1     0    2     1       0     4
## Collinsella                 0     0     8     1    20  106    71       0     0
## Coprobacillus               1     1     1     9    41  653     0       0     0
## Coprococcus                10    14    13    41    24   97    47      18    34
## Desulfovibrio               0     0     0     7     3    1     0       0     0
## Dialister                   0     1     0    32     0   77     1      16    17
## Dorea                      23    43    44    26    63  598   229      18    12
## Eggerthella                 0     0     2     0     0   70     0       2     4
## Eubacterium                15    23     0     1     1    0     0       1     0
## Faecalibacterium           92   204   145   392   989  843     0     158   235
## Gordonibacter               0     0     0     0     2    1     2       0     0
## Hespellia                   1     0     4     0     2    9    53      15     0
## Holdemania                  1     2     1     1     3    1     0       0     0
## Lachnobacterium             2     0     0     0     0    0     2       0     0
## Lactobacillus               0     0     0     9    32    0     0       0     0
## Lactococcus                 0     0     0     0     4    4     0       0     0
## Lactonifactor               3     0     0    20     9    0    14       5     6
## Marvinbryantia              1     0     0     0     0    1     0       0     0
## Megasphaera                 0     0     0    72     0    0     1       0     0
## Odoribacter                 0     0     0    10    10    0     0       7    10
## Oribacterium                2     2     0     7     0    3     0       0     0
## Oscillibacter             287   133    15    53    73   18     0     131    38
## Parabacteroides             0     0     4    46     3    0     1       3     0
## Paraprevotella              0     0     0    14   242    0     3       0     0
## Parasutterella              2     4    16     0     0    0     2       2    12
## Phascolarctobacterium      10    21    14     1     0    1    10       0     0
## Prevotella                  0     0     0     1     0    5     0       0     1
## Pseudobutyrivibrio          0     0     0     2     1    0     0       0     0
## Roseburia                   7     9   187    57   374  132   133       4    12
## Ruminococcus               86    76    55   113     0  639     1      77    38
## Slackia                     2     1     0     0     0    0     3       0     0
## Sporacetigenium             0     8     0    16    19   75     0       0    45
## Sporobacter                 0     0     1     0     0    0     0       0     1
## Streptococcus               0     0     9     1    29   72     9       5     1
## Subdoligranulum            19     8    67     8    37  934    68      28    23
## Sutterella                  0     0     0    73    74    3     0       0     0
## Tepidibacter                0     0     0     0     0    0     0       5     9
## Turicibacter                0     1     0     2    19   39    31       0     0
## Uknown                    344   320   485  1490  1067 3582  3099     588   734
## Veillonella                 0     0     0     0     0    0     0       0     0
##                       TS152.2 TS152 TS154.2 TS155.2 TS155 TS156.2 TS156 TS160.2
## Acetanaerobacterium         4     4       0       2     0       2     1       0
## Acetivibrio                 0     0       0       0     0       0     1       0
## Acidaminococcus             0     0       0       0     0       0     0       0
## Akkermansia                13     1       1       0     0       0     0       0
## Alistipes                  89    88       7       5     5      29     5       4
## Allisonella                 1     0       0       0     0       0     0       0
## Anaerostipes                0     5       0       0     0       0     0       0
## Anaerotruncus               7     8       5      77     5       3     2       0
## Asaccharobacter             0     0       0       0     0       2     0       0
## Bacteroides               264   479     219     277   279     171   352      29
## Barnesiella                 0     0       7       9     3       0     0       0
## Bifidobacterium             0     0       0       0     0       1     0       0
## Butyricicoccus              1     2       0       1     2       6     1       0
## Butyricimonas               3     0       1       1     2       0     2       0
## Catenibacterium             0     0      20       0     0       0     0       0
## Clostridium                 5     0       0      15    16       2     1       1
## Collinsella                 0     0      35       8    12      10     1       1
## Coprobacillus               0     0       0       0     0       0     0       0
## Coprococcus                36    16      20      29    14       2    15      23
## Desulfovibrio               0     0       2       0     0       1     0       1
## Dialister                   0     0       0       0     0      10    13       0
## Dorea                       2     1      11      10    25      20     5       8
## Eggerthella                 3     4       0       2     4       0     0       0
## Eubacterium                 9     7       2      23     6      19    31       2
## Faecalibacterium          105    97      23      20    43     284   152     176
## Gordonibacter               0     0       0       0     0       0     0       0
## Hespellia                  38     2       6       7     0       1     1       0
## Holdemania                  0     2       1       1     1       0     0       0
## Lachnobacterium             0     0       0       0     0       0    15       0
## Lactobacillus               1     0       0       0     0       9     0       0
## Lactococcus                 1     1       0       0     0       0     1       0
## Lactonifactor               0     0       0       0     0       0     0       0
## Marvinbryantia              0     0       0       0     0       0     0       0
## Megasphaera                 0     0       0       0     0       0     0       0
## Odoribacter                 1     5       2       1     1       0     0       4
## Oribacterium                0     0       0       0     0       0     0       0
## Oscillibacter             136   135      70     108   191      66    39      13
## Parabacteroides             3     4       5       3     2       0     1       2
## Paraprevotella              0     0       0       3     5       0     0       2
## Parasutterella              1     5       4       0     0       0     0       0
## Phascolarctobacterium      18    14       0       0     0       0     0       0
## Prevotella                  0     0      13      31    38      11     2     229
## Pseudobutyrivibrio          0     1       0       0     0       1     0       0
## Roseburia                   2     1      20      23    13      35    32       5
## Ruminococcus              102    59      51      62    40      59    26       7
## Slackia                     0     0       0       0     0       1     0       0
## Sporacetigenium             0     0       4       0     7      12     0       0
## Sporobacter                 0     0       6       3     0       5     0       0
## Streptococcus               2     0       0       0     0       4    19       0
## Subdoligranulum            68   112      24      14    18      71    37       9
## Sutterella                  0     0       0       0     3       0     8       0
## Tepidibacter                0     0       0      18     3       0     0       0
## Turicibacter                0     0       0       1     0       0     0       0
## Uknown                    373   559     342     672   599     370   265     136
## Veillonella                 0     0       0       0     0       0     0       0
##                       TS160 TS161 TS162.2 TS162 TS163.2 TS163 TS164.2 TS164
## Acetanaerobacterium       0     2       0     0       0     0       3     0
## Acetivibrio               0     0       0     0       0     0       0     0
## Acidaminococcus           0     0       0     0       0     0       0     0
## Akkermansia               0     0       0     0       0     0       0     0
## Alistipes                 6    34       2    19      74    49      22    11
## Allisonella               0     0       0     0       0     0       0     0
## Anaerostipes              0     0       0     0       0     0       0     0
## Anaerotruncus             2     0       0     3       2     2       3     4
## Asaccharobacter           0     0       0     0       0     1       0     0
## Bacteroides              42   256    1133  1427     513   307     584   628
## Barnesiella               2     0       0     0      14     9       0     0
## Bifidobacterium           0     0       0     0       0     0       1     1
## Butyricicoccus            6     0       1     1       1     1       1     2
## Butyricimonas             2     0       0     0       0     0       0     0
## Catenibacterium           0     0       0     0       0     0       0     0
## Clostridium               1     0       0     0       5     0       9     3
## Collinsella               2     1       0     0       0     0       0     0
## Coprobacillus             0     0       0     1       0     2       0     0
## Coprococcus              30     1       5     2      29     0      46    36
## Desulfovibrio             0     0       0     0       0     0       2     2
## Dialister                 0    20      35     9       8     5       0     0
## Dorea                     8     9      77    89      26    23      17     7
## Eggerthella               0     0       1     6       1     1       5     1
## Eubacterium               3     0       0     9       2    14     149   114
## Faecalibacterium        218   405     325   110     534   187     335   190
## Gordonibacter             0     0       0     0       0     0       0     0
## Hespellia                 0     0       4    90       3     2      75     1
## Holdemania                0     0       0     0       0     1       1     1
## Lachnobacterium           0     0       0     0       0     0       0    10
## Lactobacillus             0     0       0     0       0     0       0     0
## Lactococcus               0     0       0     0       0     0       2     0
## Lactonifactor             0     0       0     0       0     0       0     0
## Marvinbryantia            0     0       0     0       0     0       0     0
## Megasphaera               0     0       0     0       0     0       0     0
## Odoribacter               4     0       0     0       3     5       0     1
## Oribacterium              0     0       0     0       5     0       0     0
## Oscillibacter            23    10       3    19       4     5     154    68
## Parabacteroides           1     8       3     4       0     0       1     0
## Paraprevotella            5     0       0     0       0     0       0     0
## Parasutterella            0     0       0     0       0     0       0     0
## Phascolarctobacterium     0     0       0     0       0     0      21     5
## Prevotella              189     0       0     1       0     0       1     1
## Pseudobutyrivibrio        0     0       0     0       0     1       1     0
## Roseburia                 8    59      55    38      52    75      54    51
## Ruminococcus             10     7       0     0       0     0       0     0
## Slackia                   0     0       0     0       0     0       0     0
## Sporacetigenium           0     0       0     0       6     0       3     0
## Sporobacter               3     0       0     0       0     0       0     0
## Streptococcus             0     0       0     0       3     4       0     3
## Subdoligranulum          12     7       0     0      26    21       8    19
## Sutterella                0     4       0     0       0     0       0     0
## Tepidibacter              0     0       2     0       0     0       9     0
## Turicibacter              0     0       0     0       2     1      15    10
## Uknown                  183   209     479   235     309   314     430   425
## Veillonella               2     5       2     5       0     0       0     0
##                       TS165.2 TS165 TS166.2 TS166 TS167.2 TS167 TS168.2 TS168
## Acetanaerobacterium         0     6       0     1       0     0       0     0
## Acetivibrio                 0     0       1     0       2     2       0     4
## Acidaminococcus             0     0       0     0       0     0       4     1
## Akkermansia                 1     9       0     0       1     3       1     4
## Alistipes                 378   330       5    41      27    50     105    43
## Allisonella                 0     0       0     0       0     0       0     0
## Anaerostipes                0     0       0     0       0     0       0     0
## Anaerotruncus              35    13       3     1       0     0       8     3
## Asaccharobacter             0     0       1     1       0     0       0     0
## Bacteroides               403   174      22   723      59    87     435   198
## Barnesiella                15    23       0     0       0     0       0     0
## Bifidobacterium             0     1       0     0       2     0       0     0
## Butyricicoccus              0     0       2     3       0     2       4     4
## Butyricimonas              16     3       0     3       0     0      11     5
## Catenibacterium            61    72       0     0       0     0       0     0
## Clostridium                 7   126       1     0      82    26       0     4
## Collinsella                 4    13      25     4     139    14      13     8
## Coprobacillus               0     4      33     2      27    10      26    16
## Coprococcus                 6    21       4    19      33    11      90   124
## Desulfovibrio               0     0       0     0       0     0       0     0
## Dialister                   2    10       0     0       3    12       0     0
## Dorea                       6    27      48    16     110    31       5     7
## Eggerthella                 3     5      33     0       2     0       1     0
## Eubacterium                 0     4     113    27       0     0       0     0
## Faecalibacterium           59     7      17   159      22   128     115    94
## Gordonibacter               0     0       0     0       4     1       2     3
## Hespellia                   1     4       2    15       2     8      22     2
## Holdemania                  1     1       0     0       0     0       1     0
## Lachnobacterium             0     6       0     0       1     2       0    19
## Lactobacillus               3     2       2     0       0     1       0     0
## Lactococcus                 2     0       0     1       0     0       0     0
## Lactonifactor               0     0       0     3       0     0       8     3
## Marvinbryantia              0     0       0     0       0     0       0     0
## Megasphaera                 0     0       0     0       0    25       0     0
## Odoribacter                20     3       1     2       3     1       0     0
## Oribacterium                0     0       0     0       0     0       0     3
## Oscillibacter             126    71       2    32       9    17     116    54
## Parabacteroides             1     1       0     1       0     0       0     0
## Paraprevotella              6     5       2     6       0     0       0     0
## Parasutterella              0     0       0    19       1     2       7     3
## Phascolarctobacterium       6     5       0    36       0     0       0     0
## Prevotella                  0     0       0     1       0     1       0     5
## Pseudobutyrivibrio          0     0       0     0       0     0       0     0
## Roseburia                   5     2      50   114      26    20      24    72
## Ruminococcus               59    58      60     0      10    83      46    12
## Slackia                     0     0       0     0       0     0       0     0
## Sporacetigenium             5     0       3     0      55     4       3     0
## Sporobacter                 0     0       0     4      20     0       0     0
## Streptococcus               0     1      31     1       3     0       0     2
## Subdoligranulum            18   105      19    16     297   314     241   172
## Sutterella                  0     0       0     0       0     0       0     0
## Tepidibacter                0     0       0     4      19     0       0     0
## Turicibacter                1     7       6     2      23     1       0     0
## Uknown                    772   760    1007   347    1046   430     815   444
## Veillonella                 0     0       0     0       0     0       0     0
##                       TS169.2 TS169 TS16 TS170.2 TS178.2 TS178 TS179.2 TS179
## Acetanaerobacterium         0     1    3       0       0     0       0     0
## Acetivibrio                 1     0    7       0       0     0       0     0
## Acidaminococcus             0     0    0       0       0     2       1     0
## Akkermansia                 4     1    0       0       0     0       0     0
## Alistipes                  23    24   61       2       0     0       4     2
## Allisonella                 0     0    0       0       1     0       0     0
## Anaerostipes                0     0    0       3       0     0       0     0
## Anaerotruncus               3    24   10       8       4     3       0     2
## Asaccharobacter             0     0    0       0       0     0       0     0
## Bacteroides               193   244 1206     616      16    28      31     8
## Barnesiella                19     3    0       0       0     0       0     0
## Bifidobacterium             2     0    0       1       0     4       0     3
## Butyricicoccus              4     0    0       0       3     3       1     1
## Butyricimonas               0     0    4       0       0     0       0     1
## Catenibacterium           135   172    0       0      60   232      28    46
## Clostridium                 1     1    0      12       5     9       1     3
## Collinsella                33    10    0       9      21   106      12    44
## Coprobacillus               0     0    0      10       0     2       0     0
## Coprococcus                88    67   15       2      11     1      12     5
## Desulfovibrio               0     0    4       1       0     0       0     0
## Dialister                   0     0    0       6       0     0      32    13
## Dorea                      23    19   23      15      46   122      26    71
## Eggerthella                 0     0    4       1       0     2       1     0
## Eubacterium                87   107    0     124      69   110      35   124
## Faecalibacterium          169   292  443      42     128    72     228   131
## Gordonibacter               0     0    0       0       0     0       0     0
## Hespellia                   0     5   10       8       0     0       0     0
## Holdemania                  0     0    1       0       0     0       0     0
## Lachnobacterium             2     0   46       0       0     0       0     0
## Lactobacillus               0     0    0       1      13   302      17     2
## Lactococcus                 0     1    0       0       0     0       1     0
## Lactonifactor               0     0    6       1       1     1       0     0
## Marvinbryantia              0     1    0       0       0     0       0     0
## Megasphaera                 0     0    0      17      37    89      20     6
## Odoribacter                 2     3    0       8       1     0       0     1
## Oribacterium                0     1    0       1       0     0       1     0
## Oscillibacter              70    85  241       8       5     0      51    30
## Parabacteroides             1     0   44       1       0     1       0     0
## Paraprevotella              0     0    0       0       1     0       0     0
## Parasutterella              0     0   11       5       0     0       0     0
## Phascolarctobacterium       0     0    0       0       0     0       0     0
## Prevotella                 19    46   12       2     139     1     203     2
## Pseudobutyrivibrio          0     0    1       0       0     0       0     0
## Roseburia                  20    10   51      42      52     8      19     2
## Ruminococcus               60   121  134       4      24     1       6     3
## Slackia                     0     0    0       0       2     1       0     0
## Sporacetigenium             0     0    0       0       4    26       0    25
## Sporobacter                 0     0    0       0       2     0       2     0
## Streptococcus               2     2   34       9      10     9       3     5
## Subdoligranulum            45    66  154      24      20    32      27    82
## Sutterella                  2     0    0       0       1     0       8     0
## Tepidibacter                0    14    0      14       9     0       0     0
## Turicibacter                0     2    0       1       0     1       0     0
## Uknown                    545   511  527     302     372   342     338   718
## Veillonella                 0     2    0       1       0     1       1     0
##                       TS17 TS180 TS181.2 TS181 TS182.2 TS183.2 TS183 TS184
## Acetanaerobacterium      0     1       0     1       1       0     0     0
## Acetivibrio              0     0       0     0       0       0     1     0
## Acidaminococcus          0     0       0     0      13       0     0     0
## Akkermansia              0     0       0     3       0       0     0     0
## Alistipes                3     9       6    10      14       1    45     4
## Allisonella              0     0       0     0       0       1     1     0
## Anaerostipes             0     0       0     0       0       0     0     0
## Anaerotruncus            3     2      18    44       0       0    47     2
## Asaccharobacter          0     0       0     0       0       0     0     0
## Bacteroides            237   826     656   703     136     692  1087   573
## Barnesiella              0     0       0     0       3       0    27     0
## Bifidobacterium          0     0       1     0       0       0     0     0
## Butyricicoccus           0     2       1     2       2       6     0     1
## Butyricimonas            1     0       1     7       2       0     9     0
## Catenibacterium          0     0      99    87     139     316    84     0
## Clostridium              1     0       0    17      28       0     0     0
## Collinsella              0     1       3     8      18      18     3    35
## Coprobacillus            0     1       0     0       1       0     0    15
## Coprococcus              0     8       1    29      28       9     3    33
## Desulfovibrio            0     0       0     0       0       0     0     1
## Dialister                0     0      13    17       0       0     0     5
## Dorea                    2     6       2     4      27      33    15    15
## Eggerthella              0     0       0     1       0       0     0     0
## Eubacterium              1     0       4     0       1       0     2     0
## Faecalibacterium        44   196     230    85     527       8   103    30
## Gordonibacter            0     0       0     0       0       0     0     0
## Hespellia                0     0       5     0       0       2     8     7
## Holdemania               0     1       0     1       0       0     1     0
## Lachnobacterium          0     3       0     0       0       0     0     0
## Lactobacillus            0     0       0     0      19       0     0     0
## Lactococcus              0     0       1     0       0       0     0     0
## Lactonifactor            4     0       0     0       0      14     0     0
## Marvinbryantia           2     0       0     0       0       0     0     0
## Megasphaera              0     0      12    23      30      81    21     0
## Odoribacter              2     7       1     5       2       0     3     0
## Oribacterium             0     0       0     1       0       0     0     0
## Oscillibacter           21    17      18    17      24       2    13     2
## Parabacteroides          1     0       0     0       2       1     1     0
## Paraprevotella           0    16       0     0       0       0     0     0
## Parasutterella           0     0       9     7       0       1     0    11
## Phascolarctobacterium    2    28       0     0       0      12    36    11
## Prevotella               0     0       1     0     634       0     0     0
## Pseudobutyrivibrio       0     0       0     1       0       1     0     0
## Roseburia               17    31      13    38     101      29    33   223
## Ruminococcus            45     8      24    67      55       3    10     8
## Slackia                  0     0       0     0       0       1     0     0
## Sporacetigenium          0     0       0     1       0       3     0     5
## Sporobacter              0     0       0     0       1       0     0     0
## Streptococcus            0     0       9     5       5      38     6     2
## Subdoligranulum         14    54      26    65      96       2    27    10
## Sutterella               0     0       1     0      16       0     0     0
## Tepidibacter             0     0       0     6       0       1     0     0
## Turicibacter             1     0       1     0       9       0     0     1
## Uknown                 151   265     190   395     723     298   174   340
## Veillonella              0     1       1     1       0       1     0     1
##                       TS185.2 TS185 TS19.2 TS190 TS191 TS193.2 TS193 TS194.2
## Acetanaerobacterium         0     0      0     0     0       2     1       0
## Acetivibrio                 0     0      0     0     0       0     2       0
## Acidaminococcus             0     0      5     0     0       0     0       0
## Akkermansia                 0     0      0     0     0       0     0       0
## Alistipes                  23    17      2    25     4      19     3      18
## Allisonella                 1     1      0     0     0       0     0       0
## Anaerostipes                0     0      1     0     0       0     1       0
## Anaerotruncus               4     2     19     0     1       0     1       0
## Asaccharobacter             2     0      0     0     0       0     0       0
## Bacteroides               515   355      0   327   745     137   203     230
## Barnesiella                12    17      0     7     0       2     0       0
## Bifidobacterium             0     0      2     1     0       0     0       0
## Butyricicoccus              2     3      0     7     1      13     0       0
## Butyricimonas               0     0      0     2     0       0     0       0
## Catenibacterium             0     0      0     0     0       0     0       0
## Clostridium                 0     0      1     1     2       0     0       0
## Collinsella                 0     0    104     8     2       0     2       1
## Coprobacillus               0    36    263    13     0      46    50       0
## Coprococcus                46    15     36    90    30      33    13      10
## Desulfovibrio               0     0      0     0     0       0     1       0
## Dialister                   0     7      0     0    10       0     0       0
## Dorea                       5    22    261     0    18      16     4      11
## Eggerthella                 0     2      0     0     0       1     0       0
## Eubacterium                 0     1      0     0    44       0     1       0
## Faecalibacterium          389   488      2   272   129     202    97     264
## Gordonibacter               0     0      0     0     0       0     0       0
## Hespellia                   0     1      4     7    13       4     4       0
## Holdemania                  0     0      0     0     0       0     0       1
## Lachnobacterium             1     0      0     0     0       0     0       0
## Lactobacillus               0     0      0     0     0       0     0       0
## Lactococcus                 0     0      0     0     0       0     0       0
## Lactonifactor               3     3      1     5     3       6     3       0
## Marvinbryantia              0     0      0     0     0       0     0       0
## Megasphaera                 0     0    124     0     1       0     0       0
## Odoribacter                 0     0      0     4     0       5     1       2
## Oribacterium                0     0      0     0     0       0     0       5
## Oscillibacter              25    10      0    19     3       5     5      13
## Parabacteroides             8     0      0     1     1       0     4       1
## Paraprevotella              0     0      0     0     0       0     1       0
## Parasutterella              7     6      0     0     2       0     0       2
## Phascolarctobacterium      18     2      4     0     0      10     1       9
## Prevotella                  0     0      0     0     0      13    40       0
## Pseudobutyrivibrio          0     0      2     0     0       0     0      34
## Roseburia                 121    66     84    21    33      45    53       5
## Ruminococcus                0     1      0    47     0      38    17      30
## Slackia                     0     0      0     0     0       1     1       0
## Sporacetigenium             0     1      0     0     0       0     0       0
## Sporobacter                 0     0      0     0     0       0     0       0
## Streptococcus               4     4     13     5     1       1     0       0
## Subdoligranulum            52   115      2    75    43     221    51      48
## Sutterella                  0     0      0    11     0       0     0       0
## Tepidibacter                0     1      0     0     7       0     0       0
## Turicibacter                0     0      7     0     0       0     0       0
## Uknown                    466   573   1258   489   282     350   310     371
## Veillonella                 0     0      1     0     0       0     0       0
##                       TS194 TS195.2 TS195 TS19 TS1 TS2.2 TS20.2 TS20 TS21.2
## Acetanaerobacterium       4       0     0    3   1     3      0    0      0
## Acetivibrio               0       0     0    0   5     0      0    0      0
## Acidaminococcus           0       0     0   57   0     0      0    0     11
## Akkermansia               0       9     2    0   5     4      4   18      7
## Alistipes                28      68    71   21 120    99     58  284      0
## Allisonella               0       0     0    0   0     0      0    0      3
## Anaerostipes              0       0     0    0   0     0      0    3      0
## Anaerotruncus             0      15    14   24 103     8     13   19      7
## Asaccharobacter           0       0     0    0   1     1      0    0     19
## Bacteroides             494     456   505 1075 886   268    217  620    301
## Barnesiella               0      15    10   17  71     0     18   53      0
## Bifidobacterium           0       0     0    0   0     0      1    0      0
## Butyricicoccus            0       8     2    0  13     0      1   12      3
## Butyricimonas             0       0     0    0   7     0      0    0      0
## Catenibacterium           0       0     0    0   0     0      0    0      0
## Clostridium               0       0     0   22   0     6      0    1      6
## Collinsella               0       0     0   27   0     9     12   39     91
## Coprobacillus             9      25    73  247  53     3     55  239   1102
## Coprococcus               5       0     2   87  15     8    110  218      2
## Desulfovibrio             0       0     1    0   4     0      0    3      0
## Dialister                 0       0     0    0  20     2     89  273     14
## Dorea                     2      20    34  143   5     1     46  159     51
## Eggerthella               0       0     0    2   0     0      2    1      0
## Eubacterium               0       4     1   10   2     0      0   12      1
## Faecalibacterium        149      83   171  133 460   120    404 1546     54
## Gordonibacter             0       0     0    0   5     0      0    0      0
## Hespellia                 0       0     1  274   7     1      2    5     31
## Holdemania                1       1     0    2   0     1      0    5      0
## Lachnobacterium           0       0     1   86  32     0      4    0      0
## Lactobacillus             0       0     0    1   0     0      0    0      0
## Lactococcus               0       0     0    0   0     0      0    1      0
## Lactonifactor             0       0     0    0   0     1     11   30     17
## Marvinbryantia            0       1     1    0   0     1      1    0      0
## Megasphaera               0       0     0  158   0     0     42    1     19
## Odoribacter               3       8    10    0   6     4      1   11      0
## Oribacterium              0       0     0    0   6     0      0    4      1
## Oscillibacter            23      35    28    3 242    67     21   15     10
## Parabacteroides           4       1    14    8   7     7      9   18      1
## Paraprevotella            0       0     0    0   0     0      0    0      0
## Parasutterella            7       0     0    0   1     1      0    0      0
## Phascolarctobacterium    19      76    15   63   0     0      0    0      0
## Prevotella                0       0     0    0   0     0    842 1632      0
## Pseudobutyrivibrio       18       0     0    1   0     0     13    1      0
## Roseburia                 1      75    11  602  28    15    120  482    154
## Ruminococcus             45      29    42    1 148    38     42  680     42
## Slackia                   1       0     0    0   1     0      0    0      2
## Sporacetigenium           2       0     0    0   2     8      0    0      7
## Sporobacter               0       0     0    0   0     0      0    0      0
## Streptococcus             3       1     1    5   0     1      4    2      6
## Subdoligranulum          24      89    58    2 242    10      2   20     60
## Sutterella                1      14     1    0   0     1      0    1      2
## Tepidibacter              0       0     0    3   0     0      2    0      0
## Turicibacter              0       0     0    2   1     2      0   14      0
## Uknown                  382     772   456 1446 783   434    373 2008    944
## Veillonella               0       0     0    1   0     0      0   61      0
##                       TS21 TS22 TS23 TS25.2 TS25 TS26.2 TS26 TS27.2 TS27 TS28
## Acetanaerobacterium      1    0    0      6    2      5    2      0    7    1
## Acetivibrio              0    2    0      0    1      0    0      4   37    2
## Acidaminococcus         19    0    0      0    0      0    0      0    0    0
## Akkermansia             12    0    0      8    0      4    6     22    1    0
## Alistipes               43   18    0    151  147    212  277     30   49   33
## Allisonella              0    0    0      0    0      0    0      0    0    0
## Anaerostipes             3    1    1      0    0      0    0      0    0    0
## Anaerotruncus           10    0    1      2   24      4   14      4   22   10
## Asaccharobacter          8    0    0      2    0     35    7      1    3    1
## Bacteroides           1681   75   95    455 1708    547 4973     42 1103  427
## Barnesiella              0    0    0      0    0      0    0     12   26    0
## Bifidobacterium          0    0    0      0    0      3    0      1    0    0
## Butyricicoccus          14    1    0      1    4      1    5      0    2   13
## Butyricimonas            0    0    0      0    0      0    0      0    0    0
## Catenibacterium          0    0    0      0    0      0    0      0    0    0
## Clostridium             15    0  171     24   19      7   54      0    3    0
## Collinsella             30  258    0     23    3     89    9     41   68    0
## Coprobacillus          862    2   24     23  124    113   62      0   10    3
## Coprococcus             25   13   38     79  148     54  170     43   29    2
## Desulfovibrio            8    7    0      0    1      0    8      0    2    0
## Dialister               60    0    0      0    0      1    0      0    0    4
## Dorea                   97   76    7     10   49     63   65      6   25   46
## Eggerthella              1   21   26      7    0     20    5      2    0    1
## Eubacterium              9   10    7      0    0      0    0     45  115    6
## Faecalibacterium       522  109   27    139  133    179  312     35   35  930
## Gordonibacter            3    2    0      4    1      1    1      1    0    0
## Hespellia               14  122    8     85   60     17   87      1  176   35
## Holdemania               7    0    0      1    3      0    0      1    1    0
## Lachnobacterium          0    0    0      0    0      2   12      0    1    0
## Lactobacillus            0    0    0      0    0      0    0      0   14    0
## Lactococcus              1   11   12      0    3      0    0      0    0    0
## Lactonifactor            2    9    0      0    0      0    0      1    1    3
## Marvinbryantia           0    0    2      0    1      0    0      1    0    0
## Megasphaera            117    0    0      0    0      0    0      0    0    0
## Odoribacter              0    5    0      4   10      4   18      2   21    0
## Oribacterium             0    0    0      0    0      0    0      0    0    0
## Oscillibacter           15    2    1    298  357     58  255    113  434   70
## Parabacteroides         10   10    2      9   44      2    7      1   35    9
## Paraprevotella           0    0    0      0    0      0    0      0    0    0
## Parasutterella           0    7    0      0    0      0    0      0    0    4
## Phascolarctobacterium    0   18    7     14   96     17  100      0    0   17
## Prevotella               2    0    0      0    1      1    2      0    1    0
## Pseudobutyrivibrio     349    0    0      0    0      0    2      4    0    0
## Roseburia              245  146   49     19   11     17   12     11   79   29
## Ruminococcus           125    0    0     85  265     39  188     23 1207  307
## Slackia                  0    0    0      0    0      0    0      9   10    0
## Sporacetigenium          6    0    3      5    0      5    6      0    0   11
## Sporobacter              0    0    0      0    1      0    0      1    0    0
## Streptococcus           12   97   47      8   61      8   19     21   55    4
## Subdoligranulum        377    0    3     42   60     36   67      9   87  261
## Sutterella              43    0    0      0    0      0    0      5    1    0
## Tepidibacter             0    2    0      0    0     18    7      0    0    0
## Turicibacter             0   36    8    149  157     84   35     86  475    1
## Uknown                2447 1748  548    884 1770   1137 2505    553 2093  722
## Veillonella              0    0    0      0    5      0    0      0    0    0
##                       TS29  TS2 TS30.2 TS31.2 TS31 TS32.2 TS32 TS34 TS35 TS37.2
## Acetanaerobacterium      1    4      1     12    8      0    1    0    3      0
## Acetivibrio              2    4      0      1    4      0    0    0    0      3
## Acidaminococcus          0    0      0      0    0      2   48    0    0     19
## Akkermansia              0   80      0      4   20      1    0    0    0      9
## Alistipes                6  303    157    133  512     36   93    4   11    114
## Allisonella              0    0      0      0    0      1    4    0    0      0
## Anaerostipes             1    2      1      0    0      0    0    0    0      0
## Anaerotruncus            7   53      2      3  212      1    4    2   10     11
## Asaccharobacter          0    5      0      0    2      0    1    0    0      0
## Bacteroides             44 1316    288   1562 2649    241  759   16   21    581
## Barnesiella              0    1     21      0    0     22   46    0    0      5
## Bifidobacterium          0    0      0      0    0      2    1    0    0      0
## Butyricicoccus           1   15      0      1    4      5   13    0    3      0
## Butyricimonas            0    0      0      8   33     14   17    1    4      5
## Catenibacterium          0    0      0      0    0    492  457   91   15      0
## Clostridium              2    0      4      4    1      6    0    0    7      1
## Collinsella              4  188     24      2   47     73   20   19   49    151
## Coprobacillus            0    0      1     25  271     14    0   26    6     82
## Coprococcus              1   26      4     11   83     10    7   29  163     33
## Desulfovibrio            0    0      2      0   25     33   22    2    6      0
## Dialister                0   16      0      0    0      0    0    0    0      0
## Dorea                   11   11      0     24  106     80   65   48   42     59
## Eggerthella              0    2      4      1    3      0    0    2    0      3
## Eubacterium              1    2      0     28  123     58   18    1  114     40
## Faecalibacterium       113  394    231    239 1626    322  840  119  361    387
## Gordonibacter            0    6      0      0    0      1    0    0    0      0
## Hespellia               26  340     27     57    1      4    1    0    2      0
## Holdemania               1    2      3      3    5      0    3    0    0      1
## Lachnobacterium          0    0      0      0    0     12    0    6    0      0
## Lactobacillus            0    0      0      0    0      1    0    0    0      0
## Lactococcus              1    0      0      0    2      0    0    0    0      0
## Lactonifactor            2    1      0      0    0     10   35    5   12      9
## Marvinbryantia           0    0      0      0    0      0    0    0    0      0
## Megasphaera              0    0      0      0    0      0    0    0    1      0
## Odoribacter              0   20      0      7   16      7    5    0    0      8
## Oribacterium             0    0      0      0    0      0    0    0    0      0
## Oscillibacter            2  819     29    133  465     30   21    3   24      8
## Parabacteroides          0   32     11      6   10      2    8    1    2      9
## Paraprevotella           0    0      3      0    0     15   42    0    3      0
## Parasutterella           3    2      0     11   27      0    0    0    6     26
## Phascolarctobacterium    1    1     11      0    0      0    0    0    0     50
## Prevotella               0    0      0      0    2    443  994   10  417      0
## Pseudobutyrivibrio       0    0      0      0    1      0    1    0    0      2
## Roseburia               30  116     33     24   83     55   37   36   82    212
## Ruminococcus             3  868     12     59  495     43    1   16   21     88
## Slackia                  0    0      0      0    0      5    1    0    2      0
## Sporacetigenium          0   28      0      3    0      0    1    0   13     29
## Sporobacter              0    0      0      0    0      0    0    0    0      0
## Streptococcus            0    1      2      1   30      8    1    2    9     13
## Subdoligranulum        409  194      0     90  381     31  107   36    8     21
## Sutterella               0    0      0      0    1     79  154    0    0      1
## Tepidibacter             0    0      0      0    0      0    3    0    0      9
## Turicibacter             0   16      0      0    1      0    0    0   12      5
## Uknown                 413 2757    717   1058 3333    735  961  313  649    728
## Veillonella              0    0      0      0    0      0    0    0    0      0
##                       TS37 TS39.2 TS39 TS4.2 TS43 TS44 TS49.2 TS49  TS4 TS5.2
## Acetanaerobacterium      0      0    2     2    2    0      1    0    0     0
## Acetivibrio              0      0    1     0    0    0      0    0    0     0
## Acidaminococcus          9      0    0     0   23    0      0    0    0     0
## Akkermansia              1      6    8     0    6    0     10    0    0     1
## Alistipes               22     56  104    30   20    7     15    0   61     0
## Allisonella              0      0    0     0    2    0      0    0    0     0
## Anaerostipes             0      0    0     0    0    0      1    0    0     0
## Anaerotruncus            5      0    4     2    3    5      0    2    0     0
## Asaccharobacter          0      1    0     5    1    0      0    0    0     0
## Bacteroides            993    259 1529   772  424   42    281  256 2379    47
## Barnesiella              1     42   52    17   13    2      0    0   37     0
## Bifidobacterium          1      0    2     0    2    0      0    3    0     0
## Butyricicoccus           2      0    6     0    5    1      2    0    5     2
## Butyricimonas            2     17    8     1   13    0      0    0    1     0
## Catenibacterium          0      0    3     0    0    0    660  768    0     0
## Clostridium              6      5    6     8    0    0      2    3    2    21
## Collinsella             16     50    5     0   61  155    135   54    0    54
## Coprobacillus           20     60   29     7   87    0      0   11    3   551
## Coprococcus              6     16   24    46   16    7      6  486   71     9
## Desulfovibrio            1      4    5     1    0    0      0    0    4     0
## Dialister                0      0    0    17    0   23      0    1   58     0
## Dorea                   39     29   28    25   55   18    134  111    3   108
## Eggerthella              0      0    0     1    2    3      2    3    2     7
## Eubacterium             49      3    0     9    0   27    234  543    4     1
## Faecalibacterium        94    147  454    90   48   46      3   51  176    62
## Gordonibacter            0      0    0     0    2    0      0    0    0     0
## Hespellia                3      1    6    17    8    7     37   52  115    11
## Holdemania               1      0    2     1    0    0      0    6    2     0
## Lachnobacterium          0      0    1     0    0    0      0    0    7     0
## Lactobacillus            0      1    0     0    0    0      0    0    0     0
## Lactococcus              0      0    0     0    0    0      0    1    0     2
## Lactonifactor            2     21    1     9    6    0      2    7   26     2
## Marvinbryantia           0      0    1     0    0    0      0    0    0     0
## Megasphaera              0      4   18     0   21   14      0    0    0     0
## Odoribacter              4      2   10     7    4    0      0    0    3     0
## Oribacterium             2      1    0     2    0    0      7    2    0     0
## Oscillibacter            6     13   34    33   37   10      6    6   29     2
## Parabacteroides          4      8   10     7    5   19      3    0   44     0
## Paraprevotella           0      0    0     0    0    4      0    0    0     0
## Parasutterella           8      3    3     0   16    0      0    0    0     0
## Phascolarctobacterium   20      0    0     0    0    0      0    0    0     0
## Prevotella               0      7    0     0    0   40    164   17    0     0
## Pseudobutyrivibrio      24      0    0     0    0    0      1    0    0    14
## Roseburia              107     45  187    45   97   18     60  257   29   118
## Ruminococcus            96      7   48     0   17    0      1    0    2     0
## Slackia                  0      0    0     0    0    3      3   49    0     0
## Sporacetigenium         12     14   15     0    0    0      0    3   31     0
## Sporobacter              0      2   14     0    6    2      0    0    0     0
## Streptococcus            0      3    3     0    1   40     10   14    0     3
## Subdoligranulum         25     12   67    45   10    3      7  582   53    18
## Sutterella               0      0    0     0    0    0     14    0    0     0
## Tepidibacter             0      0   12     0    0    0      0    0    0     0
## Turicibacter             0      9    2     0    2    0      0    0    0     1
## Uknown                 369    340  607   465  575  254    783 2049 1223   927
## Veillonella              0      0    1     0    0    0      0    0    0     0
##                       TS50.2 TS50 TS51.2 TS51 TS55.2 TS55 TS56.2 TS56 TS57.2
## Acetanaerobacterium        1    0      0   11      0    0      0    0      0
## Acetivibrio                0    0      0    3      1    0      0    1      0
## Acidaminococcus            0    0     56   20      0    0      0    0      0
## Akkermansia                0    0      0  108      0    0      0    0      0
## Alistipes                  0    0     39  290      0    0      0   16      0
## Allisonella                0    0      0    0      6   10      0    0      0
## Anaerostipes               0    0      1    4      1    1      2    0      0
## Anaerotruncus              0    0      0 1013      0    0      7    0      1
## Asaccharobacter            0    0      0    0      0    0      0    0      0
## Bacteroides              441   25    382  242   1194  464      0 1000      1
## Barnesiella                0    0     20   79      0    0      0    0      0
## Bifidobacterium            1    3      2   12      0    1      1    0      0
## Butyricicoccus             1    3      1    7      4    4      0    5      0
## Butyricimonas              0    0      0    0      0    0      0    4      0
## Catenibacterium            1    0      0    0      2    0      0    0      0
## Clostridium                0   10      0    0      0    0      5    7      1
## Collinsella              246  169     11  335      0    0      0    0      0
## Coprobacillus            471  243      3   82      2    0      4   57      0
## Coprococcus                3  125      5  197     13   37      0   25      9
## Desulfovibrio              0    0      0    5      0    0      0    3      0
## Dialister                 27  177     37   59      0    0      5   47      0
## Dorea                     40  361     21  365    103   89      0   45      4
## Eggerthella                5    1      0   16      5    6      0    0      0
## Eubacterium                0    2      4    3      7    3     12    1      0
## Faecalibacterium          59   61     71   39    165  113      0  429      0
## Gordonibacter              0    0      0    0      0    0      0    0      0
## Hespellia                 20   84     24  272     24   32      2    6      0
## Holdemania                 0    0      0    8      1    1      0    0      0
## Lachnobacterium            0    0      0    0      0    0      0    0      0
## Lactobacillus              0    3      8    5      0    0      2    0      0
## Lactococcus                0   15      0    1      2    9      0    0      0
## Lactonifactor              4    4      1    6     41   31      1   16      1
## Marvinbryantia             0    0      0    0      0    0      0    0      0
## Megasphaera                1    0    206   98      0    0     29    0      0
## Odoribacter                0    0      5    1      0    0      0    4      0
## Oribacterium               0    0      0    1      0    0      0    0      0
## Oscillibacter              7    1     16  267     10    0      0    7      1
## Parabacteroides            0    0      3   13      7    1      1    0      0
## Paraprevotella             0    0      5    0    165    7      0    0      0
## Parasutterella            15    1      3    0      0    0    224   40      0
## Phascolarctobacterium      0    0      0    7      0    0      0    0      0
## Prevotella                 0    0    356  206      1    1      4    0      0
## Pseudobutyrivibrio         1  286      0    0      0    0      0   61      0
## Roseburia                 14   13     35  164    152   51    111  134    124
## Ruminococcus               0   74      0    3      0    0      0   36      5
## Slackia                    0    0      2    4      0    0      0    0      0
## Sporacetigenium            0    0      0   19      6   18      5    0      0
## Sporobacter                0    0      0    0      0    0      0    0      0
## Streptococcus              8   31      7   19      2    3      2    1      9
## Subdoligranulum            0    0     21  204      9    3      0    4      0
## Sutterella                 0    0      0    0      0    0      0    0      0
## Tepidibacter               0    0      0    1      0    0      3    0      0
## Turicibacter               0    4      0   30      2    1      2    0      0
## Uknown                  1180 3887    383 1447   1090  914   1938  706   1014
## Veillonella                2    1      4    0     88   68     11    0      0
##                       TS57  TS5 TS6.2 TS61.2 TS61 TS62.2 TS62 TS63.2 TS63
## Acetanaerobacterium      0    4     1      1    0      0    0      0    0
## Acetivibrio              0    0     0      2    0      0    0      0    1
## Acidaminococcus         36    0     0      0    0      0    0      0    0
## Akkermansia              0    3     0      0    0      0    0      1    5
## Alistipes                0    3    16      5    1     15   23     23  119
## Allisonella              0    0     0      0    0      0    0      0    0
## Anaerostipes             9    0     0      0    0      1    0      3    0
## Anaerotruncus            0    0     1      3    0      3    0      3    1
## Asaccharobacter          0    0     1      0    0      1    0      3    0
## Bacteroides            923  539   121    886   81    173  842     80  329
## Barnesiella              0    1     0      0    0      0    0      0    0
## Bifidobacterium          0    0     0      0    0      0    0      0    0
## Butyricicoccus           0    6     3      0    0      0    1      1    2
## Butyricimonas            0    0     5      0    0      0    0      0    0
## Catenibacterium          0    0     0      0    0      0    0      0    0
## Clostridium              4   67     0      0    0      0    1      1    3
## Collinsella              0  215    55     13    5      0    2      0    0
## Coprobacillus            0  723     0      6    1      3    1      1   46
## Coprococcus             12  103     5     45    4     19    7     32   54
## Desulfovibrio            0    0     0      1    0      0    0      1    0
## Dialister               13    0     0      0    0      0    0      0    0
## Dorea                   40  347    40     28   10     51   11     25   17
## Eggerthella              1   16     1      0    0     14    2      0    5
## Eubacterium              1    8     0      5    2      6    0      1    2
## Faecalibacterium        30  403   252    116    7     96   78     87   60
## Gordonibacter            0    0     0      0    0      1    0      0    0
## Hespellia               11  114    28     27    0     80   28     12   13
## Holdemania               0    1     0      2    1      1    1      2    2
## Lachnobacterium          0    0     0      1    2      0    0      0    0
## Lactobacillus            0    0     1      0    0      1    0      0    0
## Lactococcus              0    3     0      2    0      1    0      0    0
## Lactonifactor            0   36    11      0    2      2   10      0   21
## Marvinbryantia           0    0     1      0    1      0    0      0    0
## Megasphaera              0    0     0      0    0      0    0      0    0
## Odoribacter              0    0     1      6    3      3    3      0    0
## Oribacterium             0    0    10      0    0      0    0      0    0
## Oscillibacter            8   18     3     71    3     14    8     18   22
## Parabacteroides          8    2     2      3    0      2    9      0    0
## Paraprevotella           0    0     8      0    0      1    0      0    0
## Parasutterella           0    0     0     11    1      0    1      0    0
## Phascolarctobacterium    0    0     0     39    0     10   12      3    5
## Prevotella               0    2   653      0    0      5    0      0    0
## Pseudobutyrivibrio       0    8     0      0    0      3    0      0    0
## Roseburia              202  232    32     19    9    148   42     43   71
## Ruminococcus             4    8     3    144   24      2    2     39   15
## Slackia                  0    0     0      0    0      0    0      0    0
## Sporacetigenium          0    1     0     15    0     47    0      0    0
## Sporobacter              0    0     0      0    0      0    0      0    0
## Streptococcus            5    3    18      4    2     12    1      7    0
## Subdoligranulum          0  135    10     20    0      7    0      3    1
## Sutterella              49    0     0      1    0      0    0      0    0
## Tepidibacter             0    0     7      0    0      1    0      0    0
## Turicibacter             0    0     0      4    0      2    0      0    0
## Uknown                3536 4600   621    762  177   1074  180    691  382
## Veillonella              2    0     0      0    0      0    0      0    0
##                       TS65.2 TS65 TS66.2 TS66 TS67.2 TS67 TS68.2 TS68 TS69.2
## Acetanaerobacterium        0    1      1    6      0    0      3    0      0
## Acetivibrio                0    0      0    9      0    0      0    0      0
## Acidaminococcus            0    0      0    0      0    0      0    0      0
## Akkermansia                9    0      1   48      0    0      0    0      0
## Alistipes                156   46    186  283      9    9      9   21     13
## Allisonella                0    0      0    0      0    0      0    0      0
## Anaerostipes               0    0      1    0      0    0      0    0      0
## Anaerotruncus             10    3      4    1      1    0      9    0      0
## Asaccharobacter            0    3      0    1      0    2      0    2      2
## Bacteroides              301  188   1227  372    366  122    259  698    461
## Barnesiella               59   13     39   23      3    2      5    3      4
## Bifidobacterium            0    0      0    0      0    0      0    0      1
## Butyricicoccus             0    3      9    0      3    0      1    8      3
## Butyricimonas              0    0     11    8      0    0      1    0      0
## Catenibacterium            0    0      0    0      0    0      1    0      0
## Clostridium               34    6      1    1      0    1      1   10      1
## Collinsella               16    3     11    0      0    4     15   16      4
## Coprobacillus             14   20     65   13      2   24      7   27     16
## Coprococcus               17   19     41   17     34    0     18   35     22
## Desulfovibrio              4    0     10    0      0    0      0    0      0
## Dialister                  3    5     62    0      0    0      9    9      0
## Dorea                      7   29     43    5     42    1     23   20     14
## Eggerthella                0    1      2    0      1    0      0    0      0
## Eubacterium                0    1      0    1      0    0     28   31      0
## Faecalibacterium         180  119    467   68    208   66    372  364    228
## Gordonibacter              0    0      2    0      0    0      0    0      0
## Hespellia                  3    1      1   10      2    3      1   17      0
## Holdemania                 0    0      3    0      1    0      1    0      0
## Lachnobacterium            8    7     11    4      0    0      0    0      0
## Lactobacillus              0    0     14    0      0    0      0    0      0
## Lactococcus                0    0      0    0      0    0      0    1      0
## Lactonifactor              0    2     31    4      3    1      4    5      0
## Marvinbryantia             0    5      0    0      0    0      0    0      1
## Megasphaera                4    0     89    0      0    0      0    0      0
## Odoribacter                8    2     18    9      3    2      3    1      3
## Oribacterium               0    0      0    0      0    0      0    0      5
## Oscillibacter             96   26    124  274      9   13     23   19     32
## Parabacteroides            6    7     68    3      0    2      1   12      2
## Paraprevotella             0    0     11    0      0    0      0    0      1
## Parasutterella             3    7      2    0      7    4      0    1     12
## Phascolarctobacterium      3    4      1    0     15    8      0    0     16
## Prevotella                 1    0      1    0      0    0     37  131      1
## Pseudobutyrivibrio         0    1      1    0      0    0      0    0      0
## Roseburia                 23   71    133   25     65   12     88  154     72
## Ruminococcus              19   26    152   30     47   20     68   37     42
## Slackia                    0    0      0    0      0    0      0    0      0
## Sporacetigenium            0    0      0    0      0    1      0    0      0
## Sporobacter                0    0      0    2      0    0      0    0      0
## Streptococcus              2    4      0    1      0    2      1    1      1
## Subdoligranulum           63   68     85   38      8    1     15   12      5
## Sutterella                 0    0     89    0      0    0      6    0      0
## Tepidibacter               4    0      0    0      1    0      7    0      1
## Turicibacter               8    2      0    0      0    0      0    0      0
## Uknown                  1005  609   2713  799    513  174    322  333    605
## Veillonella                0    0      0    0      0    0      0    0      0
##                       TS69 TS6 TS7.2 TS70.2 TS70 TS71.2 TS71 TS72 TS73 TS74.2
## Acetanaerobacterium      0   1     1      1    2      2    0    1    0      0
## Acetivibrio              2   0     0      0    0      0    1    0    4      0
## Acidaminococcus          0   0     0      0    2      0    0    0    0      0
## Akkermansia             14   6     0      0    0      0    0    0   10      1
## Alistipes              161  36     6     19   26     22   24   20   81     35
## Allisonella              0   0     0      0    0      0    0    0    3      0
## Anaerostipes             0   1     0      0    0      0    0    0    1      0
## Anaerotruncus           12   0     0      0    4      2    0    1   14     14
## Asaccharobacter          0   0     2      0    1      0    0    2    1      0
## Bacteroides            808 835    47    218  253    448  532  285  945     42
## Barnesiella              0   2     0      8   13      1    0    0    0      4
## Bifidobacterium          0   1     0      0    0      0    0    0    0      0
## Butyricicoccus           1   3     0      1    2      5    1    2    5      2
## Butyricimonas            0   7     1      0    0      0    0    0    0      1
## Catenibacterium          0   0     0      0    0      0    0    0    0      0
## Clostridium             14   4     8      1    3      1    2    0    0      3
## Collinsella             40  14   185      1    2      2    5    1    7     20
## Coprobacillus          161  92     7      0    2      7   23    8    2      9
## Coprococcus             11  15    60     49   28     10   27   13    0     29
## Desulfovibrio            0   0     0      0    0      0    1    0    0      0
## Dialister                1   0     0      0    0     11   19    0    0      0
## Dorea                    2  54    55      7   11     12   79    4   28     24
## Eggerthella              3   1     0      0    0      0    0    0    1      1
## Eubacterium              0   1   327      0    0      2    0    0   30     34
## Faecalibacterium       140 828   369    343  137    318   93   71  188     81
## Gordonibacter            0   1     0      0    0      0    0    0    0      0
## Hespellia                0  96    37      0    0      2    8    0   36      1
## Holdemania               1   3     0      0    0      1    0    1    0      1
## Lachnobacterium          0   0     0      0    0      0    3    0    0      0
## Lactobacillus           93   6     0      0    0      0    2    0    0      0
## Lactococcus              0   0     0      0    0      0    0    0    0      1
## Lactonifactor            0  42    12      2    0      4   17    0   15      1
## Marvinbryantia           0   0     0      0    3      0    1    0    0      0
## Megasphaera              0   0     0      0    0      0    0    0    0      0
## Odoribacter              1   3     0      6   10      6    2    1    0      0
## Oribacterium             1   0     0      0    0      2    5    0    0      0
## Oscillibacter           23  36     5     33   39     30   18   13   13     46
## Parabacteroides        119   6     0      0    1      1    1    0    1      0
## Paraprevotella           0  33     0      0    0      0    0    0    0      0
## Parasutterella           0   0     0      0    0      3    1    0    0      0
## Phascolarctobacterium    0   0     1      4   16      1    9   15   47      0
## Prevotella               0 756    18      0    0      0    0    0    0     19
## Pseudobutyrivibrio       0   0     0      0    0      0    1   20    2      0
## Roseburia               71 133    42     29   12     34   16   16   71     53
## Ruminococcus             0 122    42     86   47     57   39   91   40     53
## Slackia                  2   0     0      0    0      0    0    0    0      2
## Sporacetigenium          0   0     1      5    2      0   34    0    0      0
## Sporobacter              0  17     0      0    1      0    0    0    0     33
## Streptococcus          433   9     5      3    0      6   12    9    2      0
## Subdoligranulum          1  24    12    128   60     33    6   25    0     58
## Sutterella               8   0     1      0    0      0    0    0    0      0
## Tepidibacter             0   0     3      1    2      0   22    0    0      0
## Turicibacter             0   1    13      3    0      0   24    0    0      0
## Uknown                 680 879   523    414  364    431  615  335  607    398
## Veillonella              0   0     0      0    0      0    4    0    0      0
##                       TS74 TS75.2 TS75 TS76.2 TS76 TS78.2 TS78  TS7 TS8.2
## Acetanaerobacterium      0      0    0      0    0      0    0    0     0
## Acetivibrio              0      0    0      0    0      0    0    1     0
## Acidaminococcus          0      0    0      0    0      0    3    0     0
## Akkermansia              0      0    0      0    0      0    0    0     0
## Alistipes               29      1    1     11   18      3    4  119     9
## Allisonella              0      0    0      0    0      7    0    0     0
## Anaerostipes             0      1    1      0    0      1    0    0     0
## Anaerotruncus            3      2    0     10   30      0    0    4     0
## Asaccharobacter          2      0    0      1    0      0    1    5     0
## Bacteroides            170    278  234    332  857    539   95  864    30
## Barnesiella              4      0    0      4    1      0    0    2    10
## Bifidobacterium          0      0    0      1    0      0    0    0     0
## Butyricicoccus           1      0    0      5    1      3    1   12     0
## Butyricimonas            1      0    0      0    0      0    1   11     3
## Catenibacterium          0      0    0      0    0     47    8    0     0
## Clostridium              1      0    0     11    2      0    0  127     0
## Collinsella             10      0    0     20    5      8   17   29     0
## Coprobacillus            2      0    1      2    8      1    3   36     7
## Coprococcus             17     63   47     48   20     38    5   96    10
## Desulfovibrio            0      0    0      0    0      0    0    2     5
## Dialister                0      0    0      4   19      0    0    0     1
## Dorea                   74      0    0     64   14     31   15   90     1
## Eggerthella              2      4    2      0    0      0    0    0     1
## Eubacterium             30      3    5    347    0     75   17  596     8
## Faecalibacterium       280      0    0    369  439    260  172  821    69
## Gordonibacter            0      0    0      1    0      0    0    0     0
## Hespellia               20     75    3      0    2      0    0   61     2
## Holdemania               2      0    0      0    0      0    0    0     0
## Lachnobacterium          0      0    0      5    0      0    0    0     0
## Lactobacillus            0      0    0      9    7      1    3    0     0
## Lactococcus              0      0    0     12    0      0    0    0     0
## Lactonifactor            3      0    0      2    3      0    0    8     1
## Marvinbryantia           0      0    0      0    0      0    2    0     0
## Megasphaera              0     28    2     17    3     28   19    0     0
## Odoribacter              4      4    0      1    5      0    0    4     0
## Oribacterium             0      0    0      0    2      0    0    0     4
## Oscillibacter           10      3    2     32   33      5    4  133    19
## Parabacteroides          0      0    0      0    0      0    0   13     9
## Paraprevotella           4      0    0     23    0     47    2    4     0
## Parasutterella           2     13   13      0    0      0    0    0     0
## Phascolarctobacterium    0     10    5      0    0      0    0  132     0
## Prevotella              17      0    0      0    0      0   55  668     0
## Pseudobutyrivibrio       0      0    0      1    0      0    6    0     1
## Roseburia               50    419  234    101  129     53   20  300    10
## Ruminococcus            59      2    1     27   15      0    4  144     4
## Slackia                  0      0    0      0    0      0    1    0     0
## Sporacetigenium          7      0    0      0    0      0    0   56     1
## Sporobacter              1      0    0      0    0      0    0    0     0
## Streptococcus            7      0    1     14   33     19    0    2     0
## Subdoligranulum         21      0    0     74   46     10    4   47    12
## Sutterella               0      0    0      1    6      1    7   17     4
## Tepidibacter             0      0    0      0    5      0    0    4     0
## Turicibacter             1      0    0      0    0      0    7    1     0
## Uknown                 516    343  368    533  287    370  121 1399    80
## Veillonella              0      0    0      1    3     11    1    0     0
##                       TS82.2 TS82 TS83.2 TS83 TS84.2 TS84 TS86.2 TS86 TS87.2
## Acetanaerobacterium        2    0      0    0      0    0      0    4      1
## Acetivibrio                0    0      0    0      0    0      1    0      0
## Acidaminococcus            0    0      0    0      2   23      0    1      1
## Akkermansia                0    0      0    0      0    0      0    0      0
## Alistipes                 23    0      3    1      2    0     12   11      4
## Allisonella                0    0      0    0      0    0      0    0      1
## Anaerostipes               0    0      0    0      0    0      0    0      0
## Anaerotruncus              0    0      0    0      6    1      1    5      5
## Asaccharobacter            0    0      0    0      0    0      0    0      0
## Bacteroides               59   28    565  232     12   15    304  484     27
## Barnesiella                0    0      0    0      1    1      3   10      0
## Bifidobacterium            0    0      0    1      0    0      0    0      0
## Butyricicoccus             5    4      6    2      1    1      1    4      0
## Butyricimonas              0    0      0    0      1    0      0    1      0
## Catenibacterium            0    0      0    0    164   28      0    0    110
## Clostridium                0    1      4    0      0    2      1    1      0
## Collinsella                8    1     28   17      2    2      0   37      5
## Coprobacillus              7   11     12    3      0    0      0    4      0
## Coprococcus               12   10     18    4     24   10      2   11     16
## Desulfovibrio              0    0      6    1      1    7      0    0      0
## Dialister                  0    0      0    0      0    0      0    1      0
## Dorea                     12   18     51   15     43    0     20   20     39
## Eggerthella                2    0      0    0     12    1      0    0      2
## Eubacterium                0    0      0    0     34   11     10   23     17
## Faecalibacterium         484  255    410  297     55  271    333  586     91
## Gordonibacter              0    0      0    0      0    0      0    0      0
## Hespellia                  0    1      0    0      4    0      0    3      0
## Holdemania                 0    0      0    0      0    0      0    0      0
## Lachnobacterium            0    0      0    0     11    0      0    8      0
## Lactobacillus              0    0      0    0      3    4      0    0      0
## Lactococcus                0    0      0    0      0    0      0    0      0
## Lactonifactor              2    6      6    3      0    1      0   22      2
## Marvinbryantia             0    0      0    0      0    0      0    0      0
## Megasphaera                0    0      0    0      6   29      3   55      0
## Odoribacter                2    0      0    0      0    0      1    4      3
## Oribacterium               3    0      0    0      0    1      2    0      0
## Oscillibacter             29    4      7    4     14    5     24   68     66
## Parabacteroides            2    1      2    1      0    1      0    1      2
## Paraprevotella             0    0      0    0     18   46      2   30      3
## Parasutterella             1    1      0    0      0    0      2   10      0
## Phascolarctobacterium      0    0     25    4      0    0      0    0      0
## Prevotella                33  109     37   11     59  106     51  119     24
## Pseudobutyrivibrio         0    0      0    0      0    0      0    0      0
## Roseburia                 36   22     81   60     40   22     19   79     21
## Ruminococcus              33    3     11   29      5    5     68  154     12
## Slackia                    0    0      1    0      0    0      0    0      0
## Sporacetigenium            0    0      0    0      2    6      0   10      0
## Sporobacter                0    0      0    0      1    0      2   10      0
## Streptococcus              3    4      1    0      6    0      0    1      3
## Subdoligranulum           18    2     55   30     31   17     20  254     29
## Sutterella                 0    0      3    1      0    2      0    0      0
## Tepidibacter               0    0      0    0      0    0      2    0      1
## Turicibacter               0    1      0    0      0    0      0    0      0
## Uknown                   144  224    642  303    317  186    156  638    494
## Veillonella                0    0      2    1      0    0      0    0      0
##                       TS87 TS88.2 TS88 TS89.2 TS89  TS8 TS9.2 TS90.2 TS90
## Acetanaerobacterium      0      0    0      0    0    7     1      0    1
## Acetivibrio              0      0    0      0    0    1     0      0    0
## Acidaminococcus          2      7   16      0    0    0     0      0    1
## Akkermansia              0      0    0      0    0    6     0      0    6
## Alistipes                8     13   13     44   14  421    45     83   43
## Allisonella              1      0    0      0    0    1     0      0    0
## Anaerostipes             0      0    0      0    0    0     0      0    0
## Anaerotruncus           40      1    2      7    0   16     0      4   24
## Asaccharobacter          0      1    0      0    0    1     1      0    0
## Bacteroides             89    641  919    610  361 1858   470    478  576
## Barnesiella              4      4    6     12    8   72     0     23   43
## Bifidobacterium          0      0    0      0    0    1     0      0    0
## Butyricicoccus           1      1    1     10    2   10     1      2    2
## Butyricimonas            0      0    0     10    8   16     0      6    8
## Catenibacterium        264     12   30      0    0    0     0      0    0
## Clostridium              1      1    0      1    0   20     0      0    0
## Collinsella             14      9    6      3    9    0    15      4    3
## Coprobacillus            0      3    1      0    5   53    10      0    1
## Coprococcus             36      4   10     34   11  128   115     16   87
## Desulfovibrio            0      0    3      1    0   31     0      0    0
## Dialister                0      0    0     17    0  103     7     19   10
## Dorea                   28      5   22     26   27  107     9     25   25
## Eggerthella              0      1    2      0    0    0     1      2    0
## Eubacterium             30      0    4     37   23  119     0      0    0
## Faecalibacterium        70    179   18    657  155 1995   207     93  139
## Gordonibacter            0      0    0      0    0    0     0      0    0
## Hespellia                0      4   11      0    1    0    29     20   14
## Holdemania               0      0    0      0    0    0     1      0    1
## Lachnobacterium          0      0    0      0    3    0     0      0    0
## Lactobacillus            0      0    4      0    0    0     0      0    0
## Lactococcus              0      0    2      0    0    0     0      0    0
## Lactonifactor            0      2    0      5    4    5     0      1    1
## Marvinbryantia           0      0    0      0    0    0     0      0    0
## Megasphaera              0      0    0      6   12    0     0      4    2
## Odoribacter              9      0    0     10   12   13     0      0    0
## Oribacterium             0      4    0      0    0    0     4      0    1
## Oscillibacter          143     41    5     65   37  426    30     39  109
## Parabacteroides          4      2    1      1    0   57     1      1    2
## Paraprevotella           2      0    0      6    0    0     0      0    0
## Parasutterella           0      8    4      0    0    0     5      0    1
## Phascolarctobacterium    0      9   12     12   15    0     0     12    4
## Prevotella              30      0    0      0    0    1     0      0    0
## Pseudobutyrivibrio       0      0    0      0    0    0     0      0    0
## Roseburia               82     65   83     86   33  132   150     50   28
## Ruminococcus            35     56    0     53   24  534     7    111   50
## Slackia                  1      0    0      0    2    0     0      0    0
## Sporacetigenium          0      0    0      0    0   40     0      0    1
## Sporobacter              0      0    0      6    0   15     2      0    0
## Streptococcus           15      1    1      0    1    2    20     58   47
## Subdoligranulum         61     12    0    106  167  542    31     35   72
## Sutterella               0      0    0     13   16   30     0      2    1
## Tepidibacter             0      0    0      0    1   35     0      0    1
## Turicibacter             0      0    0      0    0   11     4      0    0
## Uknown                 577    232  211    532  381 2143   312    654  356
## Veillonella              0      0    0      0    0    0     0      0    0
##                       TS91.2 TS91 TS92.2 TS92 TS94.2 TS94 TS95.2 TS95 TS96.2
## Acetanaerobacterium        0    0      0    0      1    1      0    2      0
## Acetivibrio                0    0      0    0      1    0      0    0      2
## Acidaminococcus            0    0      0    0      0    0      0    0     31
## Akkermansia                0    0      0    0      2    3      4   10      0
## Alistipes                  2    0      0    0     45   26     38   85     56
## Allisonella                6    1      0    0      0    0      0    0      0
## Anaerostipes               0    0      0    0      0    0      0    0      0
## Anaerotruncus             10    2      0    0      3    1      0  111      9
## Asaccharobacter            0    0      0    0      1    1      0    0      1
## Bacteroides              121    8    110  154    582  481    672  694    358
## Barnesiella                0    0      0    0      0    0      0    0     13
## Bifidobacterium            1    0      1    0      0    0      1    2      2
## Butyricicoccus             4    0      3    4      3    2      6    3      0
## Butyricimonas              0    0      1    1      0    0      0    0      1
## Catenibacterium          181    9      0    0     10    2     25   87      0
## Clostridium                1    0      0    0      3   10      1    0      8
## Collinsella               31    5     59    8     24   15      0    0     23
## Coprobacillus              5    2     26   15      0    1      0    0      0
## Coprococcus               26    9     26   15     21   17      4    6     25
## Desulfovibrio              0    0      0    0      0    0      0    0      0
## Dialister                  0    0     19    7     23   12     20   23      0
## Dorea                     69   22     96   21     10    0      6    7      7
## Eggerthella                0    0      3    0      0    2      0    1      0
## Eubacterium               61    8      1    0     27   12     20   77      2
## Faecalibacterium         100   51    390  455     70  143    274  181    292
## Gordonibacter              0    0      0    0      0    1      0    0      0
## Hespellia                  3    1      1    2      2    0      1    0      1
## Holdemania                 0    1      1    0      0    0      1    1      0
## Lachnobacterium            0    0      0    0      0    0      0    0      6
## Lactobacillus              0    0      0    0      2   12      3    7     32
## Lactococcus                0    1      0    0      0    0      0    0      0
## Lactonifactor             13    0      8    0      0    0      0    0      3
## Marvinbryantia             0    0      0    0      0    0      0    0      1
## Megasphaera                0    0     18   22     47   36      0    0     11
## Odoribacter                0    0      0    0      0    0      9    9     13
## Oribacterium               1    0      0    0      0    0      0    0      8
## Oscillibacter              8    0      3    0    101   33     49   93    122
## Parabacteroides            0    0      1    1      5    2      0    0      1
## Paraprevotella             6    4      0    0     16    0      0    0      0
## Parasutterella             0    0      0    0      0    0     17    6      6
## Phascolarctobacterium      0    0      0    0      0    0      0    0      0
## Prevotella               232    6      1    0      0    0      0    0      7
## Pseudobutyrivibrio         0    0      1    0      1    0      0    0      1
## Roseburia                 69   13     43   12     26   24     68   45     34
## Ruminococcus              24   13      2    0     87   33     61   69     48
## Slackia                    0    1      0    0      0    0      0    0      0
## Sporacetigenium            1    0     13    2      9    3      1    3     50
## Sporobacter                0    0      0    0      0    1      0    0     16
## Streptococcus              9    5      6    0      0    2      7    5      6
## Subdoligranulum           40   16     25   19     51  100    126  122     41
## Sutterella                 0    0      0    0      2    0      0    0      0
## Tepidibacter               0    0      0    1      0    0      2    0      5
## Turicibacter               0    0      0    0      1    0      0    0      2
## Uknown                   489  190    748  589    520  372    332  364    784
## Veillonella                0    0      0    0      0    0      0    0      0
##                       TS96 TS97.2 TS97 TS98.2 TS98
## Acetanaerobacterium      0      0    0      0    0
## Acetivibrio              0      0    0      0    0
## Acidaminococcus         15     10    6      0   11
## Akkermansia              1      0    0      0    0
## Alistipes               15      3    1      1   15
## Allisonella              0      0    0      0    0
## Anaerostipes             0      0    0      0    0
## Anaerotruncus            0      0    0      0    8
## Asaccharobacter          0      1    2      0    0
## Bacteroides            219    436  550    861  835
## Barnesiella              3      0    0      0    0
## Bifidobacterium          0      0    0      0    0
## Butyricicoccus           0      1    6      7    2
## Butyricimonas            1      0    0      0    0
## Catenibacterium          0      0    0      0    0
## Clostridium              0      3    4     11    2
## Collinsella              4     15   12      0    0
## Coprobacillus            1     41   15      8    0
## Coprococcus             17      6    2     13    0
## Desulfovibrio            0      0    1      2    0
## Dialister                0      4    2     16    0
## Dorea                    2     69   21     15   35
## Eggerthella              1      1    0      0    0
## Eubacterium              0      1    0      0    3
## Faecalibacterium       202    110  264    252  573
## Gordonibacter            0      0    0      0    0
## Hespellia                0     51   22     18    3
## Holdemania               1      0    0      0    0
## Lachnobacterium          0      0    0      0    0
## Lactobacillus           15      0    0      0    0
## Lactococcus              0      0    0      0    0
## Lactonifactor            0      1    5      1    0
## Marvinbryantia           0      0    0      0    0
## Megasphaera             22      0    0      0   28
## Odoribacter              4      0    0      0    3
## Oribacterium             2      0    0      0    0
## Oscillibacter           39      0    0      2    5
## Parabacteroides          1      0    1      0    0
## Paraprevotella           0      6    7      0    0
## Parasutterella           2      0    0      4    3
## Phascolarctobacterium    0      3    0      0    0
## Prevotella              22      0    0      0    0
## Pseudobutyrivibrio       0      0  138      0    0
## Roseburia               40     28  151     51   25
## Ruminococcus            13      0   15      1   50
## Slackia                  0      0    0      0    0
## Sporacetigenium          1      0    0      1   13
## Sporobacter              0      0    0      0    0
## Streptococcus           25      3    6      2    1
## Subdoligranulum         46      0   11     11    0
## Sutterella               0     19   30      0    0
## Tepidibacter             0      0    0      0    0
## Turicibacter             0      2    5      0    0
## Uknown                 287    634  270    199  424
## Veillonella              0      0    1      1    0
```

```r
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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-22-1.pdf)<!-- --> 

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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-25-1.pdf)<!-- --> 

Seemingly the covariate "Age_Years" does not have effect on the model as "Delivery_Mode" would,
and "Sex" to some extent. Let's take a closer look at the two latter ones:


```r
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-26-1.pdf)<!-- --> 

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.2.0 (2022-04-22)
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
 [1] fido_1.0.2                     ANCOMBC_1.99.1                
 [3] forcats_0.5.1                  stringr_1.4.0                 
 [5] dplyr_1.0.9                    purrr_0.3.4                   
 [7] readr_2.1.2                    tidyr_1.2.0                   
 [9] tibble_3.1.7                   ggplot2_3.3.6                 
[11] tidyverse_1.3.1                knitr_1.39                    
[13] MicrobiomeStat_1.1             Maaslin2_1.10.0               
[15] ALDEx2_1.28.1                  zCompositions_1.4.0-1         
[17] truncnorm_1.0-8                NADA_1.6-1.1                  
[19] survival_3.3-1                 MASS_7.3-57                   
[21] tidySummarizedExperiment_1.6.1 patchwork_1.1.1               
[23] mia_1.3.31                     MultiAssayExperiment_1.22.0   
[25] TreeSummarizedExperiment_2.1.4 Biostrings_2.64.0             
[27] XVector_0.36.0                 SingleCellExperiment_1.18.0   
[29] SummarizedExperiment_1.26.1    Biobase_2.56.0                
[31] GenomicRanges_1.48.0           GenomeInfoDb_1.32.2           
[33] IRanges_2.30.0                 S4Vectors_0.34.0              
[35] BiocGenerics_0.42.0            MatrixGenerics_1.8.1          
[37] matrixStats_0.62.0-9000        BiocStyle_2.24.0              
[39] rebook_1.6.0                  

loaded via a namespace (and not attached):
  [1] coda_0.19-4                 bit64_4.0.5                
  [3] irlba_2.3.5                 DelayedArray_0.22.0        
  [5] data.table_1.14.2           rpart_4.1.16               
  [7] doParallel_1.0.17           RCurl_1.98-1.7             
  [9] generics_0.1.3              ScaledMatrix_1.4.0         
 [11] microbiome_1.18.0           timeSeries_3062.100        
 [13] RSQLite_2.2.14              proxy_0.4-27               
 [15] bit_4.0.4                   tzdb_0.3.0                 
 [17] xml2_1.3.3                  lubridate_1.8.0            
 [19] assertthat_0.2.1            DirichletMultinomial_1.38.0
 [21] viridis_0.6.2               xfun_0.31                  
 [23] fBasics_3042.89.2           ggdist_3.1.1               
 [25] hms_1.1.1                   evaluate_0.15              
 [27] DEoptimR_1.0-11             fansi_1.0.3                
 [29] dbplyr_2.2.1                readxl_1.4.0               
 [31] igraph_1.3.2                DBI_1.1.3                  
 [33] htmlwidgets_1.5.4           tensorA_0.36.2             
 [35] hash_2.2.6.2                ellipsis_0.3.2             
 [37] energy_1.7-10               backports_1.4.1            
 [39] bookdown_0.27               permute_0.9-7              
 [41] deldir_1.0-6                sparseMatrixStats_1.8.0    
 [43] vctrs_0.4.1                 abind_1.4-5                
 [45] tidybayes_3.0.2             cachem_1.0.6               
 [47] withr_2.5.0                 robustbase_0.95-0          
 [49] checkmate_2.1.0             vegan_2.6-2                
 [51] treeio_1.20.1               getopt_1.20.3              
 [53] cluster_2.1.3               gsl_2.1-7.1                
 [55] ape_5.6-2                   dir.expiry_1.4.0           
 [57] lazyeval_0.2.2              crayon_1.5.1               
 [59] labeling_0.4.2              pkgconfig_2.0.3            
 [61] nlme_3.1-158                vipor_0.4.5                
 [63] nnet_7.3-17                 rlang_1.0.4                
 [65] spatial_7.3-15              lifecycle_1.0.1            
 [67] filelock_1.0.2              phyloseq_1.40.0            
 [69] modelr_0.1.8                rsvd_1.0.5                 
 [71] distributional_0.3.0        cellranger_1.1.0           
 [73] rngtools_1.5.2              graph_1.74.0               
 [75] Matrix_1.4-1                lpsymphony_1.24.0          
 [77] Rhdf5lib_1.18.2             boot_1.3-28                
 [79] base64enc_0.1-3             reprex_2.0.1               
 [81] beeswarm_0.4.0              png_0.1-7                  
 [83] viridisLite_0.4.0           stabledist_0.7-1           
 [85] rootSolve_1.8.2.3           bitops_1.0-7               
 [87] rhdf5filters_1.8.0          blob_1.2.3                 
 [89] DelayedMatrixStats_1.18.0   doRNG_1.8.2                
 [91] decontam_1.16.0             jpeg_0.1-9                 
 [93] DECIPHER_2.24.0             beachmat_2.12.0            
 [95] scales_1.2.0                memoise_2.0.1              
 [97] magrittr_2.0.3              plyr_1.8.7                 
 [99] zlibbioc_1.42.0             compiler_4.2.0             
[101] RColorBrewer_1.1-3          clue_0.3-61                
[103] lme4_1.1-30                 cli_3.3.0                  
[105] ade4_1.7-19                 lmerTest_3.1-3             
[107] pbapply_1.5-0               htmlTable_2.4.1            
[109] Formula_1.2-4               mgcv_1.8-40                
[111] tidyselect_1.1.2            stringi_1.7.8              
[113] highr_0.9                   yaml_2.3.5                 
[115] svUnit_1.0.6                BiocSingular_1.12.0        
[117] latticeExtra_0.6-30         ggrepel_0.9.1              
[119] grid_4.2.0                  tools_4.2.0                
[121] lmom_2.9                    parallel_4.2.0             
[123] rstudioapi_0.13             logging_0.10-108           
[125] foreign_0.8-82              foreach_1.5.2              
[127] statip_0.2.3                optparse_1.7.1             
[129] gridExtra_2.3               gld_2.6.5                  
[131] posterior_1.2.2             farver_2.1.1               
[133] Rtsne_0.16                  stable_1.1.6               
[135] RcppZiggurat_0.1.6          digest_0.6.29              
[137] BiocManager_1.30.18         Rcpp_1.0.9                 
[139] broom_1.0.0                 scuttle_1.6.2              
[141] httr_1.4.3                  Rdpack_2.3.1               
[143] colorspace_2.0-3            rvest_1.0.2                
[145] XML_3.99-0.10               fs_1.5.2                   
[147] modeest_2.4.0               splines_4.2.0              
[149] yulab.utils_0.0.5           rmutil_1.1.9               
[151] statmod_1.4.36              expm_0.999-6               
[153] tidytree_0.3.9              scater_1.24.0              
[155] Exact_3.1                   multtest_2.52.0            
[157] plotly_4.10.0               jsonlite_1.8.0             
[159] nloptr_2.0.3                CodeDepends_0.6.5          
[161] timeDate_3043.102           Rfast_2.0.6                
[163] R6_2.5.1                    Hmisc_4.7-0                
[165] pillar_1.7.0                htmltools_0.5.2            
[167] glue_1.6.2                  fastmap_1.1.0              
[169] minqa_1.2.4                 BiocParallel_1.30.3        
[171] BiocNeighbors_1.14.0        class_7.3-20               
[173] codetools_0.2-18            pcaPP_2.0-2                
[175] mvtnorm_1.1-3               utf8_1.2.2                 
[177] lattice_0.20-45             arrayhelpers_1.1-0         
[179] numDeriv_2016.8-1.1         ggbeeswarm_0.6.0           
[181] DescTools_0.99.45           interp_1.1-3               
[183] biglm_0.9-2.1               rmarkdown_2.14             
[185] biomformat_1.24.0           munsell_0.5.0              
[187] e1071_1.7-11                rhdf5_2.40.0               
[189] GenomeInfoDbData_1.2.8      iterators_1.0.14           
[191] haven_2.5.0                 reshape2_1.4.4             
[193] gtable_0.3.0                rbibutils_2.2.8            
```
</div>
