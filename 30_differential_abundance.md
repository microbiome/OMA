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

Due to the complex data characteristics of microbiome sequencing data,
differential abundance analysis of microbiome data faces many
statistical challenges [@Yang2022], including:
  
- Highly variable. The abundance of a specific taxon could range over
  several orders of magnitude.
  
- Zero-inflated. In a typical microbiome dataset, more than 70% of the
  values are zeros. Zeros could be due to either physical absence
  (structural zeros) or insufficient sampling effort (sampling zeros).
  
- Compositional. Increase or decrease in the (absolute) abundance of
  one taxon at the sampling site will lead to apparent changes in the
  relative abundances of other taxa in the sample.

As summarized in @Yang2022, to address the above statistical
chanllenegs:

- Over-dispersed count models has been proposed to address zero
  inflation, such as the negative binomial model used by edgeR
  [@Chen2016] and DESeq2 [@Love2014], the beta-binomial model used by
  corncorb [@Martin2021].

- Zero-inflated mixture models has aslo been proposed to address zero
  inflation, such as zero-inflated log-normal/normal mixture model
  used by metagenomeSeq [@Paulson2017] and RAIDA [@Sohn2015],
  zero-inflated beta-binomial model used by ZIBB [@ZIBB2018], and
  zero-inflated negative binomial model used by Omnibus
  [@Omnibus2018].

- Bayesian methods have been used to impute the zeros for methods
  working on proportion data, accounting for sampling variability and
  sequencing depth variation. Examples include ALDEx2 [@Gloor2016] and
  eBay [@Liu2020].

- Other methods use the pseudo-count approach to impute the zeros,
  such as MaAsLin2 [@Mallick2020] and ANCOMBC [@ancombc2020].

- Different strategies have been used to address compositional
  effects, including:

  - Robust normalization. For example, trimmed mean of M-values (TMM)
    normalization used by edgeR, relative log expression (RLE)
    normalization used by DESeq2 [@Love2014], cumulative sum scaling
    (CSS) normalization used by metagenomeSeq, centered log-ratio
    transformation (CLR) normalization used by ALDEx2 [@Gloor2016] and
    geometric mean of pairwise ratios (GMPR) normalization used by
    Omnibus [@Omnibus2018]. Wrench normalization [@Kumar2018] corrects
    the compositional bias by an empirical Bayes approach, which has
    been recommended in metagenomeSeq [@Paulson2017] tutorial.
  
  - Reference taxa approach used by DACOMP [@Brill2019] and RAIDA [@Sohn2015].
  
  - Analyzing the pattern of pairwise log ratios, such as ANCOM [@Mandal2015].
  
  - Bias-correction used by ANCOMBC [@ancombc2020].

Some of the popular tools for differential abundance analysis include:

- ALDEx2 [@Gloor2016] 
- ANCOMBC [@ancombc2020]
- corncob [@Martin2021]
- DESeq2 [@Love2014] 
- edgeR [@Chen2016]
- lefser [@Khlebrodova2021]
- MaAsLin2 [@Mallick2020]
- metagenomeSeq [@Paulson2017]
- limma [@Ritchie2015]
- LinDA [@Zhou2022]
- ZicoSeq [@Yang2022]
- LDM [@Hu2020]
- RAIDA [@Sohn2015]
- DACOMP [@Brill2019]
- Omnibus [@Omnibus2018]
- eBay [@Liu2020]
- ZINQ [@Ling2021]
- ANCOM [@Mandal2015]
- fastANCOM [@fastANCOM2022]
- [t-test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test)  
- [Wilcoxon test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test)  


We recommend to have a look at @Nearing2022 who compared all these
methods across 38 different datasets. Because different methods use
different approaches (parametric vs non-parametric, different
normalization techiniques, assumptions etc.), the results may differ
between methods, sometimes substantially as @Nearing2022 pointed
out. More recently @Yang2022 comprehensively evaluated these methods
via a semi-parametric framework and 106 real datasets. @Yang2022 also
concluded that different DA methods can sometimes produce discordant
results, opening to the possibility for cherry-picking tools in favor
of one’s own hypothesis. Therefore it is highly recommended to pick
several methods to assess how robust and potentially reproducible your
findings are with different methods.

In this section we demonstrate the use of four methods that can be
recommended based on recent literature (ANCOM-BC [@ancombc2020],
_ALDEx2_ [@Gloor2016], _Maaslin2_ [@Mallick2020], _LinDA_ [@Zhou2022]
and _ZicoSeq_ [@Yang2022]).

The purpose of this section is to show how to perform DAA in R, not
how to correctly do causal inference. Depending on your experimental
setup and your theory, you must determine how to specify any model
exactly.  E.g., there might be confounding factors that might drive
(the absence of) differences between the shown groups that we ignore
here for simplicity. Or your dataset is repeated sampling design,
matched-pair design or the general longitudianl design.  We will
demonstrate how to include covariates in those models. We picked a
dataset that merely has microbial abundances in a TSE object as well
as a grouping variable in the sample data. We simplify the examples by
only including two of the three groups.



```r
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
as.data.frame(colData(tse0)) %>% 
count(Geographical_location) %>%
  kable()
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

Before we jump to our analyses, we may want to perform some data
manipulation.

Let us here do aggregation to genus level, add relative abundance
assay, and perform prevalence filtering.


```r
tse <- agglomerateByRank(tse0, rank = "genus") %>%
       transformCounts(assay.type = "counts",
                       method = "relabundance",
		       MARGIN = "samples") %>%
       # subset based on the relative abundance assay		       
       subsetByPrevalentTaxa(detection = 0,
                             prevalence = 10/100,
			     assay.type = "relabundance")

# Add also clr abundances
tse <- transformCounts(tse, method="clr", pseudocount=1) # not bale to run
```

Regarding prevalence filtering, @Nearing2022 found that applying a 10%
threshold for the prevalence of the taxa generally resulted in more
robust results. Some tools have builtin arguments for that. By
applying the threshold to our input data, we can make sure it is
applied for all tools. 



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
utilize the `glm` functionality within ALDEx2.

The Benjamini-Hochberg procedure is applied by default to correct for
multiple testing. Below we show a simple example that illustrates the
workflow.



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
# Determines the median clr abundance of the feature in all samples and in
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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-1-1.pdf)<!-- --> 

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
  dplyr::select(genus, we.eBH, wi.eBH, effect, overlap) %>%
  kable()
```


\begin{tabular}{l|r|r|r|r}
\hline
genus & we.eBH & wi.eBH & effect & overlap\\
\hline
Anaerococcus & 0.0540 & 0.0150 & 0.9595 & 0.1546\\
\hline
Calditerricola & 0.0769 & 0.0299 & -0.7162 & 0.1702\\
\hline
Chitinivibrio & 0.1216 & 0.0484 & -0.7700 & 0.1776\\
\hline
Corynebacterium & 0.0280 & 0.0035 & 1.1857 & 0.1037\\
\hline
Desulfosporomusa & 0.0851 & 0.0359 & -0.8604 & 0.1733\\
\hline
Geobacillus & 0.0370 & 0.0081 & -1.0962 & 0.1293\\
\hline
Jeotgalicoccus & 0.0276 & 0.0251 & -0.9052 & 0.1676\\
\hline
Paenibacillus & 0.0837 & 0.0345 & -0.9380 & 0.1932\\
\hline
Virgibacillus & 0.1103 & 0.0443 & -0.8750 & 0.1960\\
\hline
\end{tabular}

### ANCOM-BC

The analysis of composition of microbiomes with bias correction
(ANCOM-BC) [@ancombc2020] is a recently developed method for differential
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
```


```r
# store the FDR adjusted results [test on v2.0.3] 
ancombc_result <- cbind.data.frame(taxid = out$res$taxon,
                       ancombc = as.vector(out$res$q_Geographical_locationPune))
```





```r
# store the FDR adjusted results [test on v1.2.2]
ancombc_result <- out$res %>% dplyr::select(starts_with(c("taxon", "lfc", "q")))
```

The object `out` contains all model output. Again, see the
[documentation of the
function](https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc.html)
under **Value** for details. Our question whether taxa are
differentially abundant can be answered by looking at the `res`
object, which contains dataframes with the coefficients, standard
errors, p-values and q-values. Below we show the first entries of this
dataframe.


```r
kable(head(ancombc_result))
```


\begin{tabular}{l|r|r|r|r}
\hline
taxon & lfc\_(Intercept) & lfc\_Geographical\_locationPune & q\_(Intercept) & q\_Geographical\_locationPune\\
\hline
Abyssicoccus & 0.0397 & -0.0569 & 0.8731 & 0.8468\\
\hline
Acidaminococcus & 0.6872 & -0.9022 & 0.0032 & 0.0004\\
\hline
Acinetobacter & 0.1241 & -0.1671 & 0.8972 & 0.8698\\
\hline
Actinomyces & 0.1345 & -0.1806 & 0.6548 & 0.5551\\
\hline
Actinoplanes & 0.2713 & -0.3593 & 0.2760 & 0.1476\\
\hline
Aerococcus & 0.0234 & -0.0357 & 0.8985 & 0.8698\\
\hline
\end{tabular}



### MaAsLin2 

Let us next illustrate MaAsLin2 [@Mallick2020]. This method is based on
generalized linear models and flexible for different study designs
and covariate structures. For details, check their
[Biobakery tutorial](https://github.com/biobakery/biobakery/wiki/maaslin2).


```r
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
```

Which genera are identified as differentially abundant? (leave out "head" to see all).


```r
kable(head(filter(maaslin2_out$results, qval <= 0.05)))
```


\begin{tabular}{l|l|l|r|r|r|l|r|r|r}
\hline
feature & metadata & value & coef & stderr & pval & name & qval & N & N.not.zero\\
\hline
Fructobacillus & Geographical\_location & Nashik & 0.0080 & 0.0011 & 0 & Geographical\_locationNashik & 0 & 47 & 9\\
\hline
Desulfosporomusa & Geographical\_location & Nashik & 0.0373 & 0.0059 & 0 & Geographical\_locationNashik & 0 & 47 & 13\\
\hline
Geobacillus & Geographical\_location & Nashik & 0.1294 & 0.0207 & 0 & Geographical\_locationNashik & 0 & 47 & 27\\
\hline
Pullulanibacillus & Geographical\_location & Nashik & 0.0395 & 0.0062 & 0 & Geographical\_locationNashik & 0 & 47 & 9\\
\hline
Chitinivibrio & Geographical\_location & Nashik & 0.0274 & 0.0045 & 0 & Geographical\_locationNashik & 0 & 47 & 10\\
\hline
Thermoanaerobacter & Geographical\_location & Nashik & 0.0089 & 0.0015 & 0 & Geographical\_locationNashik & 0 & 47 & 10\\
\hline
\end{tabular}

This will create a folder that is called like in the output specified
above. It contains also figures to visualize difference between
significant taxa.


### LinDA 

Lastly, we cover linear models for differential abundance analysis of
microbiome compositional data (@Zhou2022). This is very similar to
ANCOMBC with few differences: 1) LinDA correct for the compositional
bias differently using the mode of all regression coefficients. 2) it
is faster (100x-1000x than ANCOMBC and according to the authors); 3)
it supports hierarchical models. The latest ANCOMBC versions are also
supporting hierarchical models. Nevertheless, LinDA seems a promising
tool that achieves a very good power/fdr trade-off together with
ANCOMBC according to the review. The speed improvements might make it
critical especially for datasets that have higher sample or feature
set sizes.



```r
meta <- as.data.frame(colData(tse)) %>% dplyr::select(Geographical_location)
linda.res <- linda(
  as.data.frame(assay(tse)), 
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
linda_out <- linda.res$output$Geographical_locationPune
```


```r
# to scan the table for genera where H0 could be rejected:
kable(head(filter(as.data.frame(linda_out), reject)))
```


\begin{tabular}{l|r|r|r|r|r|r|l|r}
\hline
  & baseMean & log2FoldChange & lfcSE & stat & pvalue & padj & reject & df\\
\hline
Acidaminococcus & 1194.7 & -1.9084 & 0.3579 & -5.332 & 0.0000 & 0.0000 & TRUE & 45\\
\hline
Aciditerrimonas & 393.5 & -0.6655 & 0.2184 & -3.048 & 0.0039 & 0.0165 & TRUE & 45\\
\hline
Actinomadura & 836.9 & -1.7985 & 0.3596 & -5.001 & 0.0000 & 0.0001 & TRUE & 45\\
\hline
Agromyces & 938.1 & -1.7530 & 0.3910 & -4.483 & 0.0001 & 0.0004 & TRUE & 45\\
\hline
Aminivibrio & 416.9 & -0.7489 & 0.2352 & -3.185 & 0.0026 & 0.0128 & TRUE & 45\\
\hline
Amycolatopsis & 556.1 & -0.9299 & 0.3320 & -2.801 & 0.0075 & 0.0302 & TRUE & 45\\
\hline
\end{tabular}

### ZicoSeq

Subsequently, we add a linear model and permutation-based method, see
details at [tutorial](https://cran.r-project.org/web/packages/GUniFrac/vignettes/ZicoSeq.html).

This approach has been assessed to exhibit high power and a low false
discovery rate, which has the following components:

  1. Winsorization to decrease the influence of outliers;
  
  1. Posterior sampling based on a beta mixture prior to address
  sampling variability and zero inflation;
  
  1. Reference-based multiple-stage normalization to address
  compositional effects;



```r
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
```

```
## 0  features are filtered!
## The data has  47  samples and  262  features will be tested!
## On average,  1  outlier counts will be replaced for each feature!
## Fitting beta mixture ...
## Finding the references ...
## Permutation testing ...
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## Completed!
```

```r
zicoseq_out <- cbind.data.frame(p.raw=zicoseq.obj$p.raw, p.adj.fdr=zicoseq.obj$p.adj.fdr) 
```

```r
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))
```


\begin{tabular}{l|r|r}
\hline
  & p.raw & p.adj.fdr\\
\hline
Actinomadura & 0.001 & 0.0001\\
\hline
Agromyces & 0.001 & 0.0038\\
\hline
Aneurinibacillus & 0.001 & 0.0071\\
\hline
Anoxybacillus & 0.001 & 0.0001\\
\hline
Chitinispirillum & 0.001 & 0.0100\\
\hline
Chitinivibrio & 0.001 & 0.0001\\
\hline
\end{tabular}



```r
## x-axis is the effect size: R2 * direction of coefficient
ZicoSeq.plot(ZicoSeq.obj = zicoseq.obj,
             meta.dat = meta,
	     pvalue.type ='p.adj.fdr')
```

![](30_differential_abundance_files/figure-latex/ZicoSeqplot-1.pdf)<!-- --> 



### Comparison of methods

The different methods yield somewhat different results but they could
be also expected to overlap to a substantial degree. As an exercirse,
you can compare the outcomes between the different methods in terms of
effect sizes, significances, or other aspects that are comparable
between the methods.




## Confounding variables

Confounders are common in experimental research. In general, these can be
classified into 3 types:

- Biological confounder, such as age, sex, etc. 

- Technical confounder that caused by data collection, storage, DNA
  extraction, sequencing process, etc.

- Confounder caused by experimental models, such as cage effect,
  sample background, etc.

Adjusting confounder is necessary and important to reach a valid
conclusion. To perform causal inference, it is crucial that the method
is able to include covariates in the model. This is not possible with
e.g. the Wilcoxon test. Other methods such as DESeq2, edgeR, ANCOMBC,
LDM, Aldex2, Corncob, MaAsLin2, ZicoSeq, fastANCOM and ZINQ allow
this. Below we show how to include a confounder/covariate in ANCOMBC,
LinDA and ZicoSeq.


### ANCOMBC


```r
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
```


```r
tab <- ancombc_cov$res %>% dplyr::select(starts_with(c("taxon", "lfc", "q")))
kable(head(tab))
```


\begin{tabular}{l|r|r|r|r|r|r|r|r}
\hline
taxon & lfc\_(Intercept) & lfc\_Geographical\_locationPune & lfc\_AgeElderly & lfc\_AgeMiddle\_age & q\_(Intercept) & q\_Geographical\_locationPune & q\_AgeElderly & q\_AgeMiddle\_age\\
\hline
Abyssicoccus & 0.0395 & -0.0947 & 0.0895 & 0.0127 & 0.8841 & 0.9271 & 0.9118 & 0.9948\\
\hline
Acidaminococcus & 0.7025 & -0.7557 & -0.2431 & -0.1576 & 0.0032 & 0.0388 & 0.8895 & 0.9774\\
\hline
Acinetobacter & 0.0202 & -1.0765 & 1.4995 & 1.1535 & 0.9880 & 0.6353 & 0.7718 & 0.7886\\
\hline
Actinomyces & 0.1762 & 0.1239 & -0.4581 & -0.4485 & 0.5453 & 0.8893 & 0.5720 & 0.6120\\
\hline
Actinoplanes & 0.3027 & -0.1600 & -0.2766 & -0.3340 & 0.2111 & 0.7623 & 0.8167 & 0.6120\\
\hline
Aerococcus & 0.0021 & -0.1441 & 0.1355 & 0.2462 & 0.9904 & 0.8209 & 0.8895 & 0.7886\\
\hline
\end{tabular}

### LinDA



```r
linda_cov <- linda(
  as.data.frame(assay(tse, "counts")), 
  as.data.frame(colData(tse)), 
  formula = '~ Geographical_location + Age', 
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
linda.res <- linda_cov$output$Geographical_locationPune
```


```r
kable(head(filter(linda.res, reject==T)))
```


\begin{tabular}{l|r|r|r|r|r|r|l|r}
\hline
  & baseMean & log2FoldChange & lfcSE & stat & pvalue & padj & reject & df\\
\hline
Acidaminococcus & 1140.5 & -1.438 & 0.4390 & -3.276 & 0.0021 & 0.0182 & TRUE & 43\\
\hline
Actinomadura & 804.7 & -1.520 & 0.4477 & -3.394 & 0.0015 & 0.0139 & TRUE & 43\\
\hline
Anaerococcus & 939.0 & 6.795 & 1.3978 & 4.861 & 0.0000 & 0.0004 & TRUE & 43\\
\hline
Aneurinibacillus & 988.0 & -1.687 & 0.5426 & -3.110 & 0.0033 & 0.0256 & TRUE & 43\\
\hline
Anoxybacillus & 1712.8 & -2.727 & 0.5574 & -4.893 & 0.0000 & 0.0004 & TRUE & 43\\
\hline
Brachybacterium & 403.0 & 2.060 & 0.7303 & 2.821 & 0.0072 & 0.0450 & TRUE & 43\\
\hline
\end{tabular}



### ZicoSeq


```r
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
```

```
## 0  features are filtered!
## The data has  47  samples and  262  features will be tested!
## On average,  1  outlier counts will be replaced for each feature!
## Fitting beta mixture ...
## Finding the references ...
## Permutation testing ...
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## ...................................................................................................
## Completed!
```

```r
zicoseq_out <- cbind.data.frame(p.raw=zicoseq.obj$p.raw,
                                p.adj.fdr=zicoseq.obj$p.adj.fdr) 
```


```r
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))
```


\begin{tabular}{l|r|r}
\hline
  & p.raw & p.adj.fdr\\
\hline
Actinomadura & 0.001 & 0.0004\\
\hline
Agromyces & 0.001 & 0.0030\\
\hline
Aneurinibacillus & 0.001 & 0.0036\\
\hline
Anoxybacillus & 0.001 & 0.0001\\
\hline
Chitinispirillum & 0.001 & 0.0089\\
\hline
Chitinivibrio & 0.001 & 0.0001\\
\hline
\end{tabular}



## Tree-based methods

Let us next cover phylogeny-aware methods to perform group-wise
associations.

### Group-wise associations testing based on balances

For testing associations based on balances, check the philr
R/Bioconductor package.


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
 [3] foreach_1.5.2                   GUniFrac_1.7                   
 [5] ANCOMBC_2.2.0                   MicrobiomeStat_1.1             
 [7] Maaslin2_1.7.3                  ALDEx2_1.32.0                  
 [9] zCompositions_1.4.0-1           truncnorm_1.0-9                
[11] NADA_1.6-1.1                    survival_3.5-5                 
[13] MASS_7.3-60                     lubridate_1.9.2                
[15] forcats_1.0.0                   stringr_1.5.0                  
[17] dplyr_1.1.2                     purrr_1.0.1                    
[19] readr_2.1.4                     tidyr_1.3.0                    
[21] tibble_3.2.1                    ggplot2_3.4.2                  
[23] tidyverse_2.0.0                 knitr_1.42                     
[25] tidySummarizedExperiment_1.10.0 patchwork_1.1.2                
[27] mia_1.9.2                       MultiAssayExperiment_1.26.0    
[29] TreeSummarizedExperiment_2.1.4  Biostrings_2.68.0              
[31] XVector_0.40.0                  SingleCellExperiment_1.22.0    
[33] SummarizedExperiment_1.30.1     Biobase_2.60.0                 
[35] GenomicRanges_1.52.0            GenomeInfoDb_1.36.0            
[37] IRanges_2.34.0                  S4Vectors_0.38.1               
[39] BiocGenerics_0.46.0             MatrixGenerics_1.12.0          
[41] matrixStats_0.63.0-9003         BiocStyle_2.28.0               
[43] rebook_1.9.0                   

loaded via a namespace (and not attached):
  [1] bitops_1.0-7                DirichletMultinomial_1.42.0
  [3] RColorBrewer_1.1-3          doParallel_1.0.17          
  [5] httr_1.4.5                  numDeriv_2016.8-1.1        
  [7] backports_1.4.1             tools_4.3.0                
  [9] utf8_1.2.3                  R6_2.5.1                   
 [11] vegan_2.6-4                 lazyeval_0.2.2             
 [13] mgcv_1.8-42                 rhdf5filters_1.12.1        
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
 [37] ggbeeswarm_0.7.2            fansi_1.0.4                
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
 [65] gtable_0.3.3                cachem_1.0.8               
 [67] xfun_0.39                   rbibutils_2.2.13           
 [69] S4Arrays_1.0.1              Rfast_2.0.7                
 [71] coda_0.19-4                 pcaPP_2.0-3                
 [73] modeest_2.4.0               timeDate_4022.108          
 [75] iterators_1.0.14            statmod_1.5.0              
 [77] gmp_0.7-1                   TH.data_1.1-2              
 [79] ellipsis_0.3.2              nlme_3.1-162               
 [81] phyloseq_1.44.0             bit64_4.0.5                
 [83] filelock_1.0.2              fBasics_4022.94            
 [85] irlba_2.3.5.1               vipor_0.4.5                
 [87] rpart_4.1.19                colorspace_2.1-0           
 [89] DBI_1.1.3                   Hmisc_5.0-1                
 [91] nnet_7.3-19                 ade4_1.7-22                
 [93] Exact_3.2                   tidyselect_1.2.0           
 [95] emmeans_1.8.5               timeSeries_4021.105        
 [97] bit_4.0.5                   compiler_4.3.0             
 [99] graph_1.78.0                htmlTable_2.4.1            
[101] BiocNeighbors_1.18.0        expm_0.999-7               
[103] DelayedArray_0.26.2         plotly_4.10.1              
[105] bookdown_0.33               checkmate_2.2.0            
[107] scales_1.2.1                DEoptimR_1.0-13            
[109] spatial_7.3-16              digest_0.6.31              
[111] minqa_1.2.5                 rmarkdown_2.21             
[113] base64enc_0.1-3             htmltools_0.5.5            
[115] pkgconfig_2.0.3             lme4_1.1-33                
[117] sparseMatrixStats_1.12.0    lpsymphony_1.28.0          
[119] highr_0.10                  stabledist_0.7-1           
[121] fastmap_1.1.1               rlang_1.1.1                
[123] htmlwidgets_1.6.2           DelayedMatrixStats_1.22.0  
[125] farver_2.1.1                energy_1.7-11              
[127] zoo_1.8-12                  jsonlite_1.8.4             
[129] BiocParallel_1.34.1         BiocSingular_1.16.0        
[131] RCurl_1.98-1.12             magrittr_2.0.3             
[133] Formula_1.2-5               scuttle_1.10.1             
[135] GenomeInfoDbData_1.2.10     Rhdf5lib_1.22.0            
[137] munsell_0.5.0               Rcpp_1.0.10                
[139] ape_5.7-1                   viridis_0.6.3              
[141] RcppZiggurat_0.1.6          CVXR_1.0-11                
[143] stringi_1.7.12              rootSolve_1.8.2.3          
[145] stable_1.1.6                zlibbioc_1.46.0            
[147] plyr_1.8.8                  parallel_4.3.0             
[149] ggrepel_0.9.3               lmom_2.9                   
[151] splines_4.3.0               hash_2.2.6.2               
[153] multtest_2.56.0             hms_1.1.3                  
[155] igraph_1.4.2                reshape2_1.4.4             
[157] ScaledMatrix_1.8.1          rmutil_1.1.10              
[159] XML_3.99-0.14               evaluate_0.20              
[161] BiocManager_1.30.20         nloptr_2.0.3               
[163] tzdb_0.3.0                  getopt_1.20.3              
[165] clue_0.3-64                 rsvd_1.0.5                 
[167] xtable_1.8-4                Rmpfr_0.9-2                
[169] e1071_1.7-13                tidytree_0.4.2             
[171] viridisLite_0.4.2           class_7.3-22               
[173] gsl_2.1-8                   lmerTest_3.1-3             
[175] memoise_2.0.1               beeswarm_0.4.0             
[177] cluster_2.1.4               timechange_0.2.0           
```
</div>


