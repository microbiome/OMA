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

Due to the complex data characteristics of microbiome sequencing data, differential abundance analysis of microbiome data faces many statistical challenges [@Yang2022], including:
  
- Highly variable. The abundance of a specific taxon could range over several orders of magnitude. 
  
- Zero-inflated. In a typical microbiome dataset, more than 70% of the values are zeros. Zeros could be due to either physical absence (structural zeros) or insufficient sampling effort (sampling zeros).
  
- Compositional. Increase or decrease in the (absolute) abundance of one taxon at the sampling site will lead to apparent changes in the relative abundances of other taxa in the sample. 

As summarized in @Yang2022, to address the above statistical chanllenegs:

- Over-dispersed count models has been proposed to address zero inflation, such as the negative binomial model used by edgeR [@Chen2016] and DESeq2 [@Love2014], the beta-binomial model used by corncorb [@Martin2021]. 

- Zero-inflated mixture models has aslo been proposed to address zero inflation, such as zero-inflated log-normal/normal mixture model used by metagenomeSeq [@Paulson2017] and RAIDA [@Sohn2015], zero-inflated beta-binomial model used by ZIBB [@ZIBB2018], and zero-inflated negative binomial model used by Omnibus [@Omnibus2018]. 

- Bayesian methods have been used to impute the zeros for methods working on proportion data, accounting for sampling variability and sequencing depth variation. Examples include ALDEx2 [@Gloor2016] and eBay [@Liu2020]. 

- Other methods use the pseudo-count approach to impute the zeros, such as MaAsLin2 [@Mallick2020] and ANCOMBC [@ancombc2020].

- Different strategies have been used to address compositional effects, including:

  - Robust normalization. For example, trimmed mean of M-values (TMM) normalization used by edgeR, relative log expression (RLE) normalization used by DESeq2 [@Love2014], cumulative sum scaling (CSS) normalization used by metagenomeSeq, centered log-ratio transformation (CLR) normalization used by ALDEx2 [@Gloor2016] and geometric mean of pairwise ratios (GMPR) normalization used by Omnibus [@Omnibus2018]. Wrench normalization [@Kumar2018] corrects the compositional bias by an empirical Bayes approach, which has been recommended in metagenomeSeq [@Paulson2017] tutorial.
  - Reference taxa approach used by DACOMP [@Brill2019] and RAIDA [@Sohn2015]. 
  - Analyzing the pattern of pairwise log ratios, such as ANCOM [@Mandal2015].
  - Bias-correction used by ANCOMBC [@ancombc2020].

The most popular tools, without going into evaluating whether or not they perform well for this task, are:  

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


We recommend to have a look at @Nearing2022
who compared all these listed methods across 38
different datasets. Because different methods use different approaches
(parametric vs non-parametric, different normalization techiniques, assumptions
etc.), results can differ between methods. 
Unfortunately, as @Nearing2022 point out, they can differ substantially. Also, more recently @Yang2022 comprehensively evaluated these methods via a Semi-parametric framework and 106 real datasets. @Yang2022 also pointed out, different DAA tools could sometimes produce quite discordant results, opening to the possibility of cherry-picking the tool in favor of one’s own hypothesis. Therefore, it is highly recommended to pick several methods to get an idea about how robust and potentially reproducible your findings are depending on the method. In this section we demonstrate 4 methods that can be recommended based on recent literature (ANCOM-BC [@ancombc2020], ALDEx2 [@Gloor2016], Maaslin2 [@Mallick2020], LinDA [@Zhou2022] and ZicoSeq [@Yang2022]) and we will compare the results between them.

Note that the purpose of this section is to show how to perform DAA in R, not
how to correctly do causal inference. Depending on your experimental setup
and your theory, you must determine how to specify any model exactly. 
E.g., there might be confounding factors that might drive (the absence of)
differences between the shown groups that we ignore here for simplicity. Or your dataset is repeated sampling design, matched-pair design or the general longitudianl design.
However, we will show how you could include covariates in those models.
Furthermore, we picked a dataset that merely has microbial abundances in a TSE
object as well as a grouping variable in the sample data. We simplify the
analysis by only including 2 of the 3 groups. 



```r
library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(knitr)
library(tidyverse)
library(phyloseq)
library(ALDEx2)
library(Maaslin2)
library(MicrobiomeStat)

if(!require(ANCOMBC)){
  BiocManager::install("ANCOMBC")
  library(ANCOMBC)

}
if(!require(GUniFrac)){
  install.packages("GUniFrac")
  library(GUniFrac)
}


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
  dplyr::select(genus, we.eBH, wi.eBH, effect, overlap) %>%
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
# make phyloseq object
otu <- otu_table(assay(tse), taxa_are_rows = T)
taxa <- tax_table(as.matrix(rowData(tse)))
meta <- sample_data(as.data.frame(colData(tse)))
phy <- phyloseq(otu,taxa,meta)
# perform the analysis 
out <- ancombc(
  phyloseq = phy,
  formula = "Geographical_location", 
  p_adj_method = "fdr", 
  lib_cut = 0, 
  group = "Geographical_location", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, 
  global = TRUE # multi group comparison will be deactivated automatically 
)
# store the FDR adjusted results [test on v2.0.3] 
ancombc_out <- cbind.data.frame(taxid = out$res$q_val$taxon, ancombc = as.vector(out$res$q_val$Geographical_locationPune))
# store the FDR adjusted results [test on v1.2.2] 
# ancombc_out <- out$res$q_val %>% rownames_to_column('taxid') %>% dplyr::rename(ancombc = 2)
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
kable(head(ancombc_out))
```


\begin{tabular}{l|r}
\hline
taxid & ancombc\\
\hline
OTU2 & 0.5501\\
\hline
OTU15 & 0.0253\\
\hline
OTU22 & 0.0000\\
\hline
OTU53 & 0.8949\\
\hline
OTU69 & 0.0000\\
\hline
OTU76 & 0.0000\\
\hline
\end{tabular}



### MaAsLin2 

Next, we will illustrate the use of MaAsLin2 [@Mallick2020]. The method is based on
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
maaslin2_out <- Maaslin2(
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
kable(head(filter(maaslin2_out$results, qval <= 0.05)))
```


\begin{tabular}{l|l|l|r|r|r|l|r|r|r}
\hline
feature & metadata & value & coef & stderr & pval & name & qval & N & N.not.zero\\
\hline
OTU1053 & Geographical\_location & Nashik & 0.0080 & 0.0011 & 0 & Geographical\_locationNashik & 0 & 47 & 9\\
\hline
OTU860 & Geographical\_location & Nashik & 0.0373 & 0.0059 & 0 & Geographical\_locationNashik & 0 & 47 & 13\\
\hline
OTU1075 & Geographical\_location & Nashik & 0.1295 & 0.0207 & 0 & Geographical\_locationNashik & 0 & 47 & 27\\
\hline
OTU1980 & Geographical\_location & Nashik & 0.0395 & 0.0062 & 0 & Geographical\_locationNashik & 0 & 47 & 9\\
\hline
OTU611 & Geographical\_location & Nashik & 0.0274 & 0.0045 & 0 & Geographical\_locationNashik & 0 & 47 & 10\\
\hline
OTU2335 & Geographical\_location & Nashik & 0.0089 & 0.0015 & 0 & Geographical\_locationNashik & 0 & 47 & 10\\
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
meta <- as.data.frame(colData(tse)) %>% dplyr::select(Geographical_location)
linda.res <- linda(
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
linda_out <- linda.res$output$Geographical_locationPune
# to scan the table for genera where H0 could be rejected:
kable(head(filter(as.data.frame(linda_out), reject)))
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

### ZicoSeq

Subsequently, we add a linear model and permutation-based method, see details at [tutorial](https://cran.r-project.org/web/packages/GUniFrac/vignettes/ZicoSeq.html). This approach has been assessed to exhibit high power and a low false discovery rate, which has the following components: 

  - 1) Winsorization to decrease the influence of outliers;
  - 2) Posterior sampling based on a beta mixture prior to address sampling variability and zero inflation;
  - 3) Reference-based multiple-stage normalization to address compositional effects;



```r
set.seed(123)
otu.tab <- as.matrix(assay(tse))
meta <- as.data.frame(colData(tse)) 
zicoseq.obj <- GUniFrac::ZicoSeq(meta.dat = meta, 
                                 feature.dat = otu.tab,
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
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))
```


\begin{tabular}{l|r|r}
\hline
  & p.raw & p.adj.fdr\\
\hline
OTU76 & 0.001 & 0.0001\\
\hline
OTU127 & 0.001 & 0.0038\\
\hline
OTU207 & 0.001 & 0.0094\\
\hline
OTU211 & 0.001 & 0.0001\\
\hline
OTU610 & 0.001 & 0.0094\\
\hline
OTU611 & 0.001 & 0.0001\\
\hline
\end{tabular}



```r
## x-axis is the effect size: R2 * direction of coefficient
ZicoSeq.plot(ZicoSeq.obj = zicoseq.obj, meta.dat = meta, pvalue.type ='p.adj.fdr')
```

![](30_differential_abundance_files/figure-latex/ZicoSeqplot-1.pdf)<!-- --> 



### Comparison of the methods

When we compare the methods in the context of a research question, we could
look at e.g. at whether they agree based on the applied decision criterion
(e.g. adjusted p value <= 0.05). That is what we illustrate here. First we will 
look at how many taxa were identified by each method to begin with. In the next
step we will look at the intersection of identified taxa. To achieve that, we
first create a dataframe that summarises the decision criterion for each method
and shows a score from 0 to 3 indicating how many methods agreed on a particular
taxon.



```r
summ <- full_join(rownames_to_column(aldex_out, "taxid") %>% dplyr::select(taxid, aldex2 = wi.eBH),
    ancombc_out,by = "taxid") %>%
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

# This is how it looks like:
kable(head(summ))
```


\begin{tabular}{l|l|l|l|l|l|r}
\hline
taxid & aldex2 & ancombc & maaslin2 & LinDA & ZicoSeq & score\\
\hline
OTU2 & FALSE & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
OTU15 & FALSE & TRUE & TRUE & TRUE & FALSE & 3\\
\hline
OTU22 & FALSE & TRUE & TRUE & TRUE & FALSE & 3\\
\hline
OTU53 & FALSE & FALSE & FALSE & FALSE & FALSE & 0\\
\hline
OTU69 & FALSE & TRUE & FALSE & FALSE & FALSE & 1\\
\hline
OTU76 & FALSE & TRUE & TRUE & TRUE & TRUE & 4\\
\hline
\end{tabular}

Now we can answer our questions:


```r
# how many genera were identified by each method?
summarise(summ, across(where(is.logical), sum)) %>%
  kable()
```


\begin{tabular}{r|r|r|r|r}
\hline
aldex2 & ancombc & maaslin2 & LinDA & ZicoSeq\\
\hline
9 & 165 & 67 & 75 & 18\\
\hline
\end{tabular}

```r
# which genera are identified by all methods?
filter(summ, score == 5) %>% kable()
```


\begin{tabular}{l|l|l|l|l|l|r}
\hline
taxid & aldex2 & ancombc & maaslin2 & LinDA & ZicoSeq & score\\
\hline
OTU611 & TRUE & TRUE & TRUE & TRUE & TRUE & 5\\
\hline
OTU773 & TRUE & TRUE & TRUE & TRUE & TRUE & 5\\
\hline
OTU860 & TRUE & TRUE & TRUE & TRUE & TRUE & 5\\
\hline
OTU1075 & TRUE & TRUE & TRUE & TRUE & TRUE & 5\\
\hline
OTU2529 & TRUE & TRUE & TRUE & TRUE & TRUE & 5\\
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
tse <- transformCounts(tse, method="clr", pseudocount=1) # not bale to run

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
  ancombc_plots[[6]] +
  plot_layout(nrow = 1)
```

![](30_differential_abundance_files/figure-latex/daplotting-2.pdf)<!-- --> 



### Confounding variables

Confounders are common for microbiome studies. In general, it can be classified into 3 types:

- Biological confounder, such as age, sex, etc. 
- Technical confounder that caused by data collection, storage, DNA extraction, sequencing process, etc. 
- Confounder caused by experimental models, such as cage effect, sample background, etc. 

Adjusting confounder is necessary and important to reach a valid conclusion. To perform causal inference, it is crucial that the method is able to include covariates in the model. This is not possible with e.g. the Wilcoxon test. Other methods such as DESeq2, edgeR, ANCOMBC, LDM, Aldex2, Corncob, MaAsLin2, ZicoSeq, fastANCOM and ZINQ allow this. Below we show how to include a confounder/covariate in ANCOMBC, LinDA and ZicoSeq.


#### ANCOMBC


```r
# make phyloseq object
otu <- otu_table(assay(tse), taxa_are_rows = T)
taxa <- tax_table(as.matrix(rowData(tse)))
meta <- sample_data(as.data.frame(colData(tse)))
phy <- phyloseq(otu,taxa,meta)
# perform the analysis 
ancombc_cov <- ancombc(
  phyloseq = phy,
  formula = "Geographical_location + Age", 
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
```


\begin{tabular}{l|r|r|r|r}
\hline
taxon & (Intercept) & Geographical\_locationPune & AgeElderly & AgeMiddle\_age\\
\hline
OTU2 & 0.3110 & 0.9049 & 0.8293 & 0.8582\\
\hline
OTU15 & 0.0122 & 0.0188 & 0.4545 & 0.9543\\
\hline
OTU22 & 0.9749 & 0.0000 & 0.7095 & 0.3344\\
\hline
OTU53 & 0.2281 & 0.5025 & 0.4272 & 0.4503\\
\hline
OTU69 & 0.8700 & 0.0000 & 0.3927 & 0.3892\\
\hline
OTU76 & 0.1172 & 0.0000 & 0.7337 & 0.8582\\
\hline
\end{tabular}

#### LinDA


```r
otu.tab <- as.data.frame(assay(tse))
meta <- as.data.frame(colData(tse))
linda_cov <- linda(
  otu.tab, 
  meta, 
  formula = '~ Geographical_location + Age', 
  alpha = 0.05, 
  prev.filter = 0, 
  mean.abund.filter = 0)
```

```
## 0  features are filtered!
## The filtered data has  47  samples and  275  features will be tested!
## Pseudo-count approach is used.
## Fit linear models ...
## Completed.
```

```r
linda.res <- linda_cov$output$Geographical_locationPune
kable(head(filter(linda.res, reject==T)))
```


\begin{tabular}{l|r|r|r|r|r|r|l|r}
\hline
  & baseMean & log2FoldChange & lfcSE & stat & pvalue & padj & reject & df\\
\hline
OTU15 & 1137.8 & -1.4861 & 0.4407 & -3.372 & 0.0016 & 0.0146 & TRUE & 43\\
\hline
OTU22 & 363.8 & -0.7457 & 0.2603 & -2.865 & 0.0064 & 0.0442 & TRUE & 43\\
\hline
OTU76 & 802.8 & -1.5674 & 0.4509 & -3.477 & 0.0012 & 0.0116 & TRUE & 43\\
\hline
OTU127 & 910.8 & -1.3446 & 0.4880 & -2.756 & 0.0086 & 0.0498 & TRUE & 43\\
\hline
OTU194 & 936.8 & 6.7471 & 1.3990 & 4.823 & 0.0000 & 0.0005 & TRUE & 43\\
\hline
OTU207 & 985.6 & -1.7352 & 0.5444 & -3.187 & 0.0027 & 0.0223 & TRUE & 43\\
\hline
\end{tabular}



#### ZicoSeq


```r
set.seed(123)
otu.tab <- as.matrix(assay(tse))
meta <- as.data.frame(colData(tse)) 
zicoseq.obj <- GUniFrac::ZicoSeq(meta.dat = meta, 
                                 feature.dat = otu.tab,
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
## The data has  47  samples and  275  features will be tested!
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
kable(head(filter(zicoseq_out, p.adj.fdr<0.05)))
```


\begin{tabular}{l|r|r}
\hline
  & p.raw & p.adj.fdr\\
\hline
OTU76 & 0.001 & 0.0008\\
\hline
OTU127 & 0.001 & 0.0052\\
\hline
OTU207 & 0.001 & 0.0064\\
\hline
OTU211 & 0.001 & 0.0001\\
\hline
OTU610 & 0.001 & 0.0194\\
\hline
OTU611 & 0.001 & 0.0001\\
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
 [3] foreach_1.5.2                   GUniFrac_1.7                   
 [5] ANCOMBC_2.2.0                   MicrobiomeStat_1.1             
 [7] Maaslin2_1.7.3                  ALDEx2_1.32.0                  
 [9] zCompositions_1.4.0-1           truncnorm_1.0-9                
[11] NADA_1.6-1.1                    survival_3.5-5                 
[13] MASS_7.3-60                     phyloseq_1.44.0                
[15] lubridate_1.9.2                 forcats_1.0.0                  
[17] stringr_1.5.0                   dplyr_1.1.2                    
[19] purrr_1.0.1                     readr_2.1.4                    
[21] tidyr_1.3.0                     tibble_3.2.1                   
[23] ggplot2_3.4.2                   tidyverse_2.0.0                
[25] knitr_1.42                      tidySummarizedExperiment_1.10.0
[27] patchwork_1.1.2                 mia_1.9.2                      
[29] MultiAssayExperiment_1.26.0     TreeSummarizedExperiment_2.1.4 
[31] Biostrings_2.68.0               XVector_0.40.0                 
[33] SingleCellExperiment_1.22.0     SummarizedExperiment_1.30.1    
[35] Biobase_2.60.0                  GenomicRanges_1.52.0           
[37] GenomeInfoDb_1.36.0             IRanges_2.34.0                 
[39] S4Vectors_0.38.1                BiocGenerics_0.46.0            
[41] MatrixGenerics_1.12.0           matrixStats_0.63.0-9003        
[43] BiocStyle_2.28.0                rebook_1.9.0                   

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
 [81] bit64_4.0.5                 filelock_1.0.2             
 [83] fBasics_4022.94             irlba_2.3.5.1              
 [85] vipor_0.4.5                 rpart_4.1.19               
 [87] colorspace_2.1-0            DBI_1.1.3                  
 [89] Hmisc_5.0-1                 nnet_7.3-19                
 [91] ade4_1.7-22                 Exact_3.2                  
 [93] tidyselect_1.2.0            emmeans_1.8.5              
 [95] timeSeries_4021.105         bit_4.0.5                  
 [97] compiler_4.3.0              graph_1.78.0               
 [99] htmlTable_2.4.1             BiocNeighbors_1.18.0       
[101] expm_0.999-7                DelayedArray_0.26.1        
[103] plotly_4.10.1               bookdown_0.33              
[105] checkmate_2.2.0             scales_1.2.1               
[107] DEoptimR_1.0-13             spatial_7.3-16             
[109] digest_0.6.31               minqa_1.2.5                
[111] rmarkdown_2.21              htmltools_0.5.5            
[113] pkgconfig_2.0.3             base64enc_0.1-3            
[115] lme4_1.1-33                 sparseMatrixStats_1.12.0   
[117] lpsymphony_1.28.0           highr_0.10                 
[119] stabledist_0.7-1            fastmap_1.1.1              
[121] rlang_1.1.1                 htmlwidgets_1.6.2          
[123] DelayedMatrixStats_1.22.0   farver_2.1.1               
[125] energy_1.7-11               zoo_1.8-12                 
[127] jsonlite_1.8.4              BiocParallel_1.34.0        
[129] BiocSingular_1.16.0         RCurl_1.98-1.12            
[131] magrittr_2.0.3              Formula_1.2-5              
[133] scuttle_1.10.1              GenomeInfoDbData_1.2.10    
[135] Rhdf5lib_1.22.0             munsell_0.5.0              
[137] Rcpp_1.0.10                 ape_5.7-1                  
[139] viridis_0.6.3               RcppZiggurat_0.1.6         
[141] CVXR_1.0-11                 stringi_1.7.12             
[143] rootSolve_1.8.2.3           stable_1.1.6               
[145] zlibbioc_1.46.0             plyr_1.8.8                 
[147] parallel_4.3.0              ggrepel_0.9.3              
[149] lmom_2.9                    splines_4.3.0              
[151] hash_2.2.6.2                multtest_2.56.0            
[153] hms_1.1.3                   igraph_1.4.2               
[155] reshape2_1.4.4              ScaledMatrix_1.8.1         
[157] rmutil_1.1.10               XML_3.99-0.14              
[159] evaluate_0.20               BiocManager_1.30.20        
[161] nloptr_2.0.3                tzdb_0.3.0                 
[163] getopt_1.20.3               clue_0.3-64                
[165] rsvd_1.0.5                  xtable_1.8-4               
[167] Rmpfr_0.9-2                 e1071_1.7-13               
[169] tidytree_0.4.2              viridisLite_0.4.2          
[171] class_7.3-22                gsl_2.1-8                  
[173] lmerTest_3.1-3              memoise_2.0.1              
[175] beeswarm_0.4.0              cluster_2.1.4              
[177] timechange_0.2.0           
```
</div>


