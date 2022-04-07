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
identifying that a bacterial taxon is different between e.g. a patient
group with disease *X* vs a healthy control group might lead to
important insights into the pathophysiology. Changes in the microbiota
might be causal or a consequence of the disease. Either way, it can
help to understand the system as a whole. Be aware that this approach
has also been criticized recently [@quinnCritiqueDifferentialAbundance2021].


### Examples and tools

There are many tools to perform DAA. The most popular tools, without going into
evaluating whether or not they perform well for this task, are:  
- [ALDEx2](https://bioconductor.org/packages/release/bioc/html/ALDEx2.html)   
- [ANCOM-BC](https://bioconductor.org/packages/release/bioc/html/ANCOMBC.html)  
- [corncob](https://cran.r-project.org/web/packages/corncob/index.html)  
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)  
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)  
- [LEFse](https://bioconductor.org/packages/release/bioc/html/lefser.html)  
- [MaAsLin2](https://www.bioconductor.org/packages/release/bioc/html/Maaslin2.html)  
- [metagenomeSeq](https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html)  
- [limma voom](https://bioconductor.org/packages/release/bioc/html/limma.html)  
- [t-test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/t.test)  
- [Wilcoxon test](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test)  


We recommend to have a look at @nearingMicrobiomeDifferentialAbundance2021 who compared all these listed methods across 38
different datasets. Because different methods have different approaches
(parametric vs non-parametric, different normalization techiniques etc.) to
perform the same task (test differential abundance), results can differ between
methods. Unfortunately, as @nearingMicrobiomeDifferentialAbundance2021 point out, they differ disturbingly
much. Therefore, it is highly recommended to pick several methods to get an
idea about how robust and potentially reproducible your findings are depending
on the method. In this section we demonstrate 3 methods that can be recommended
based on this recent review (ANCOM-BC, ALDEx2 and Maaslin2) and we will compare
the results between them.
Note that the purpose of this section is to show how to perform DAA in R, not
how to correctly do causal inference. E.g. there might be confounding factors
that might drive (the absence of) differences between the shown groups that we
ignore for simplicity. However, we will show how you could include covariates
in those models. Furthermore, we picked a dataset that merely has
microbial abundances in a TSE object as well as a grouping variable in the
sample data. We simplify the analysis by only including 2 of the 3 groups. 





```r
library(mia)
library(patchwork)
library(tidySummarizedExperiment)
library(ANCOMBC)
library(ALDEx2)
library(Maaslin2)
library(knitr)
library(tidyverse)

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
@nearingMicrobiomeDifferentialAbundance2021 found that applying a 10% threshold for the prevalence of
the taxa generally resulted in more robust results. Some tools have builtin
arguments for that. By applying the threshold to our input data, we can make
sure it is applied for all tools. Below we show how to do this in `mia`:


```r
tse <- subsetByPrevalentTaxa(tse, detection = 0, prevalence = 0.1)
```

### ALDEx2

In this section, we will show how to perform a simple ALDEx2 analysis. If you
would choose to pick a single method, this method could be recommended to use.
According to the developers experience, it tends to identify the common
features identified by other methods. This statement is in line with a recent
independent evaluation by @nearingMicrobiomeDifferentialAbundance2021.  Please also have a look at
the more extensive [vignette](https://bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html) 
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
rare taxa and grey ones are abundant taxa.


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
the authors recommend an effect 
size cutoff of 1 rather than only interpreting the pvalue. Again, this is not 
the case for any feature.



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
(ANCOM-BC) [@linAnalysisCompositionsMicrobiomes2020] 
is a recently developed method for differential abundance testing. It is based
on an 
earlier published approach [@mandalAnalysisCompositionMicrobiomes2015]. 
The previous version of ANCOM was among the methods that produced the 
most consistent results and is probably a conservative approach
[@nearingMicrobiomeDifferentialAbundance2021]. 
However, the new ANCOM-BC method operates quite differently compared to the
former ANCOM method.


As the only method, ANCOM-BC incorporates the so called *sampling fraction*
into the model. The latter term could be empirically estimated by the ratio of
the library size to the microbial load. According to the authors, variations in
this sampling fraction would bias differential abundance analyses if ignored.
Furthermore, this method provides p-values, and confidence intervals for each
taxon. It also controls the FDR and it is computationally simple to implement. 


As we will see below, to obtain results, all that is needed is to pass 
a phyloseq object to the `ancombc()` function. Therefore, below we first
convert our `TreeSE` object to a `phyloseq` object. Then, we specify the formula.
In this formula, other covariates could potentially be included to adjust for
confounding. We show this further below. 
Please check the [function documentation](https://rdrr.io/github/FrederickHuangLin/ANCOMBC/man/ancombc.html) 
to learn about the additional arguments that we specify below.


```r
# currently, ancombc requires the phyloseq format, but we can easily convert:
pseq <- makePhyloseqFromTreeSummarizedExperiment(tse)

# perform the analysis 
out = ancombc(
  phyloseq = pseq, 
  formula = "pheno", 
  p_adj_method = "fdr", 
  zero_cut = 1, # no prev filtering necessary anymore 
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
abundant between the groups. Below we show the first 6 entries of this
dataframe:  


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

Lastly, we will illustrate how to use MaAsLin2, which is the next generation of
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
  mutate(
    across(c(aldex2, maaslin2), ~ .x <= 0.05),
    # the following line would be necessary without prevalence filtering 
    # as some methods output NA
    #across(-genus, function(x) ifelse(is.na(x), FALSE, x)),
    score = rowSums(across(c(aldex2, ancombc, maaslin2)))
  )

# This is how it looks like:
kable(head(summ))
```


\begin{tabular}{l|l|l|l|r}
\hline
genus & aldex2 & ancombc & maaslin2 & score\\
\hline
Acetanaerobacterium & FALSE & TRUE & TRUE & 2\\
\hline
Acetivibrio & FALSE & FALSE & FALSE & 0\\
\hline
Acidaminococcus & FALSE & TRUE & TRUE & 2\\
\hline
Akkermansia & FALSE & FALSE & FALSE & 0\\
\hline
Alistipes & TRUE & TRUE & TRUE & 3\\
\hline
Allisonella & FALSE & FALSE & FALSE & 0\\
\hline
\end{tabular}

Now we can answer our questions:


```r
# how many genera were identified by each method?
summarise(summ, across(where(is.logical), sum)) %>%
  kable()
```


\begin{tabular}{r|r|r}
\hline
aldex2 & ancombc & maaslin2\\
\hline
9 & 22 & 16\\
\hline
\end{tabular}

```r
# which genera are identified by all methods?
filter(summ, score == 3) %>% kable()
```


\begin{tabular}{l|l|l|l|r}
\hline
genus & aldex2 & ancombc & maaslin2 & score\\
\hline
Alistipes & TRUE & TRUE & TRUE & 3\\
\hline
Barnesiella & TRUE & TRUE & TRUE & 3\\
\hline
Catenibacterium & TRUE & TRUE & TRUE & 3\\
\hline
Lactobacillus & TRUE & TRUE & TRUE & 3\\
\hline
Megasphaera & TRUE & TRUE & TRUE & 3\\
\hline
Oscillibacter & TRUE & TRUE & TRUE & 3\\
\hline
Parabacteroides & TRUE & TRUE & TRUE & 3\\
\hline
Phascolarctobacterium & TRUE & TRUE & TRUE & 3\\
\hline
\end{tabular}

We see that each method identified at least 9 genera as differentially
abundant. Eight of those that were identified by ALDEx2,
were also identified by both of the other methods. We could plot the data for 
any method or for those taxa that were identified by all methods:




```r
plot_data <- data.frame(t(assay(tse)))
plot_data$pheno <- colData(tse)$pheno
# create a plot for each genus where the score is indicated in the title
plots <- pmap(select(summ, genus, score), function(genus, score) {
  ggplot(plot_data, aes_string("pheno", genus)) +
    geom_boxplot(aes(fill = pheno), outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    ggtitle(glue::glue("Robustness score {score}")) +
    theme_bw() +
    theme(legend.position = "none")
})

# now we can show only those genera that have at least score 3 (or 2 or 1)
robust_plots <- plots[summ$score == 3] 


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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> 

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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-11-2.pdf)<!-- --> 



### Confounding variables

To perform causal inference, it is crucial that the method is able to include
covariates in the model. This is not possible with e.g. the Wilcoxon test.
Other methods such as both ANCOM methods, ALDEx2, DESeq2, MaAsLin2 and others
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
pseq <- makePhyloseqFromTreeSummarizedExperiment(tse)
out_cov = ancombc(
  phyloseq = pseq, 
  formula = "pheno + age", # here we add age to the model
  p_adj_method = "fdr", 
  zero_cut = 0.90, 
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
data of @sprockettMicrobiotaAssemblyStructure2020 available
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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-20-1.pdf)<!-- --> 

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

![](30_differential_abundance_files/figure-latex/unnamed-chunk-23-1.pdf)<!-- --> 

Seemingly the covariate "Age_Years" does not have effect on the model as "Delivery_Mode" would,
and "Sex" to some extent. Let's take a closer look at the two latter ones:


```r
plot(posterior, par="Lambda", focus.cov = rownames(X)[c(2,4)])
```

![](30_differential_abundance_files/figure-latex/unnamed-chunk-24-1.pdf)<!-- --> 

## Session Info {-}

<button class="rebook-collapse">View session info</button>
<div class="rebook-content">
```
R version 4.1.3 (2022-03-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.4 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

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
 [1] microbiomeDataSets_1.1.5       fido_1.0.0                    
 [3] forcats_0.5.1                  stringr_1.4.0                 
 [5] dplyr_1.0.8                    purrr_0.3.4                   
 [7] readr_2.1.2                    tidyr_1.2.0                   
 [9] tibble_3.1.6                   ggplot2_3.3.5                 
[11] tidyverse_1.3.1                knitr_1.38                    
[13] Maaslin2_1.8.0                 ALDEx2_1.26.0                 
[15] zCompositions_1.4.0-1          truncnorm_1.0-8               
[17] NADA_1.6-1.1                   survival_3.3-1                
[19] MASS_7.3-56                    ANCOMBC_1.4.0                 
[21] tidySummarizedExperiment_1.4.1 patchwork_1.1.1               
[23] mia_1.3.19                     MultiAssayExperiment_1.20.0   
[25] TreeSummarizedExperiment_2.1.4 Biostrings_2.62.0             
[27] XVector_0.34.0                 SingleCellExperiment_1.16.0   
[29] SummarizedExperiment_1.24.0    Biobase_2.54.0                
[31] GenomicRanges_1.46.1           GenomeInfoDb_1.30.1           
[33] IRanges_2.28.0                 S4Vectors_0.32.4              
[35] BiocGenerics_0.40.0            MatrixGenerics_1.6.0          
[37] matrixStats_0.61.0-9003        BiocStyle_2.22.0              
[39] rebook_1.4.0                  

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3                coda_0.19-4                  
  [3] bit64_4.0.5                   irlba_2.3.5                  
  [5] DelayedArray_0.20.0           data.table_1.14.2            
  [7] KEGGREST_1.34.0               RCurl_1.98-1.6               
  [9] generics_0.1.2                ScaledMatrix_1.2.0           
 [11] microbiome_1.16.0             RSQLite_2.2.12               
 [13] bit_4.0.4                     tzdb_0.3.0                   
 [15] httpuv_1.6.5                  xml2_1.3.3                   
 [17] lubridate_1.8.0               assertthat_0.2.1             
 [19] DirichletMultinomial_1.36.0   viridis_0.6.2                
 [21] xfun_0.30                     hms_1.1.1                    
 [23] ggdist_3.1.1                  promises_1.2.0.1             
 [25] evaluate_0.15                 DEoptimR_1.0-11              
 [27] fansi_1.0.3                   dbplyr_2.1.1                 
 [29] readxl_1.4.0                  igraph_1.3.0                 
 [31] DBI_1.1.2                     htmlwidgets_1.5.4            
 [33] tensorA_0.36.2                hash_2.2.6.2                 
 [35] ellipsis_0.3.2                backports_1.4.1              
 [37] bookdown_0.25                 permute_0.9-7                
 [39] sparseMatrixStats_1.6.0       vctrs_0.4.0                  
 [41] abind_1.4-5                   tidybayes_3.0.2              
 [43] cachem_1.0.6                  withr_2.5.0                  
 [45] robustbase_0.95-0             checkmate_2.0.0              
 [47] vegan_2.5-7                   treeio_1.18.1                
 [49] getopt_1.20.3                 cluster_2.1.3                
 [51] ExperimentHub_2.2.1           ape_5.6-2                    
 [53] dir.expiry_1.2.0              lazyeval_0.2.2               
 [55] crayon_1.5.1                  pkgconfig_2.0.3              
 [57] labeling_0.4.2                nlme_3.1-157                 
 [59] vipor_0.4.5                   rlang_1.0.2                  
 [61] lifecycle_1.0.1               filelock_1.0.2               
 [63] BiocFileCache_2.2.1           phyloseq_1.38.0              
 [65] modelr_0.1.8                  rsvd_1.0.5                   
 [67] AnnotationHub_3.2.2           cellranger_1.1.0             
 [69] distributional_0.3.0          graph_1.72.0                 
 [71] Matrix_1.4-1                  lpsymphony_1.22.0            
 [73] Rhdf5lib_1.16.0               reprex_2.0.1                 
 [75] beeswarm_0.4.0                png_0.1-7                    
 [77] viridisLite_0.4.0             bitops_1.0-7                 
 [79] rhdf5filters_1.6.0            blob_1.2.2                   
 [81] DelayedMatrixStats_1.16.0     decontam_1.14.0              
 [83] DECIPHER_2.22.0               beachmat_2.10.0              
 [85] scales_1.1.1                  memoise_2.0.1                
 [87] magrittr_2.0.3                plyr_1.8.7                   
 [89] zlibbioc_1.40.0               compiler_4.1.3               
 [91] RColorBrewer_1.1-3            cli_3.2.0                    
 [93] ade4_1.7-18                   pbapply_1.5-0                
 [95] mgcv_1.8-40                   tidyselect_1.1.2             
 [97] stringi_1.7.6                 highr_0.9                    
 [99] yaml_2.3.5                    BiocSingular_1.10.0          
[101] svUnit_1.0.6                  ggrepel_0.9.1                
[103] grid_4.1.3                    tools_4.1.3                  
[105] parallel_4.1.3                rstudioapi_0.13              
[107] foreach_1.5.2                 logging_0.10-108             
[109] optparse_1.7.1                gridExtra_2.3                
[111] posterior_1.2.1               farver_2.1.0                 
[113] Rtsne_0.15                    RcppZiggurat_0.1.6           
[115] digest_0.6.29                 BiocManager_1.30.16          
[117] shiny_1.7.1                   Rcpp_1.0.8.3                 
[119] broom_0.7.12                  scuttle_1.4.0                
[121] later_1.3.0                   BiocVersion_3.14.0           
[123] AnnotationDbi_1.56.2          httr_1.4.2                   
[125] Rdpack_2.3                    colorspace_2.0-3             
[127] rvest_1.0.2                   XML_3.99-0.9                 
[129] fs_1.5.2                      splines_4.1.3                
[131] yulab.utils_0.0.4             tidytree_0.3.9               
[133] scater_1.22.0                 multtest_2.50.0              
[135] plotly_4.10.0                 xtable_1.8-4                 
[137] jsonlite_1.8.0                nloptr_2.0.0                 
[139] CodeDepends_0.6.5             Rfast_2.0.6                  
[141] R6_2.5.1                      mime_0.12                    
[143] pillar_1.7.0                  htmltools_0.5.2              
[145] glue_1.6.2                    fastmap_1.1.0                
[147] BiocParallel_1.28.3           BiocNeighbors_1.12.0         
[149] interactiveDisplayBase_1.32.0 codetools_0.2-18             
[151] pcaPP_1.9-74                  mvtnorm_1.1-3                
[153] utf8_1.2.2                    lattice_0.20-45              
[155] arrayhelpers_1.1-0            curl_4.3.2                   
[157] ggbeeswarm_0.6.0              biglm_0.9-2.1                
[159] rmarkdown_2.13                biomformat_1.22.0            
[161] munsell_0.5.0                 rhdf5_2.38.1                 
[163] GenomeInfoDbData_1.2.7        iterators_1.0.14             
[165] haven_2.4.3                   reshape2_1.4.4               
[167] gtable_0.3.0                  rbibutils_2.2.7              
```
</div>
