--- 
title: "Orchestrating Microbiome Analysis"
documentclass: book
site: bookdown::bookdown_site
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
github-repo: microbiome/OMA
---



---
date: "**Authors:** Leo Lahti [aut], Sudarshan Shetty [aut], Felix GM Ernst [aut, cre]<br/>
  **Version:** 0.98.9<br/>
  **Modified:** 2021-04-10<br/>
  **Compiled:** 2021-07-21<br/>
  **Environment:** R version 4.1.0 (2021-05-18), Bioconductor 3.14<br/>
  **License:** CC BY-NC-SA 3.0 US<br/>
  **Copyright:** <br/>
  **Source:** https://github.com/microbiome/OMA"
url: "https://github.com/microbiome/OMA"
---

# Preface {-}

This website is a book on microbiome analysis in the Bioconductor universe and
is showing common principles and workflows of performing microbiome analysis.

The book was borne out of necessity, while updating tools for microbiome
analysis to work with common classes of the Bioconductor project handling 
count data of various sorts. It is heavily influenced by similar resources, such
as the [Orchestrating Single-Cell Analysis with Bioconductor](https://www.nature.com/articles/s41592-019-0654-x)
book, [phyloseq tutorials](http://joey711.github.io/phyloseq/tutorials-index)
and [microbiome tutorials](https://microbiome.github.io/tutorials/).

We focus on microbiome analysis tools, new, updated and established methods.
In the *Introduction* section, we show how to work with the key data 
infrastructure `TreeSummarizedExperiment` and related classes, how this 
framework relates to other infrastructure and how to load microbiome analysis 
data to work with in the context of this framework.

The second section, *Focus Topics*, is all about the steps for analyzing
microbiome data, beginning with the most common steps and progressing to
more specialized methods in subsequent sections.

The third section, *Appendix*, contains the rest of things we didn't find 
another place for, yet.


