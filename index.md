--- 
title: "Microbiome Analysis"
documentclass: book
site: bookdown::bookdown_site
bibliography: [book.bib, packages.bib]
biblio-style: apalike
link-citations: yes
github-repo: microbiome/MiaBook
---



---
date: "**Authors:** Leo Lahti [aut], Sudarshan Shetty [aut], Felix GM Ernst [aut, cre]<br/>
  **Version:** 0.98.0003<br/>
  **Modified:** 2020-12-06<br/>
  **Compiled:** 2020-12-10<br/>
  **Environment:** R version 4.0.3 (2020-10-10), Bioconductor 3.12<br/>
  **License:** CC BY-NC-SA 3.0 US<br/>
  **Copyright:** <br/>
  **Source:** https://github.com/microbiome/MiaBook"
url: "https://github.com/microbiome/MiaBook"
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


