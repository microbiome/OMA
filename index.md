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
  **Modified:** 2021-09-10<br/>
  **Compiled:** 2021-10-05<br/>
  **Environment:** R version 4.1.1 (2021-08-10), Bioconductor 3.14<br/>
  **License:** CC BY-NC-SA 3.0 US<br/>
  **Copyright:** <br/>
  **Source:** https://github.com/microbiome/OMA"
url: "https://github.com/microbiome/OMA"
---

# Welcome {-}


<a href="https://bioconductor.org"><img src="https://github.com/Bioconductor/BiocStickers/raw/master/Bioconductor/Bioconductor-serial.gif" width="200" alt="Bioconductor Sticker" align="right" style="margin: 0 1em 0 1em" /></a>


You are reading the online book, **Orchestrating Microbiome Analysis
with Bioconductor** [@OMA], where we walk through common strategies and
workflows in microbiome data science.

The book shows through concrete examples how you can use some of the
latest advances in the R/Bioconductor framework for microbiome data
manipulation, analysis, and reproducible reporting. The book was borne
out of necessity, while updating microbiome analysis tools to work
with Bioconductor classes that provide support for multi-modal data
collections. These techniques are generic and widely applicable in
other contexts, too.

This work is heavily influenced by similar resources, such as the
Orchestrating Single-Cell Analysis with Bioconductor [@Amezquita2020natmeth],
[phyloseq tutorials](http://joey711.github.io/phyloseq/tutorials-index) [@Callahan2016] and
[microbiome tutorials](https://microbiome.github.io/tutorials/) [@Shetty2019].
This book extends the previous efforts of such related work, teaching the grammar
of Bioconductor workflows and thus supporting the
adoption of general data analysis skills in the analysis of large,
hierarchical multi-modal data collections.

We focus on microbiome analysis tools, new, updated and established methods.
In the *Introduction* section, we show how to work with the key data 
infrastructure `TreeSummarizedExperiment` and related classes, how this 
framework relates to other infrastructure and how to load microbiome analysis 
data to work with in the context of this framework.

This book is organized into three parts. We start by introducing the
material and link to further resources for learning R and
Bioconductor. We will describe the microbiome data containers, the
TreeSummarizedExperiment class. The second section, *Focus Topics*, is
all about the steps for analyzing microbiome data, beginning with the
most common steps and progressing to more specialized methods in
subsequent sections. Third, *Workflows*, provides case studies for the
various datasets used throughout the book. Finally, *Appendix*, links
to further resources.


--------------

The book is written in RMarkdown with the bookdown R package. OMA is a
collaborative effort, where several individuals have contributed
methods, workflows and improvements as acknowledged in the Appendix.
This online resource is **free to use** with the [Creative Commons
Attribution-NonCommercial
3.0](https://creativecommons.org/licenses/by-nc/3.0/us/) License.




<script type="text/javascript">
// This block adds image to the front page
title=document.getElementById('header');
title.innerHTML = title.innerHTML + 

'<img src="https://user-images.githubusercontent.com/60338854/128359392\
-6feef8df-30e9-4ea0-ae3b-4bb619d746ed.png" alt="Microbiome" width="50%"/>' +

'<p style="font-size:12px">Figure source: Moreno-Indias <i>et al</i>. (2021) \
<a href="https://doi.org/10.3389/fmicb.2021.635781">Statistical and \
Machine Learning Techniques in Human Microbiome Studies: Contemporary \
Challenges and Solutions</a>. Frontiers in Microbiology 12:11.</p>'
</script>
