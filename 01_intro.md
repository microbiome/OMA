# (PART) Introduction {-}

# Introduction {#introduction}

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




## Community

This online resource and its associated ecosystem of microbiome data
science tools are a result of a community-driven development process,
and welcoming to new contributors. You can find more information on
how to find us online and join the developer community through the
project homepage at
[microbiome.github.io](https://microbiome.github.io).



## Acknowledgements

This work would not have been possible without the countless
contributions and interactions within the broader research
community. In particular, we express our gratitude to the entire
Bioconductor community for developing this high-quality open research
software repository for life science analytics, continuously pushing
the limits in emerging application fields.

The base ecosystem of data containers, packages, and tutorials for was
set up as a collaborative effort by Tuomas Borman, Henrik Eckermann,
Chouaib Benchraka, Chandler Ross, Shigdel Rajesh, Yağmur Şimşek,
Giulio Benedetti, Sudarshan Shetty, Felix Ernst, [Leo
Lahti](http://www.iki.fi/Leo.Lahti). The framework is based on the
_TreeSummarizedExperiment_ data container created by Ruizhu Huang and
others [@R-TreeSummarizedExperiment], and the idea of using this as a
basis for microbiome data science was initially advanced by the
groundwork of Domenick Braccia, Héctor Corrada Bravo and others, and
brought together with other microbiome data science developers
[@Shetty2019]. Ample demonstration data resources have been made
available as the
[curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/)
project by Edoardo Pasolli, Lucas Schiffer, Levi Waldron and others
[@Pasolli2017] added support for this framework. A number of further
contributors have advanced the ecosystem, and are acknowledged in the
individual packages, [pull
requests](https://github.com/microbiome/OMA/graphs/contributors),
[issues](https://github.com/microbiome/OMA/issues), and other work.

The work has drawn inspiration from many sources, most
notably from the work on _phyloseq_ by Paul McMurdie and Susan Holmes
[@McMurdie2013] who have pioneered the work on rigorous and
reproducible microbiome data science in R/Bioconductor. The phyloseq
framework continues to provide a strong of complementary packages and
methods in this field, and we thrive to support full
interoperability. Open source books by Susan Holmes and Wolfgang
Huber, Modern Statistics for Modern Biology [@Holmes2019] and by
Garret Grolemund and Hadley Wickham, the R for Data Science
[@Grolemund2017], and Richard McElreath's Statistical Rethinking and
the associated online resources by Solomon Kurz [@McElreath2020] are
key references that advanced reproducible data science training and
dissemination. The Orchestrating Single-Cell Analysis with
Bioconductor, or _OSCA_ book by Robert Amezquita, Aaron Lun, Stephanie
Hicks, and Raphael Gottardo [@Amezquita2020] has implemented closely
related work on the _SummarizedExperiment_ data container and its
derivatives in the field of single cell sequencing studies. Many
approaches used in this book have been derived from the OSCA
framework, with the necessary adjustments for microbiome research.

<<<<<<< HEAD
- Getting started with R? Try [swirl](https://swirlstats.com/)

- [Orchestrating Microbiome Analysis in R and Bioconductor](microbiome.github.io/OMA)

=======
>>>>>>> origin/master








