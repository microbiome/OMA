
# Contributors {-}

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


### *Chandler Ross, PhD* {-}

Chandler Ross is a PhD researcher at Turku Data science group research team.His research interest is mainly into Eco-Evolutionary Dynamics, Statistical Modelling and Time Series.

### *Chouaib Benchraka, Msc* {-}

Chouaib Benchraka is a research assistant at Turku data science research group.He is interested working with Machine Learning and Data Science researches.

### *Felix Ernst, PhD* {-}

Felic Ernst is among the developers who participated in the development of [miaverse](https://microbiome.github.io) which is open source project developing R/Bioc methods, benchmarking data, and educational material for microbiome research based on the SummarizedExperiment class and its derivatives.

### *Giulio Benedetti, Research Assistant* {-}

Guilio Benedetti is a research assistant at Turku data science research group. His research  interest  is mostly related to  Data Science.

### *Henrik Eckermann, PhD* {-}

Henerik Eckermann is PhD researcher at the research group Turku data science. His research interest is in Human microbiome bioinformatics.

### *Leo Lahti, Professor* {-}

Leo Lahti is the group leader for the research team at Turku data science and professor in Data Science in University of Turku, Finland. His research team focuses on computational analysis and modeling of complex natural and social systems. Lahti obtained doctoral degree (DSc) in statistical machine learning and bioinformatics from Aalto University in Finland (2010), developing probabilistic data integration methods for high-throughput life science data.

### *Shigdel Rajesh, PhD* {-}

Shigdel Rajesh was a visiting scholar to Turku data science group.His research interest is mainly related to bioinformatics, microbiome data analysis. 

### *Sudarshan Shetty, PhD* {-}

Sudarshan Shetty is one of the key developers in the [miaverse](https://microbiome.github.io) R/Bioconductor package for microbiome bioinformatics based on (Tree)SummerizedExperiment. 

### *Tuomas Borman, PhD* {-}

Tuomas Borman is a PhD researcher and one of the key developers of [miaverse](https://microbiome.github.io). He participated in setting up the base ecosystem of data containers, packages, and tutorials with other team members.

### *Yağmur Şimşek, *{-}

Yağmur Şimşek was a research assistant at the Turku data science reseach group when she contributed to this work.


### Acknowledgments {-}

This work would not have been possible without the countless
contributions and interactions with other researchers, developers, and
users. We express our gratitude to the entire Bioconductor community
for developing this high-quality open research software repository for
life science analytics, continuously pushing the limits in emerging
fields [@Gentleman2004, @Huber2015]. The developers and contributors
of this online tutorial are listed in Chapter \@ref(contributors).

The base ecosystem of data containers, packages, and tutorials was set
up as a collaborative effort by Tuomas Borman, Henrik Eckermann,
Chouaib Benchraka, Chandler Ross, Shigdel Rajesh, Yağmur Şimşek,
Giulio Benedetti, Sudarshan Shetty, Felix Ernst, and [Leo
Lahti](http://www.iki.fi/Leo.Lahti).

The work has been supported by the COST Action network on Statistical
and Machine Learning Techniques for Human Microbiome Studies
([ML4microbiome](https://www.ml4microbiome.eu/)) [@MorenoIndias2021].

The framework is based on the _TreeSummarizedExperiment_ data
container created by Ruizhu Huang and others
[@R_TreeSummarizedExperiment], and on the MultiAssayExperiment by
Marcel Ramos et al. [@Ramos2017]. The idea of using these containers
as a basis for microbiome data science was initially advanced by the
groundwork of Domenick Braccia, Héctor Corrada Bravo and others, and
subsequently brought together with other microbiome data science
developers [@Shetty2019].

Ample demonstration data resources have been made available as the
[curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/)
project by Edoardo Pasolli, Lucas Schiffer, Levi Waldron and others
[@Pasolli2017] adding important support.
A number of other contributors have advanced the ecosystem
further, and will be acknowledged in the individual
packages, [pull
requests](https://github.com/microbiome/OMA/graphs/contributors),
[issues](https://github.com/microbiome/OMA/issues), and other work.

The work has drawn inspiration from many sources, most notably from
the work on _phyloseq_ by Paul McMurdie and Susan Holmes
[@McMurdie2013] who pioneered the work on rigorous and reproducible
microbiome data science ecosystems in R/Bioconductor. The phyloseq
framework continues to provide a vast array of complementary packages
and methods for microbiome studies, and we aim to support full
interoperability.

The open source books by Susan Holmes and Wolfgang Huber, Modern
Statistics for Modern Biology [@Holmes2019] and by Garret Grolemund
and Hadley Wickham, the R for Data Science [@Grolemund2017], and
Richard McElreath's Statistical Rethinking and the associated online
resources by Solomon Kurz [@McElreath2020] are key references that
advanced reproducible data science training and dissemination. The
Orchestrating Single-Cell Analysis with Bioconductor, or _OSCA_ book
by Robert Amezquita, Aaron Lun, Stephanie Hicks, and Raphael Gottardo
[@Amezquita2020natmeth] has implemented closely related work on the
_SummarizedExperiment_ data container and its derivatives in the field
of single cell sequencing studies. Many approaches used in this book
have been derived from the [OSCA
framework](https://bioconductor.org/books/release/OSCA/), with various
adjustments and extensions dedicated to microbiome data science.

