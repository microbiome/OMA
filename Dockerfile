ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}
COPY . /opt/pkg

# The devel image pins CRAN to 2026-10-28, before treelayout and tidytree
# 0.4.9 that ggtree 4.3.2 needs (Bioconductor/bioconductor_docker#146).
# Remove when the image's CRAN snapshot is newer.
RUN Rscript -e 'install.packages(c("treelayout", "tidytree"), repos = "https://cloud.r-project.org")'

# Install book package 
RUN Rscript -e 'repos <- BiocManager::repositories(); options(repos = repos); devtools::install("/opt/pkg", dependencies = TRUE, build_vignettes = FALSE); sessioninfo::session_info(installed.packages()[,"Package"], include_base = TRUE)'

## Build/install using same approach than BBS
RUN R CMD INSTALL /opt/pkg
RUN R CMD build --keep-empty-dirs --no-resave-data /opt/pkg
