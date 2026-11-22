ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}
# Copy DESCRIPTION alone first so that the dependency layers stay cached
# until DESCRIPTION changes.
COPY DESCRIPTION /opt/pkg/DESCRIPTION

# The devel image pins CRAN to 2026-10-28, before treelayout and tidytree
# 0.4.9 that ggtree 4.3.2 needs (Bioconductor/bioconductor_docker#146).
# Remove when the image's CRAN snapshot is newer.
RUN Rscript -e 'install.packages(c("treelayout", "tidytree"), repos = "https://cloud.r-project.org")'

# Install the book's dependencies
RUN Rscript -e 'repos <- BiocManager::repositories(); options(repos = repos); devtools::install_deps("/opt/pkg", dependencies = TRUE); sessioninfo::session_info(installed.packages()[,"Package"], include_base = TRUE)'

COPY . /opt/pkg

## Build/install using same approach than BBS
RUN R CMD INSTALL /opt/pkg
RUN R CMD build --keep-empty-dirs --no-resave-data /opt/pkg
