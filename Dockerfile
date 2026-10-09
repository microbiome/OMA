ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}
COPY . /opt/pkg

# Install book package 
RUN Rscript -e 'repos <- BiocManager::repositories(); options(repos = repos); devtools::install("/opt/pkg", dependencies = TRUE, build_vignettes = FALSE); sessioninfo::session_info(installed.packages()[,"Package"], include_base = TRUE)'

## Build/install using same approach than BBS
RUN R CMD INSTALL /opt/pkg
RUN R CMD build --keep-empty-dirs --no-resave-data /opt/pkg
