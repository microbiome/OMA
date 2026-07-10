ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}
COPY . /opt/pkg

# Install system Python, pip, and venv via Ubuntu's package manager – needed by mofa and book packages
RUN apt-get update && apt-get install -y python3.12 python3.12-venv python3-pip && Rscript -e 'if(!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_local(path = "/opt/pkg/", dependencies=TRUE, build_vignettes=FALSE, upgrade=TRUE)' && rm -rf /var/lib/apt/lists/*

RUN pip3 install mofapy2==0.7.3 --break-system-packages

# Force R/reticulate to use this specific Python installation
ENV RETICULATE_PYTHON="/usr/bin/python3.12"

## Build/install using same approach than BBS
RUN R CMD INSTALL /opt/pkg
RUN R CMD build --keep-empty-dirs --no-resave-data /opt/pkg
