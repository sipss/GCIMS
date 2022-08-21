FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_15

RUN apt-get update && apt-get install -y \
    qpdf \
    ghostscript \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra

RUN Rscript -e 'remotes::install_cran(c("rcmdcheck", "BiocManager", "covr", "pkgdown", "DT"))' && \
    Rscript -e 'BiocManager::install("BiocCheck")'

COPY DESCRIPTION /workdir/DESCRIPTION
WORKDIR /workdir

RUN Rscript -e 'remotes::install_deps(dependencies = TRUE)'
