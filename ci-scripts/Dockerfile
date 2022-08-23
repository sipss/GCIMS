FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_15

RUN apt-get update && apt-get install -y \
    libmagick++-dev \
    qpdf \
    ghostscript \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-extra

COPY DESCRIPTION /workdir/DESCRIPTION
WORKDIR /workdir

RUN Rscript -e 'BiocManager::install(c("rcmdcheck", "BiocManager", "covr", "pkgdown", "DT", "magick", "BiocCheck", remotes::dev_package_deps(".", dependencies = TRUE)$package))'
