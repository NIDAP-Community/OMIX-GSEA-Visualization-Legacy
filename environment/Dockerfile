# hash:sha256:d873abd175858260ef7960f560107b263dea4c08c96e3c860ad8fbd63a1b5f36
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/omix-r-pathway:r4.4.3-bioconductor3.20-v2

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "ComplexHeatmap" \
    ))' # Original versions: 2.22.0
