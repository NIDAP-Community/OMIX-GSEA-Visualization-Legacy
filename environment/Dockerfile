# hash:sha256:7324ace058c103f46134c25023f1e323194e17f8d5e129c4ab255da029a5eac9
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/omix-r-pathway:r4.4.3-bioconductor3.20-v1

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "ComplexHeatmap" \
    ))' # Original versions: 2.22.0
