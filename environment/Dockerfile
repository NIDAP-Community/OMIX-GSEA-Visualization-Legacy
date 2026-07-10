# hash:sha256:ce1882f58426497c7b638f61b8e26c5d3285c55be3bc4678e8ae66e46807c4b8
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/mosuite:v0.3.2

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("DT", "0.34.0")' \
    && Rscript -e 'remotes::install_version("RColorBrewer", "1.1-3")' \
    && Rscript -e 'remotes::install_version("htmlwidgets", "1.6.4")' \
    && Rscript -e 'remotes::install_version("optparse", "1.8.2")' \
    && Rscript -e 'remotes::install_version("patchwork", "1.3.2")' \
    && Rscript -e 'remotes::install_version("plotly", "4.12.0")'
