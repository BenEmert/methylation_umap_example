# 1. Base Image
FROM rocker/r-ver:4.3.2

# 2. Install System Dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# 3. Set Project Root
WORKDIR /home/methylation_umap_example

# 4. Install Package Managers
RUN R -e "install.packages('BiocManager'); install.packages('renv')"

# 5. Restore R Environment
COPY renv.lock .
RUN RENV_PATHS_CACHE=/tmp/renv_cache RENV_CONFIG_SANDBOX_ENABLED=FALSE \
    R -e "renv::restore(prompt=FALSE)" \
    && rm -rf /tmp/renv_cache

# Clean up R documentation, tests, and temporary files to save space
RUN rm -rf /tmp/downloaded_packages \
    && strip /usr/local/lib/R/site-library/*/libs/*.so \
    && rm -rf /usr/local/lib/R/site-library/*/help \
    && rm -rf /usr/local/lib/R/site-library/*/doc \
    && rm -rf /usr/local/lib/R/site-library/*/tests

# 6. Copy Code and Create Directories
COPY code/ code/
COPY shinyApp/ shinyApp/
RUN mkdir -p data/raw/idats data/analyzed plots

# 7. Default Command
CMD ["/bin/bash"]
