# Base image with Bioconductor 3.18 (R 4.3)
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system libraries
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory to the PROJECT ROOT
WORKDIR /home/methylation_umap_example

# Copy lockfile and restore R packages
COPY renv.lock .
RUN Rscript -e "install.packages('renv'); renv::restore(prompt=FALSE)"

COPY code/ code/

# 3. Create empty directories for mount points
RUN mkdir -p data/raw/idats data/analyzed plots

# 4. Default command
CMD ["/bin/bash"]