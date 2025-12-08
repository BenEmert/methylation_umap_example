# Reproducible Methylation Analysis Pipeline (GSE140686)

This repository contains a fully reproducible pipeline for analyzing DNA methylation data from **GSE140686**. The pipeline handles data downloading (IDAT files), preprocessing, normalization, and dimensionality reduction (UMAP/t-SNE) to generate publication-quality figures.

The analysis is implemented in R and containerized with Docker to ensure exact reproducibility of results across different computing environments.

## Repository Structure

```text
.
├── Dockerfile              # Instructions to build the reproducible container
├── renv.lock               # Exact package versions used in the analysis
├── README.md               # This file
├── code/                   # Analysis scripts
│   ├── run_pipeline.R      # Master wrapper script
│   ├── download_GEO_IDAT.R # Data acquisition
│   ├── process_IDAT.R      # Normalization & filtering
│   ├── prepare_metadata.R  # Metadata cleaning
│   └── make_plots.R        # Figure generation
├── data/                   # (Auto-generated) Input and processed data
└── plots/                  # (Auto-generated) Final figure PDFs

```
1. Running the Pipeline Locally (No Docker)
If you prefer to run the analysis directly on your machine, you can use renv to restore the exact R package versions used in this project.

Prerequisites:
* R (version 4.3 or higher)
* Git
* Note for Linux users: You may need to install system libraries (e.g., libcurl4-openssl-dev, libxml2-dev) before installing R packages. See the Dockerfile for a reference list.

# Instructions
