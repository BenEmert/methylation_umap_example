# Reproducible Methylation Analysis Pipeline (GSE140686)

This repository contains a fully reproducible pipeline for analyzing DNA methylation data from [GSE140686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140686). The pipeline handles data downloading (IDAT files), preprocessing, normalization, and dimensionality reduction (UMAP/t-SNE).

The analysis is implemented in R and containerized with Docker to ensure reproducibility of results across different computing environments.

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

## Running the Pipeline with Docker (Recommended)
Using Docker is the most reliable way to reproduce the analysis, as it encapsulates the operating system, system libraries, and R environment.

**Prerequisites**
- Docker Desktop (or Docker Engine on Linux)


1. **Getting the Docker Image**

Pull the pre-built image from Docker hub: 

```bash
docker pull bemert/methylation_umap_pipeline:latest
```

If you want to build the Docker image yourself, see 
[Building the Docker Image](#building-the-docker-image)

Once you have the image (either built locally or pulled from a registry), use the following command to run the analysis.

```
docker run --rm \
  -v $(pwd)/data:/home/rstudio/project/data \
  -v $(pwd)/plots:/home/rstudio/project/plots \
  methylation_umap_pipeline Rscript code/run_pipeline.R
```

**What this command does:**
- `--rm`: Automatically removes the container after the script finishes.
- `-v $(pwd)/data:...`: Maps your local data folder to the container
- `methylation_umap_pipeline`: The name of the Docker image.
- `Rscript code/run_pipeline.R`: The command to execute the pipeline.

## Running the Pipeline Locally (No Docker)
If you prefer to run the analysis directly on your machine, you can use renv to restore the exact R package versions used in this project.

**Prerequisites:**
- R (version 4.3 or higher)
- Git
- *Note for Linux users: You may need to install system libraries (e.g., libcurl4-openssl-dev, libxml2-dev) before installing R packages. See the Dockerfile for a reference list.*

1. **Clone the repository:**
    ```
    git clone [https://github.com/BenEmert/methylation_umap_example.git](https://github.com/BenEmert/methylation_umap_example.git)
    cd methylation_umap_example
    ```

2. **Restore the R environment:** Open R in the project root directory and run:
    ```
    if (!require("renv")) install.packages("renv")
    renv::restore()
    ```
    This will automatically install all required packages specified in renv.lock.

3. **Run the analysis:** You can execute the entire pipeline using the  wrapper script from your terminal:

    `Rscript code/run_pipeline.R`

    Alternatively, you can source `code/run_pipeline.R` from within an RStudio session.

## Outputs
Upon successful execution, the pipeline will populate the following directories:
- data/raw/: Contains downloaded IDAT files and metadata.
- data/analyzed/: Contains intermediate RDS files (normalized beta values, M-values) and the clean metadata CSV.
- plots/: Contains the final visualizations, including:
    - minfi_QC_plot.png
    - dx_legend.pdf
    - umap_GPL13534_M_random_state_comparison.pdf
    - umap_GPL13534_M_n_neighbors_comparison.pdf
    - tsne_GPL13534_M_perplexity_comparison.pdf


## Building the Docker Image
If you want to build the Docker image yourself (e.g., to verify the build process or modify the environment), follow these steps.

1. **Navigate to the project root:**
```
cd /path/to/YOUR_REPO_NAME
```

2. **Build the image:** Run the following command. This may take 20-40 minutes the first time as it compiles the R packages from source for Linux.

```
docker build -t methylation_umap_pipeline .
```

3. **Verify the build:** You can verify the image was created by running:

```
docker images
```

You should see methylation_umap_pipeline in the list. You can now proceed to step 2 ("Running the Pipeline with Docker").