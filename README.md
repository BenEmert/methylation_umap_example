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
  -v $(pwd)/data:/home/methylation_umap_example/data \
  -v $(pwd)/plots:/home/methylation_umap_example/plots \
  bemert/methylation_umap_pipeline Rscript code/run_pipeline.R
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


## Interactive Visualization (Shiny App)

This repository includes a Shiny application to interactively explore the UMAP and t-SNE embeddings, allowing you to tune hyperparameters (e.g., neighbors, perplexity) and filter samples by diagnosis on your local computer.

**Prerequisites**
- The methylation analysis pipeline must be run first so that processed data exists in `data/analyzed/`.
- If running locally, ensure you have restored the R environment with `renv::restore()`.

**How to Run (Docker)**

To run the Shiny app inside the container, use the command below. 

*   `-p 3838:3838`: Maps the container's port to your local machine.
*   `-v $(pwd)/data:...`: **Crucial.** Mounts your local processed data so the app can read it.

```bash
docker run --rm -p 3838:3838 \
  -v $(pwd)/data:/home/methylation_umap_example/data \
  methylation_umap_pipeline \
  Rscript -e "shiny::runApp('shinyApp', port = 3838, host = '0.0.0.0')"
```

**How to Run (Locally)**
You can launch the app directly from an R console at the project root:

```r
shiny::runApp("shinyApp")
```

Alternatively, if using RStudio:
1. Open `shinyApp/app.R`
2. Click the **Run App** button in the source editor.

**Features**
- **Dual View**: Compare two different embedding settings (e.g., UMAP vs t-SNE) side-by-side.
- **Dynamic Tuning**: Adjust `n_neighbors`, `min_dist`, and `perplexity` on the fly.
- **Filtering**: Select specific diagnoses to highlight or subset.
- **Download**: Export the coordinate data for your custom views.


## Cloud Demo (Google Colab)

If you prefer not to install R or run the pipeline locally, you can explore the data using our Python-based dashboard on Google Colab. This notebook downloads the necessary data from the cloud and launches an interactive interface similar to the Shiny app.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1mhxDb12dY-lz1B4hIdR-VUvqNSNUu_mg?usp=sharing)

**Features:**
- **Zero Setup:** Runs entirely in your browser.
- **Interactive:** Uses `Panel` and `Plotly` to replicate the UMAP/t-SNE tuning.
- **Self-Contained:** Automatically fetches processed data.


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
