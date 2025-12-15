suppressPackageStartupMessages({
  library(minfi)
  library(GEOquery)
  library(stringr)
  library(data.table)
  library(limma)
  library(sva)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(parallel)
})

options(timeout = max(300, getOption("timeout"))) # 5 mins timeout

# Relative Paths
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
idat_dir <- file.path(data_dir, "raw", "idats")
processed_dir <- file.path(data_dir, "analyzed")

# Create directories
sapply(c(data_dir, idat_dir, processed_dir), function(x) {
  if(!dir.exists(x)) dir.create(x, recursive = TRUE)
})

# --- DOWNLOAD METADATA (Hugging Face) ---
message("Fetching Supplementatry Table 4 (Koelsche et al 2021)")
meta_url <- "https://huggingface.co/datasets/bemert/GSE140686_GPL13534/resolve/main/41467_2020_20603_MOESM4_ESM.xlsx"
meta_dir <- file.path(data_dir, "raw", "GSE140686")
meta_file <- file.path(meta_dir, "41467_2020_20603_MOESM4_ESM.xlsx")

if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE)

if (!file.exists(meta_file)) {
  message(paste("Downloading:", meta_file))
  tryCatch({
    download.file(meta_url, meta_file, mode = "wb", quiet = TRUE)
    message("Metadata downloaded successfully.")
  }, error = function(e) {
    stop("Failed to download metadata from Hugging Face: ", e$message)
  })
} else {
  message("Metadata already exists. Skipping download")
}

meta_df <- read_excel(meta_file)

# --- GEO METADATA ---
message("Fetching GSE140686 Series Matrices from Hugging Face...")

base_hf_url <- "https://huggingface.co/datasets/bemert/GSE140686_GPL13534/resolve/main/"
matrix_files <- c(
  "GSE140686-GPL13534_series_matrix.txt.gz",
  "GSE140686-GPL21145_series_matrix.txt.gz"
)

meta_list <- list()

for (fn in matrix_files) {
  url <- paste0(base_hf_url, fn)
  dest_matrix <- file.path(meta_dir, fn)
  
  if (!file.exists(dest_matrix)) {
    message(paste("Downloading:", fn))
    tryCatch({
      download.file(url, dest_matrix, mode = "wb", quiet = TRUE, method = "libcurl")
    }, error = function(e) {
      stop(paste("Failed to download", fn, "from Hugging Face:", e$message))
    })
  }
  
  message(paste("Parsing:", fn))
  gse_obj <- getGEO(filename = dest_matrix, getGPL = FALSE)
  
  if (inherits(gse_obj, "list")) {
    meta_list[[fn]] <- pData(gse_obj[[1]])
  } else {
    meta_list[[fn]] <- pData(gse_obj)
  }
}

message("Merging platform metadata...")
geo_meta <- bind_rows(meta_list)

geo_meta$basename_map <- str_extract(
  basename(as.character(geo_meta$supplementary_file)), 
  "\\d{10,12}_R\\d{2}C\\d{2}"
)

matched_data <- geo_meta %>% 
  filter(basename_map %in% meta_df$IDAT) %>%
  rename(DNA_Type = "material preparation:ch1") %>%
  select(geo_accession, supplementary_file, basename_map, platform_id, DNA_Type)

message(paste("Found", nrow(matched_data), "matching samples."))

download_idats_parallel <- function(url_list, dest_dir) {
  
  num_cores <- if (.Platform$OS.type == "windows") 1 else max(1, parallel::detectCores() - 2)
  
  message(sprintf("Starting download using %d parallel cores...", num_cores))
  
  download_worker <- function(url) {
    fn <- basename(url)
    dest <- file.path(dest_dir, fn)
    
    if (file.exists(dest) && file.size(dest) > 0) {
      return(data.frame(Filename = fn, Status = "Skipped (Exists)", URL = url))
    }
    
    try_result <- tryCatch({    
      download.file(url, dest, mode = "wb", quiet = FALSE, method = "libcurl")
      "Success"
    }, error = function(e) {
      message(paste("Failed to download", url, ":", e$message))
      if(file.exists(dest)) unlink(dest)
      return("Failed")
    })
    
    if(try_result == "Success") cat(".") 
    
    return(data.frame(Filename = fn, Status = try_result, URL = url))
  }
  
  results_list <- parallel::mclapply(url_list, download_worker, mc.cores = num_cores)
  
  # Combine results
  results_list <- results_list[!sapply(results_list, is.null)]
  return(do.call(rbind, results_list))
}

final_urls <- matched_data %>%
  mutate(
    # Convert ftp:// to https://
    supplementary_file = str_replace(as.character(supplementary_file), "^ftp://", "https://"),
    
    url_root = str_remove(supplementary_file, "_(Red|Grn|Green)\\.idat(\\.gz)?$"),
    grn_suffix = if_else(str_detect(supplementary_file, "_Green"), "_Green.idat.gz", "_Grn.idat.gz")
  ) %>%
  mutate(
    Red = paste0(url_root, "_Red.idat.gz"), 
    Grn = paste0(url_root, grn_suffix), 
    .keep = "none"
  ) %>%
  pivot_longer(cols = everything(), values_to = "url") %>%
  pull(url)

message(paste("\nDownloading", length(final_urls), "IDAT files..."))
message(paste("Downloading", length(final_urls), "IDAT files..."))
download_log <- download_idats_parallel(final_urls, idat_dir)

# --- RETRY LOGIC ---
failed_downloads <- download_log %>% filter(Status == "Failed")
if (nrow(failed_downloads) > 0) {
  message(sprintf("\nRetrying %d failed downloads sequentially...", nrow(failed_downloads)))
  
  # Simple sequential retry for the stubborn files
  retry_log <- lapply(failed_downloads$URL, function(url) {
    fn <- basename(url)
    dest <- file.path(idat_dir, fn)
    tryCatch({
      download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
      data.frame(Filename = fn, Status = "Success", URL = url)
    }, error = function(e) {
      data.frame(Filename = fn, Status = "Failed", URL = url)
    })
  }) %>% bind_rows()
  
  download_log <- download_log %>% filter(Status != "Failed") %>% bind_rows(retry_log)
}

final_failures <- download_log %>% filter(Status == "Failed")
num_final_failed <- nrow(final_failures)

message("\n--- Download Summary ---")
message(paste("Total Files:", nrow(download_log)))
message(paste("Successful:", sum(download_log$Status == "Success" | download_log$Status == "Skipped (Exists)")))
message(paste("Failed:", num_final_failed))

if (num_final_failed > 0) {
  warning("The following files failed to download:")
  print(final_failures$Filename)
} else {
  message("\nAll downloads completed successfully.")
}

targets <- matched_data %>%
  inner_join(meta_df, by = c("basename_map" = "IDAT")) %>%
  mutate(
    Basename = file.path(idat_dir, paste0(geo_accession, "_", basename_map)),
    Array_Platform = ifelse(platform_id == "GPL13534", "450k", "EPIC")
  ) %>%
  select(Basename, geo_accession, IDAT = basename_map, Array_Platform, DNA_Type = DNA, Diagnosis, Methylation_Class = `Methylation Class Name`, Batch, Supplier, Color = Colour)

targets <- targets %>%
  mutate(Methylation_Class = str_remove(Methylation_Class, "methylation class "))

write.csv(targets, file = file.path(data_dir, "raw/minfi_targets_reference.csv"), row.names = FALSE)
message("Download pipeline complete.")
