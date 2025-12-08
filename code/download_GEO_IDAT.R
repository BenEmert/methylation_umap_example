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
})

options(timeout = max(120, getOption("timeout"))) # 20 mins for large downloads

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
message("Fetching Supplementatry Table 4 (Koelsche et al 2021) from Hugging Face...")
meta_url <- "https://huggingface.co/datasets/bemert/GSE140686_GPL13534/resolve/main/41467_2020_20603_MOESM4_ESM.xlsx"
meta_dir <- file.path(data_dir, "raw", "GSE140686")
meta_file <- file.path(meta_dir, "41467_2020_20603_MOESM4_ESM.xlsx")

if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE)

tryCatch({
  download.file(meta_url, meta_file, mode = "wb", quiet = TRUE)
  message("Metadata downloaded successfully.")
}, error = function(e) {
  stop("Failed to download metadata from Hugging Face: ", e$message)
})

meta_df <- read_excel(meta_file)

# --- GEO METADATA ---
message("Fetching GSE140686 GEO records...")
gse_list <- getGEO("GSE140686", GSEMatrix = TRUE, getGPL = FALSE)

if (inherits(gse_list, "list")) {
  geo_meta <- do.call(bind_rows, lapply(gse_list, pData))
} else {
  geo_meta <- pData(gse_list)
}

geo_meta$basename_map <- str_extract(
  basename(as.character(geo_meta$supplementary_file)), 
  "\\d{10,12}_R\\d{2}C\\d{2}"
)

matched_data <- geo_meta %>% 
  filter(basename_map %in% meta_df$IDAT) %>%
  rename(DNA_Type = "material preparation:ch1") %>%
  select(geo_accession, supplementary_file, basename_map, platform_id, DNA_Type)

message(paste("Found", nrow(matched_data), "matching samples."))

download_idats <- function(url_list, dest_dir) {
  results_list <- vector("list", length(url_list))
  
  for (i in seq_along(url_list)) {
    url <- url_list[i]
    fn <- basename(url)
    dest <- file.path(dest_dir, fn)
    
    status <- "Unknown"
    
    if (file.exists(dest) && file.size(dest) > 0) {
      status <- "Skipped (Exists)"
    } else {
      try_result <- tryCatch({
        download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
        list(status = "Success", msg = NA_character_)
      }, error = function(e) {
        if(file.exists(dest)) unlink(dest)
        list(status = "Failed", msg = as.character(e$message))
      })
      status <- try_result$status
    }
    
    results_list[[i]] <- data.frame(Filename = fn, Status = status, URL = url)
    if (i %% 10 == 0) message(sprintf("Progress: %.1f%%", (i / length(url_list)) * 100))
  }
  return(do.call(rbind, results_list))
}

final_urls <- matched_data %>%
  mutate(
    url_root = str_remove(as.character(supplementary_file), "_(Red|Grn|Green)\\.idat(\\.gz)?$"),
    grn_suffix = if_else(str_detect(supplementary_file, "_Green"), "_Green.idat.gz", "_Grn.idat.gz")
  ) %>%
  mutate(Red = paste0(url_root, "_Red.idat.gz"), Grn = paste0(url_root, grn_suffix), .keep = "none") %>%
  pivot_longer(cols = everything(), values_to = "url") %>%
  pull(url)

message(paste("Downloading", length(final_urls), "IDAT files..."))
download_log <- download_idats(final_urls, idat_dir)

# --- RETRY LOGIC ---
failed_downloads <- download_log %>% filter(Status == "Failed")
if (nrow(failed_downloads) > 0) {
  message(sprintf("Retrying %d failed downloads...", nrow(failed_downloads)))
  retry_log <- download_idats(failed_downloads$URL, idat_dir)
  download_log <- download_log %>% filter(Status != "Failed") %>% bind_rows(retry_log)
}

final_failures <- download_log %>% filter(Status == "Failed")
num_final_failed <- nrow(final_failures)

message("\n--- Download Summary ---")
message(paste("Total Files:", nrow(download_log)))
message(paste("Successful:", sum(download_log$Status == "Success" | download_log$Status == "Skipped (Exists)")))
message(paste("Failed:", num_final_failed))

# --- TARGETS CREATION ---
targets <- matched_data %>%
  inner_join(meta_df, by = c("basename_map" = "IDAT")) %>%
  mutate(
    Basename = file.path(idat_dir, paste0(geo_accession, "_", basename_map)),
    Array_Platform = ifelse(platform_id == "GPL13534", "450k", "EPIC")
  ) %>%
  select(Basename, geo_accession, IDAT = basename_map, Array_Platform, DNA_Type = DNA, Diagnosis, Batch, Supplier)

write.csv(targets, file = file.path(data_dir, "raw/minfi_targets_reference.csv"), row.names = FALSE)
message("Download pipeline complete.")