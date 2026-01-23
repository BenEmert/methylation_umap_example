# code/create_small_data.R

suppressPackageStartupMessages({
  library(matrixStats)
  library(dplyr)
  library(arrow) 
})

# --- CONFIGURATION ---
N_TOP_VAR <- 50000 

project_dir <- getwd()
data_dir <- file.path(project_dir, "data", "analyzed")

message("Loading full M-values for variance calculation...")
file_mvals <- file.path(data_dir, "GSE140686_Final_Mvals.rds")

if (file.exists(file_mvals)) {
  mvals <- readRDS(file_mvals)

  message(sprintf("Calculating variance for %d probes...", ncol(mvals)))
  probe_vars <- colVars(mvals)

  top_indices <- order(probe_vars, decreasing = TRUE)[1:min(N_TOP_VAR, length(probe_vars))]
  top_probes <- colnames(mvals)[top_indices]

  message(sprintf("Selected top %d variable probes.", length(top_probes)))

  # Create filtered M-values Parquet
  message("Creating filtered M-values dataframe...")
  mvals_lite <- mvals[, top_probes]

  df_mvals <- as.data.frame(mvals_lite)
  df_mvals$ID <- rownames(mvals_lite)
  df_mvals <- df_mvals %>% select(ID, everything())

  outfile_m <- file.path(data_dir, "GSE140686_Final_Mvals_50k_probes.parquet")
  write_parquet(df_mvals, outfile_m)
  message(paste("Saved:", outfile_m))

  # Clear memory
  rm(mvals, mvals_lite, df_mvals); gc()
} else {
   warning("M values file not found. Run process_IDAT.R first. Skipping.")
   stop()
}
exists("mva")
# Create filtered Beta-values Parquet
# Using the same set of probes chosen based on variance in M-values
# Selecting probes based on variance in Beta-values, bounded [0,1],
# might bias towards probes with average values near ~0.5. 
message("Loading full Beta-values...")
file_betas <- file.path(data_dir, "GSE140686_Final_Betas_corrected.rds")

if (file.exists(file_betas) && exists("top_probes")) {
  betas <- readRDS(file_betas)
  
  message("Subsetting Beta-values...")
  betas_lite <- betas[, top_probes]
  
  df_betas <- as.data.frame(betas_lite)
  df_betas$ID <- rownames(betas_lite)
  df_betas <- df_betas %>% select(ID, everything())
  
  outfile_b <- file.path(data_dir, "GSE140686_Final_Betas_corrected_50k.parquet")
  write_parquet(df_betas, outfile_b)
  message(paste("Saved:", outfile_b))
  
  rm(betas, betas_lite, df_betas); gc()
} else {
  warning("Beta values file not found. Skipping.")
}

message("Done.")