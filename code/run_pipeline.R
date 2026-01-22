#!/usr/bin/env Rscript

# Master script to run the methylation analysis pipeline
# Location: code/run_pipeline.R
# Usage: Rscript code/run_pipeline.R (Run from Project Root)

# Define the directory where scripts live
# Since we are running from root, the scripts are in "code"
script_dir <- "code"

run_step <- function(script_name) {
  script_path <- file.path(script_dir, script_name)
  
  message(paste0("\n======================================================\n",
                 "RUNNING: ", script_name,
                 "\n======================================================"))
  
  if (!file.exists(script_path)) {
    stop(paste("Script not found:", script_path))
  }
  
  tryCatch({
    source(script_path, local = new.env())
    message(paste0("SUCCESS: ", script_name))
  }, error = function(e) {
    stop(paste0("FATAL ERROR in ", script_name, ":\n", e$message))
  })
}

# Run steps
run_step("download_GEO_IDAT.R")
run_step("process_IDAT.R")
run_step("prepare_metadata.R")
run_step("make_plots.R")
run_step("create_excel_summary.R")

message("\n======================================================\n",
        "PIPELINE COMPLETED SUCCESSFULLY\n",
        "Outputs available in /data/analyzed and /plots\n",
        "======================================================")