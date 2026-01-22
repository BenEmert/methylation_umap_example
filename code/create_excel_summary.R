# create_excel_summary.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(openxlsx)
})

# Define Paths
# ------------------------------------------------------------------------------
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
processed_dir <- file.path(data_dir, "analyzed")

input_meta <- file.path(processed_dir, "metadata_complete.csv")
output_file <- file.path(processed_dir, "supplementary_file1.xlsx")

# Main Script
# ------------------------------------------------------------------------------
message("=== Generating Excel Summary ===")

# Load Data
message("   LOADING Metadata...")
if (!file.exists(input_meta)) stop(paste("Metadata not found at:", input_meta))
metadata <- fread(input_meta, data.table = FALSE)

# Create Summary
message("   Creating summary data...")
summary_data <- metadata %>%
  count(Dx, Diagnosis, WHO_differentiation, Color_dx, name = "count") %>%
  arrange(Dx) %>%
  relocate(count, .before = Color_dx)

# Create Excel Workbook
message("   Creating Excel workbook...")
wb <- createWorkbook()

# Add Sheet 1: Full metadata
addWorksheet(wb, "sample_metadata")
writeData(wb, "sample_metadata", metadata)

# Add Sheet 2: Summary
addWorksheet(wb, "color_key")
writeData(wb, "color_key", summary_data)

# Apply cell styling
message("   Applying cell colors...")
for (i in 1:nrow(summary_data)) {
  # Style for Color_dx
  style_dx <- createStyle(fgFill = summary_data$Color_dx[i])
  addStyle(
    wb,
    "color_key",
    style = style_dx,
    rows = i + 1,
    cols = which(names(summary_data) == "Color_dx")
  )

  # Style for Color_WHO
  style_who <- createStyle(fgFill = summary_data$Color_WHO[i])
  addStyle(
    wb,
    "Diagnosis_Summary",
    style = style_who,
    rows = i + 1,
    cols = which(names(summary_data) == "Color_WHO")
  )
}

# Save Workbook
saveWorkbook(wb, output_file, overwrite = TRUE)

message(paste("=== Successfully created Excel file:", output_file, "==="))
