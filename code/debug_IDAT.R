# code/debug_idat.R

# 1. The URL that failed (FTP)
ftp_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4180nnn/GSM4180927/suppl/GSM4180927_200325530180_R05C02_Grn.idat.gz"

# 2. Convert to HTTPS
https_url <- gsub("^ftp://", "https://", ftp_url)
dest_file <- "test_idat.idat.gz"

message("=== DEBUGGING IDAT DOWNLOAD ===")
message(paste("Original (FTP):", ftp_url))
message(paste("Testing (HTTPS):", https_url))

# 3. Attempt Download
tryCatch({
  download.file(https_url, dest_file, method = "libcurl", mode = "wb", quiet = TRUE)
  
  if (file.exists(dest_file) && file.size(dest_file) > 1000) {
    message("\n[PASS] HTTPS download successful!")
    message(paste("   File size:", file.size(dest_file), "bytes"))
    unlink(dest_file) # Cleanup
  } else {
    message("\n[FAIL] File downloaded but is empty.")
  }
  
}, error = function(e) {
  message("\n[FAIL] Download crashed.")
  message(paste("   Error:", e$message))
  
  if (grepl("403", e$message)) {
    message("   Note: If you see a 403 here, NCBI might still be blocking the User-Agent, even on HTTPS.")
  }
})

message("=== DEBUG COMPLETE ===")