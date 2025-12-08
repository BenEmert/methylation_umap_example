url_to_test <- "https://huggingface.co/datasets/bemert/GSE140686_GPL13534/resolve/main/GSE140686-GPL13534_series_matrix.txt.gz"
dest_file <- "test_hf_matrix.txt.gz"

message("=== DEBUGGING HUGGING FACE DOWNLOAD ===")

# 2. Test Download (Using libcurl)
message("\n2. Testing download from Hugging Face...")
tryCatch({
  # We use libcurl here as it's the standard for your pipeline
  download.file(url_to_test, dest_file, method = "libcurl", mode = "wb", quiet = TRUE)
  
  if (file.exists(dest_file) && file.size(dest_file) > 1000) {
    message("   [PASS] Download successful!")
    message(paste("   File size:", file.size(dest_file), "bytes"))
    
    # Optional: Read first few lines to verify it's text/gzip content and not an error HTML page
    con <- gzfile(dest_file, "rt")
    first_line <- readLines(con, n = 1)
    close(con)
    message(paste("   First line verification:", substr(first_line, 1, 50), "..."))
    
    unlink(dest_file)
  } else {
    message("   [FAIL] File downloaded but seems empty or too small.")
  }
  
}, error = function(e) {
  message("   [FAIL] Download crashed.")
  message(paste("   Error:", e$message))
  
  # Check for 404 specifically (File not found / Typo in URL)
  if (grepl("404", e$message)) {
    message("   HINT: A 404 error usually means the file hasn't been uploaded to the repo yet, or the filename in the script doesn't match exactly.")
  }
})

message("\n=== DEBUG COMPLETE ===")