suppressPackageStartupMessages({
  library(minfi)
  library(limma)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(dplyr)
  library(stringr)
  library(data.table)
  library(readxl)
})

project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
processed_dir <- file.path(data_dir, "analyzed")
targets_file <- file.path(data_dir, "raw/minfi_targets_reference.csv")

if(!file.exists(targets_file)) stop("Targets file not found. Run download script first.")
targets <- read.csv(targets_file)

message("Loading data by platform...")

targets_450k <- targets[targets$Array_Platform == "450k", ]
targets_EPIC <- targets[targets$Array_Platform == "EPIC", ]

message(sprintf("%d 450k samples", nrow(targets_450k)))
message(sprintf("%d EPIC samples", nrow(targets_EPIC)))
message( "Loading data by platform...")

RGset_450k <- read.metharray.exp(targets = targets_450k, verbose = TRUE)
RGset_EPIC <- read.metharray.exp(targets = targets_EPIC, verbose = TRUE, force = TRUE)

message("Merging platforms into Virtual 450k Array...")
RGset_Virtual <- combineArrays(RGset_450k, RGset_EPIC, outType = "IlluminaHumanMethylation450k")

# Clean up memory
rm(RGset_450k, RGset_EPIC); gc()

# ---------------------------------------------------------
# Probe QC
# ---------------------------------------------------------
message("Calculating QC Metrics...")
detP <- detectionP(RGset_Virtual)
frac_detP_pass <- colMeans(detP < 0.05) 

# Median Log2 Intensities 
MSet_Raw <- preprocessRaw(RGset_Virtual)
qc_metrics <- getQC(MSet_Raw)

# Apply Filtering Criteria
pass_qc <- (frac_detP_pass > 0.95) & 
           (qc_metrics$mMed > 8) & 
           (qc_metrics$uMed > 8)

message(paste("Samples passing QC:", sum(pass_qc), "/", length(pass_qc)))

plot_dir <- file.path(project_dir, "plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

png(file.path(plot_dir, "minfi_QC_plot.png"))
plot(qc_metrics$mMed, qc_metrics$uMed, 
     xlim=c(6,14), ylim=c(6,14), pch=16, cex=0.5,
     xlab="Median Log2 Meth", ylab="Median Log2 Unmeth",
     main="QC: Median Intensities")
abline(v = 8, h = 8, col = "red", lty = 2)
points(qc_metrics$mMed[!pass_qc], qc_metrics$uMed[!pass_qc], col="red", pch=16)
dev.off()

#By visual inspection, the 1 sample that doesn't pass filtering criteria
#because frac_detP_pass is 0.948861, groups with the rest of the samples
#This sample should be kept. Round to the nearest hundredth to let it pass. 
#There is a separate group of samples with abnormal ratio of log2(meth) to log2(unmeth)
#although these sample could be removed, I will leave them in to be consistent with the published studies

frac_detP_pass_rounded <- round(frac_detP_pass, 2) 

pass_qc_rounded <- (frac_detP_pass_rounded >= 0.95) & 
           (qc_metrics$mMed > 8) & 
           (qc_metrics$uMed > 8)

message(paste("Samples passing QC after rounding frac_detP_pass to hundreth decimal :", sum(pass_qc_rounded), "/", length(pass_qc_rounded)))

# Filter the original RGset
RGset_Filtered <- RGset_Virtual[, pass_qc_rounded]
targets <- targets[pass_qc_rounded, ]
detP <- detP[, pass_qc_rounded]

# Clean up memory (MSet_Raw is no longer needed)
rm(MSet_Raw); gc()

# ---------------------------------------------------------
# Normalization: "GenomeStudio-Style"
# ---------------------------------------------------------
message("Running 'Illumina' Normalization...")
MSet_Norm <- preprocessIllumina(RGset_Filtered, bg.correct = TRUE, normalize = "controls")

message("Mapping to Genome (MethylSet -> GenomicMethylSet)...")
GMSet_Norm <- mapToGenome(MSet_Norm)
GMSet_Norm <- addSnpInfo(GMSet_Norm)

rm(MSet_Norm, RGset_Filtered); gc()

# Match probeset from detP with the Normalized Object
# GMSet_Norm may have fewer probes due to preprocessIllumina dropping some probes
common_probes <- intersect(rownames(GMSet_Norm), rownames(detP))

# Subset both objects to the intersection
GMSet_Norm <- GMSet_Norm[common_probes, ]
detP <- detP[common_probes, ]

# Probe Filtering
# ---------------------------------------------------------
message("Filtering probes...")

# Remove probe if it fails (p > 0.05) in > 5% of samples
GMSet_Norm <- GMSet_Norm[rowMeans(detP > 0.05) < 0.05, ]
rm(detP); gc()

# Drop loci with SNPs
GMSet_Norm <- dropLociWithSnps(GMSet_Norm, snps = c("CpG", "SBE"))

# Drop sex Chromosomes
ann <- getAnnotation(GMSet_Norm)
keep_autosomes <- !(ann$chr %in% c("chrX", "chrY", "X", "Y"))
GMSet_Norm <- GMSet_Norm[keep_autosomes, ]

# Remove Cross-Reactive Probes 
# Download the Chen et al. list. Available from https://epigen.ccm.sickkids.ca/
message("Downloading Chen et al. blacklist...")
blacklist_url <- "https://epigen.ccm.sickkids.ca/probes/cross-reactive%20probe_Chen%20YA_Epigenetics.xlsx"
dest_file <- file.path(data_dir, "raw", "chen_blacklist.xlsx")

if (!file.exists(dest_file)) {
  tryCatch({
    download.file(blacklist_url, dest_file, mode = "wb", quiet = TRUE, method = "libcurl")
  }, error = function(e) {
    warning("Could not download blacklist. Proceeding without this filter.")
  })
}

if (file.exists(dest_file)) {
    cg_blacklist <- read_excel(dest_file, sheet = 'cg## cross-reactive targets', .name_repair = "minimal")
    ch_blacklist <- read_excel(dest_file, sheet = 'ch## cross-reactive targets', .name_repair = "minimal")
    #Filter for probes not mapping uniquely to the human genome allowing for up to one mismatch
    cg_blacklist <- cg_blacklist %>% 
        filter(`49` > 0 | `50` > 0)
    ch_blacklist <- ch_blacklist %>% 
        filter(`49` > 0 | `50` > 0)
    chen_blacklist_ids <- c(cg_blacklist$TargetID, ch_blacklist$TargetID)
    
    keep_unique <- !(rownames(GMSet_Norm) %in% chen_blacklist_ids)
    GMSet_Norm <- GMSet_Norm[keep_unique, ]
    message(paste("Removed", sum(!keep_unique), "cross-reactive probes."))
}

# ---------------------------------------------------------
# Batch Correction: Separate Channel
# ---------------------------------------------------------
message("Performing Separate Channel Batch Correction...")

# Prepare Design Matrix
targets$Batch_Correction_Type <- targets$DNA_Type
targets$Batch_Correction_Type[targets$Batch_Correction_Type == "EDTA blood"] <- "KRYO"
targets$Diagnosis <- as.factor(targets$Diagnosis)
design <- model.matrix(~Diagnosis, data = targets)
batch1 <- targets$Array_Platform
batch2 <- targets$Batch_Correction_Type

# temp file paths
temp_M_file <- file.path(processed_dir, "temp_M_corrected.rds")
temp_U_file <- file.path(processed_dir, "temp_U_corrected.rds")

message("Processing Methylated Channel (M)...")
M_matrix <- getMeth(GMSet_Norm)
log2_M <- log2(M_matrix + 1)
rm(M_matrix); gc() 

log2_M_corrected <- removeBatchEffect(log2_M, batch=batch1, batch2=batch2, design=design)

# Save result to disk
rm(log2_M); gc() 
saveRDS(log2_M_corrected, temp_M_file)
rm(log2_M_corrected); gc() 

message("Processing Unmethylated Channel (U)...")
U_matrix <- getUnmeth(GMSet_Norm)
rm(GMSet_Norm); gc()

log2_U <- log2(U_matrix + 1)
rm(U_matrix); gc()

log2_U_corrected <- removeBatchEffect(log2_U, batch=batch1, batch2=batch2, design=design)

# Save result to disk
rm(log2_U); gc()
saveRDS(log2_U_corrected, temp_U_file)
rm(log2_U_corrected); gc()

message("Calculating Final Betas...")

# Load corrected values
log2_M_corrected <- readRDS(temp_M_file)
M_corrected <- 2^log2_M_corrected - 1
M_corrected[M_corrected < 0] <- 0
rm(log2_M_corrected); gc()
# unlink(temp_M_file) # Delete temp file

log2_U_corrected <- readRDS(temp_U_file)
U_corrected <- 2^log2_U_corrected - 1
U_corrected[U_corrected < 0] <- 0
rm(log2_U_corrected); gc()
# unlink(temp_U_file)

# Calculate Betas
denom <- M_corrected + U_corrected + 100
Betas <- M_corrected / denom
# Makerows = samples, columns = probes
Betas <- t(Betas)
rm(denom, U_corrected, M_corrected); gc()

alpha <- 0.001 
M_vals <- log2((Betas + alpha) / (1 - Betas + alpha))

# --- SAVING ---
saveRDS(Betas, file.path(processed_dir, "GSE140686_Final_Betas_corrected.rds"))
saveRDS(M_vals, file.path(processed_dir, "GSE140686_Final_Mvals.rds"))
fwrite(targets, file.path(processed_dir, "GSE140686_Final_Targets.csv"), row.names = FALSE)

message("Processing complete.")