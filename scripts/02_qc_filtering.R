\
#!/usr/bin/env Rscript

# 02_qc_filtering.R
# - Add mitochondrial percentage
# - QC violin plots per sample
# - Filter cells (nFeature_RNA and percent.mt)
# - Save filtered objects

suppressPackageStartupMessages({
  library(Seurat)
})

source("scripts/_helpers.R")
set_analysis_seed(1234)

in_file <- "results/objects/01_inputs_decontx_downsampled.rds"
objs <- readRDS(in_file)

K  <- objs$K
KP <- objs$KP
W  <- objs$W
WP <- objs$WP

# Parameters (from the original notes)
min_features <- 200
max_features <- 6000
max_percent_mt <- 10
mt_pattern <- "^mt-"   # mouse mitochondrial genes

qc_one <- function(obj, sample_name) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)

  safe_pdf(sprintf("results/figures/qc/%s_QCVlnPlot.pdf", sample_name), width = 10, height = 5)
  VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()

  obj <- subset(
    obj,
    subset = nFeature_RNA > min_features & nFeature_RNA < max_features & percent.mt < max_percent_mt
  )
  return(obj)
}

message("QC plotting + filtering...")
K  <- qc_one(K,  "K")
KP <- qc_one(KP, "KP")
W  <- qc_one(W,  "W")
WP <- qc_one(WP, "WP")

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(K = K, KP = KP, W = W, WP = WP), file = "results/objects/02_inputs_filtered.rds")

save_session_info()
message("Done: results/objects/02_inputs_filtered.rds")
