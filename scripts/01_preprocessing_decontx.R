\
#!/usr/bin/env Rscript

# 01_preprocessing_decontx.R
# - Read 10x matrices for each sample
# - Ambient RNA removal with celda::decontX
# - Create Seurat objects and downsample to equal cells/sample
# - Save intermediate objects

suppressPackageStartupMessages({
  library(Seurat)              # Seurat v5
  library(celda)               # decontX
  library(SingleCellExperiment)
})

source("scripts/_helpers.R")
set_analysis_seed(1234)

# Input directories (expected layout; see data/raw/README.md)
input_dirs <- list(
  K  = "data/raw/K/",
  KP = "data/raw/KP/",
  W  = "data/raw/W/",
  WP = "data/raw/WP/"
)

# Parameters
cells_per_sample <- 16000

read_decontx_seurat <- function(data_dir, project) {
  counts <- Read10X(data.dir = data_dir)
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- decontX(sce)
  seurat_obj <- CreateSeuratObject(round(decontXcounts(sce)), project = project)
  return(seurat_obj)
}

message("Reading data and running decontX...")
K  <- read_decontx_seurat(input_dirs$K,  "K")
KP <- read_decontx_seurat(input_dirs$KP, "KP")
W  <- read_decontx_seurat(input_dirs$W,  "W")
WP <- read_decontx_seurat(input_dirs$WP, "WP")

# Downsample to equal cell numbers
message("Downsampling to ", cells_per_sample, " cells per sample...")
K  <- K[,  sample(colnames(K),  size = cells_per_sample, replace = FALSE)]
KP <- KP[, sample(colnames(KP), size = cells_per_sample, replace = FALSE)]
W  <- W[,  sample(colnames(W),  size = cells_per_sample, replace = FALSE)]
WP <- WP[, sample(colnames(WP), size = cells_per_sample, replace = FALSE)]

# Quick checks
message("Cell counts per sample:")
print(table(K$orig.ident))
print(table(KP$orig.ident))
print(table(W$orig.ident))
print(table(WP$orig.ident))

# Save
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(K = K, KP = KP, W = W, WP = WP), file = "results/objects/01_inputs_decontx_downsampled.rds")

save_session_info()
message("Done: results/objects/01_inputs_decontx_downsampled.rds")
