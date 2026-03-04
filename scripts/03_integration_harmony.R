\
#!/usr/bin/env Rscript

# 03_integration_harmony.R
# - Merge samples
# - Normalize, HVG, scale, PCA
# - Harmony integration by orig.ident
# - UMAP + clustering
# - Save integrated object and basic figures/tables

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
  library(harmony)
})

source("scripts/_helpers.R")
set_analysis_seed(1234)

# Input
objs <- readRDS("results/objects/02_inputs_filtered.rds")
K  <- objs$K
KP <- objs$KP
W  <- objs$W
WP <- objs$WP

# Add group labels
K$group  <- "K"
KP$group <- "KP"
W$group  <- "W"
WP$group <- "WP"

# Merge
all <- merge(K, y = c(KP, W, WP), project = "PIPseq")
DefaultAssay(all) <- "RNA"

# Parallel settings (safe defaults)
plan("multicore", workers = 4)
options(future.globals.maxSize = 36000 * 1024^2)

# Standard Seurat workflow
all <- NormalizeData(all)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(all)
all <- RunPCA(all, npcs = 30)

# Harmony integration
all$group <- factor(all$group, levels = c("W", "WP", "K", "KP"))
all@meta.data$group <- as.factor(all@meta.data$group)

all <- RunHarmony(
  object = all,
  reduction = "pca",
  group.by.vars = "orig.ident",
  assay.use = "RNA",
  reduction.save = "harmony"
)

# UMAP + clustering
all <- RunUMAP(all, reduction = "harmony", dims = 1:30)
all <- FindNeighbors(all, reduction = "harmony", dims = 1:30)
all <- FindClusters(all, resolution = 0.4)

# Figures
safe_pdf("results/figures/umap/UMAP_clustering_harmony.pdf", width = 6, height = 5)
DimPlot(all, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

safe_pdf("results/figures/umap/UMAP_clustering_harmony_bygroup.pdf", width = 18, height = 5)
DimPlot(all, reduction = "umap", label = TRUE, raster = FALSE, split.by = "group")
dev.off()

# Cells per cluster per sample
cells_per_cluster <- table(Idents(all), all$orig.ident)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.table(cells_per_cluster, "results/tables/number_of_cells_per_cluster.txt", sep = "\t", quote = FALSE)

# Save object
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
saveRDS(all, file = "results/objects/03_harmony_umap_clustered.rds")

save_session_info()
message("Done: results/objects/03_harmony_umap_clustered.rds")
