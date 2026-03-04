\
#!/usr/bin/env Rscript

# 04_markers_and_annotation.R
# - Find cluster markers
# - Export marker tables (all + top10 per cluster)
# - Rename clusters to annotated cell types (as in the original notes)
# - Plot annotated UMAPs
# - Save annotated object

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source("scripts/_helpers.R")
set_analysis_seed(1234)

all <- readRDS("results/objects/03_harmony_umap_clustered.rds")

# Seurat v5 layers: join prior to marker tests
all <- JoinLayers(all)
DefaultAssay(all) <- "RNA"

# Markers
all.markers <- FindAllMarkers(all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.table(all.markers, "results/tables/all_cluster_marker_genes.csv", sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)

topmarkers <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

write.table(topmarkers, "results/tables/top10_cluster_marker_genes.csv", sep = ",", col.names = TRUE, row.names = TRUE, quote = FALSE)

# Cell type annotation (from the original notes)
new.cluster.ids <- c(
  "Macrophages",
  "CXCR2-high Neutrophils",
  "Atp1b1-high Fibroblasts",
  "Kdr-high BECs",
  "Dendritic cells",
  "T, B, NK cells",
  "Vwf-high BECs",
  "SMCs & Pericytes",
  "Fibroblasts",
  "Neutrophils",
  "CCR2-high Monocytes",
  "RBCs",
  "Monocytes",
  "T, B, NK cells",
  "Schwann cells",
  "Lyve1-high Macs",
  "Neurons",
  "LECs",
  "Plasma cells",
  "Mast cells",
  "Lyve1-high Macs",
  "Lyve1-high Macs",
  "Lyve1-high Macs"
)

names(new.cluster.ids) <- levels(all)
all <- RenameIdents(all, new.cluster.ids)

# Order levels for plotting/consistency
levels(all) <- c(
  "RBCs",
  "Schwann cells",
  "Neurons",
  "Mast cells",
  "T, B, NK cells",
  "Plasma cells",
  "CXCR2-high Neutrophils",
  "Neutrophils",
  "Monocytes",
  "CCR2-high Monocytes",
  "Dendritic cells",
  "Macrophages",
  "Lyve1-high Macs",
  "LECs",
  "Kdr-high BECs",
  "Vwf-high BECs",
  "SMCs & Pericytes",
  "Fibroblasts",
  "Atp1b1-high Fibroblasts"
)

# Plots (use default colors; journals often recolor in figure production anyway)
safe_pdf("results/figures/umap/UMAP_annotated.pdf", width = 10, height = 5)
DimPlot(all, reduction = "umap", label = TRUE, raster = FALSE)
dev.off()

safe_pdf("results/figures/umap/UMAP_annotated_bygroup.pdf", width = 30, height = 5)
DimPlot(all, reduction = "umap", label = TRUE, raster = FALSE, split.by = "group")
dev.off()

# Cells per annotated cell type per sample
cells_per_celltype <- table(Idents(all), all$orig.ident)
write.table(cells_per_celltype, "results/tables/cells_per_celltype_annotated.txt", sep = "\t", quote = FALSE)

# Save annotated object
dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
saveRDS(all, file = "results/objects/04_harmony_umap_annotated.rds")

save_session_info()
message("Done: results/objects/04_harmony_umap_annotated.rds")
