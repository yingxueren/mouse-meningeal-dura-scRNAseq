\
#!/usr/bin/env Rscript

# 05_differential_expression.R
# - Differential expression per cell type between conditions using MAST
# - Outputs one CSV per cell type per comparison

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
})

source("scripts/_helpers.R")
set_analysis_seed(1234)

all <- readRDS("results/objects/04_harmony_umap_annotated.rds")

# Parallel settings
plan("multicore", workers = 4)
options(future.globals.maxSize = 36000 * 1024^2)

# Join layers before DE (Seurat v5)
all <- JoinLayers(all)
DefaultAssay(all) <- "RNA"

# Make combined identity: celltype_condition
all$celltype <- Idents(all)
all$celltype.condition <- paste(all$celltype, all$group, sep = "_")
Idents(all) <- "celltype.condition"

out_dir <- "results/tables/DEG"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_deg <- function(celltypes, cond1, cond2, suffix) {
  for (ct in celltypes) {
    ident1 <- paste0(ct, "_", cond1)
    ident2 <- paste0(ct, "_", cond2)

    try({
      de <- FindMarkers(
        all,
        ident.1 = ident1,
        ident.2 = ident2,
        test.use = "MAST",
        min.pct = 0.1,
        logfc.threshold = 0.25,
        verbose = FALSE
      )
      write.csv(de, file = file.path(out_dir, paste0(ct, "_", suffix, ".csv")))
    }, silent = TRUE)
  }
}

celltypes <- levels(factor(all$celltype))

# Comparisons (from the original notes)
run_deg(celltypes, "WP", "W",  "WPvW")
run_deg(celltypes, "K",  "W",  "KvW")
run_deg(celltypes, "KP", "K",  "KPvK")

# Example single comparison output retained in original notes
# (kept here as a sanity-check example, but redundant with the loop above if "LECs" exists)
if (all("LECs_WP" %in% levels(Idents(all)), "LECs_W" %in% levels(Idents(all)))) {
  de_lec <- FindMarkers(
    all,
    ident.1 = "LECs_WP",
    ident.2 = "LECs_W",
    test.use = "MAST",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  write.csv(de_lec, file = file.path(out_dir, "LECs_WPvW.csv"))
}

save_session_info()
message("Done: DEGs written to results/tables/DEG/")
