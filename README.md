# Mouse meningeal dura scRNA-seq

This repository contains a clean, GitHub-ready analysis pipeline for mouse meningeal dura scRNA-seq data using **Seurat v5**.

## Experimental groups
- **W**: Wild type  
- **WP**: Wild type + PLX5622  
- **K**: FIRE-KO  
- **KP**: FIRE-KO + PLX5622  

## Workflow overview
1. Read 10x matrices for each sample  
2. Ambient RNA correction with `decontX`  
3. Downsample to equal cell numbers per sample (16,000)  
4. QC metrics + filtering  
5. Merge samples  
6. Normalize, HVG, scale, PCA  
7. Harmony integration (batch by `orig.ident`)  
8. UMAP, neighbors, clustering  
9. Marker genes and manual cell-type annotation  
10. Differential expression per cell type (MAST)

## How to run
Run scripts in order from the repo root:

```bash
Rscript scripts/01_preprocessing_decontx.R
Rscript scripts/02_qc_filtering.R
Rscript scripts/03_integration_harmony.R
Rscript scripts/04_markers_and_annotation.R
Rscript scripts/05_differential_expression.R
```

Outputs are written to `results/` (figures, tables, objects).

## Reproducibility
Recommended: use **renv** to lock package versions.

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

A `session_info.txt` will be written by the pipeline.

## Notes
- Raw data are intentionally excluded; see `data/raw/README.md` for the expected layout.
- The pipeline uses relative paths to be portable.
