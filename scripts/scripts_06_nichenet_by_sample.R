#!/usr/bin/env Rscript

# 06_nichenet_by_sample.R
# NicheNet cell–cell interaction analysis for one sample (K/KP/W/WP).
#
# Inputs:
#   - results/objects/04_harmony_umap_annotated.rds  (Seurat object with Idents = cell types and meta $group)
#   - data/metadata/mouse_ref_nichnet.RData         (NicheNet reference objects: lr_network, ligand_target_matrix, weighted_networks, ...)
#   - data/metadata/LEC_DEGs_1313.txt               (receiver gene set; one gene symbol per line)
#
# Outputs (written to results/):
#   - results/figures/nichenet/K/  (PDFs including circos plots used in Fig. 2h and Ext. Fig. 6d)
#   - results/tables/nichenet/K/   (tables used for plotting / manuscript source data)
#
# Notes:
# - This script reproduces the sample-specific analyses, but with repo-relative paths and clean outputs.
# - For KP/W/WP, copy this script and change SAMPLE at the top, or parameterize via commandArgs().

suppressPackageStartupMessages({
  library(nichenetr)
  library(tidyverse)
  library(circlize)
  library(Seurat)
})

source("scripts/_helpers.R")
set_analysis_seed(1234)

# ----------------------------
# Configuration
# ----------------------------
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Accept: --sample K|KP|W|WP  (also allow: --sample=K)
sample_arg <- NULL
for (a in args) {
  if (grepl("^--sample=", a)) {
    sample_arg <- sub("^--sample=", "", a)
  }
}
if (is.null(sample_arg)) {
  idx <- match("--sample", args)
  if (!is.na(idx) && length(args) >= idx + 1) {
    sample_arg <- args[idx + 1]
  }
}

if (is.null(sample_arg) || !(sample_arg %in% c("K","KP","W","WP"))) {
  stop("Usage: Rscript scripts/06_nichenet_by_sample.R --sample K|KP|W|WP")
}

SAMPLE <- sample_arg
receiver <- "LECs"

sender_celltypes <- c(
  "Mast cells",
  "T, B, NK cells",
  "Plasma cells",
  "Neutrophils",
  "CXCR2-high Neutrophils",
  "Monocytes",
  "CCR2-high Monocytes",
  "Dendritic cells",
  "Macrophages",
  "Lyve1-high Macs"
)

# Output directories
fig_dir <- file.path("results/figures/nichenet", SAMPLE)
tbl_dir <- file.path("results/tables/nichenet", SAMPLE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load inputs
# ----------------------------
all <- readRDS("results/objects/04_harmony_umap_annotated.rds")
stopifnot("group" %in% colnames(all@meta.data))

# Subset to one sample and relevant cell types (receiver + senders)
seuratObj <- subset(
  all,
  subset = group == SAMPLE,
  idents = c(sender_celltypes, receiver)
)

seuratObj@meta.data$celltype <- Idents(seuratObj)

# NicheNet reference objects (provided separately; not tracked in git)
load("data/metadata/mouse_ref_nichnet.RData")

# Receiver gene set of interest (e.g., LEC DEGs used for the manuscript panel)
genes <- read.delim("data/metadata/LEC_DEGs_1313.txt", header = FALSE)
geneset_oi <- as.character(genes$V1)

# Ensure Seurat v5 layers are joined for downstream functions
seuratObj <- JoinLayers(seuratObj)

# ----------------------------
# NicheNet analysis
# ----------------------------

### Because we combined the expressed genes of each sender cell type, in this example, we will perform one NicheNet analysis by pooling all ligands from all cell types together. Later on during the interpretation of the output, we will check which sender cell type expresses which ligand.
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

####Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(40, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

####To see which cell population expresses which of these top-ranked ligands, you can run the following:
seuratObj_1 <- subset(seuratObj, idents="LECs", invert=T)
pdf(file = file.path(fig_dir, paste0("DotPlot_best_upstream_ligands_by_celltype_", SAMPLE, ".pdf")), width = 12,height = 5)
DotPlot(seuratObj_1,assay="RNA", features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
dev.off()

####Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
pdf(file = file.path(fig_dir, paste0("Heatmap_ligand_target_network_", SAMPLE, ".pdf")), width = 40,height = 5)
p_ligand_target_network
dev.off()

#### Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
pdf(file = file.path(fig_dir, paste0("Heatmap_ligand_receptor_network_", SAMPLE, ".pdf")))
p_ligand_receptor_network
dev.off()

#### circos plot
#Calculate average ligand expression in sender cells
avg_expression_ligands = AverageExpression(seuratObj, features = best_upstream_ligands, assay="RNA")

sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
  }) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)
#### note: CXCR2-high Neutrophils doesn't have ligands and are removed at this step above

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)


Neutrophils_specific_ligands = sender_ligand_assignment$"Neutrophils" %>% names() %>% setdiff(general_ligands)
Monocytes_specific_ligands = sender_ligand_assignment$"Monocytes" %>% names() %>% setdiff(general_ligands)
Mast_cells_specific_ligands = sender_ligand_assignment$"Mast cells" %>% names() %>% setdiff(general_ligands)
Dendritic_cells_specific_ligands = sender_ligand_assignment$"Dendritic cells" %>% names() %>% setdiff(general_ligands)
Macrophages_specific_ligands = sender_ligand_assignment$"Macrophages" %>% names() %>% setdiff(general_ligands)
T_B_NK_cells_specific_ligands = sender_ligand_assignment$"T, B, NK cells" %>% names() %>% setdiff(general_ligands)
CCR2_high_Monocytes_specific_ligands = sender_ligand_assignment$"CCR2-high Monocytes" %>% names() %>% setdiff(general_ligands)
Lyve1_high_Macs_specific_ligands = sender_ligand_assignment$"Lyve1-high Macs" %>% names() %>% setdiff(general_ligands)
Plasma_cells_specific_ligands = sender_ligand_assignment$"Plasma cells" %>% names() %>% setdiff(general_ligands)

### Monocytes_specific_ligands, CCR2_high_Monocytes_specific_ligands have 0 ligands

ligand_type_indication_df = tibble(
  ligand_type = c(
		  rep("Mast_cells_specific", times = Mast_cells_specific_ligands %>% length()),
		  rep("T_B_NK_cells_specific", times = T_B_NK_cells_specific_ligands %>% length()),
		  rep("Plasma_cells_specific", times = Plasma_cells_specific_ligands %>% length()),
		  rep("Neutrophils_specific", times = Neutrophils_specific_ligands %>% length()),
		  rep("Monocytes_specific", times = Monocytes_specific_ligands %>% length()),
		  rep("CCR2_high_Monocytes_specific", times = CCR2_high_Monocytes_specific_ligands %>% length()),
		  rep("Dendritic_cells_specific", times = Dendritic_cells_specific_ligands %>% length()),
		  rep("Macrophages_specific", times = Macrophages_specific_ligands %>% length()),
		  rep("Lyve1_high_Macs_specific", times = Lyve1_high_Macs_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(Mast_cells_specific_ligands,T_B_NK_cells_specific_ligands,Plasma_cells_specific_ligands,Neutrophils_specific_ligands,Monocytes_specific_ligands,CCR2_high_Monocytes_specific_ligands,Dendritic_cells_specific_ligands,Macrophages_specific_ligands,Lyve1_high_Macs_specific_ligands, general_ligands))

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
### export the lignda target strength table
#ligand_target_strength <- as.data.frame(active_ligand_target_links_df)


active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "DEGs") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
### Remove rows with NA
circos_links = active_ligand_target_links_df %>% filter(!is.na(weight))
table(circos_links$ligand_type)

### export the lignda target strength table
ligand_target_strength <- as.data.frame(circos_links)
write.table(ligand_target_strength, file.path(tbl_dir, paste0(SAMPLE, "_ligand_target_strength_immune_sender.csv")), sep=",")


#####Prepare plot
grid_col_ligand =c("General" = "lightsteelblue4",
		   "Dendritic_cells_specific"="goldenrod4",
		   "Lyve1_high_Macs_specific"="cadetblue",
		   "Macrophages_specific"="firebrick1",
		   "Mast_cells_specific"="powderblue",
		   "Neutrophils_specific"="yellow",
		   "Plasma_cells_specific"="ivory4",
		   "T_B_NK_cells_specific"="lightpink"       
		   )
grid_col_target =c(
            "DEGs" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
#############
######################################################################
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency


target_order = circos_links$target %>% unique()
ligand_order = c(Mast_cells_specific_ligands,T_B_NK_cells_specific_ligands,Plasma_cells_specific_ligands,Neutrophils_specific_ligands,Dendritic_cells_specific_ligands,Macrophages_specific_ligands,Lyve1_high_Macs_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)


width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mast_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "T_B_NK_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Plasma_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophils_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Dendritic_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophages_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Lyve1_high_Macs_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "DEGs") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
  )

pdf(file = file.path(fig_dir, paste0("ExtFig6d_circos_ligand_target_", SAMPLE, ".pdf")), width = 16, height = 14)
#circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow",annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #
circos.clear()
dev.off()

######Visualize ligand-receptor interactions of the prioritized ligands in a circos plot
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "DEGs") %>% inner_join(ligand_type_indication_df)



grid_col_ligand =c("General" = "lightsteelblue4",
			   "Dendritic_cells_specific"="goldenrod4",
                   "Lyve1_high_Macs_specific"="cadetblue",
                   "Macrophages_specific"="firebrick1",
                   "Mast_cells_specific"="powderblue",
                   "Neutrophils_specific"="yellow",
			"Plasma_cells_specific"="ivory4",
                   "T_B_NK_cells_specific"="lightpink")
grid_col_receptor =c(
            "DEGs" = "darkred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

### export the lignda target strength table
ligand_receptor_strength <- as.data.frame(circos_links)
write.table(ligand_receptor_strength, file.path(tbl_dir, paste0(SAMPLE, "_ligand_receptor_strength_immune_sender.csv")), sep=",")



ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score

transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency

receptor_order = circos_links$receptor %>% unique()
ligand_order = c(Mast_cells_specific_ligands,T_B_NK_cells_specific_ligands,Plasma_cells_specific_ligands,Neutrophils_specific_ligands,Dendritic_cells_specific_ligands,Macrophages_specific_ligands,Lyve1_high_Macs_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5
gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mast_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "T_B_NK_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Plasma_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophils_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Dendritic_cells_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophages_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Lyve1_high_Macs_specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "DEGs") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
  )

pdf(file = file.path(fig_dir, paste0("Fig2h_circos_ligand_receptor_", SAMPLE, ".pdf")), width = 12, height = 10)
#circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()
dev.off()


# ----------------------------
# Final bookkeeping
# ----------------------------
save_session_info()
message("Done. Figures: ", fig_dir, " | Tables: ", tbl_dir)
