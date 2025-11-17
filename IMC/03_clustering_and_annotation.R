################################################################################
# IMC Analysis Pipeline
# Step 03: Clustering and Cell Type Annotation
# 
# Description:
#   - Manual cell type annotation based on marker expression
#   - Define cell type-specific color palettes
#   - Generate ComplexHeatmap for annotation verification
#   - Create barplots showing cell type proportions
#   - Export line graphs of marker expression per cell type
#
# Input:
#   - ./data/spe_batch_corrected.rds (from Step 02)
#
# Output:
#   - ./data/spe_annotated.rds
#   - Annotation heatmaps and proportion plots
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(SingleCellExperiment)
library(SpatialExperiment)
library(dittoSeq)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(scales)

# Set parameters
ANNOTATION_PARAMS <- list(
    seed = 900223
)

# Create output directory
if (!dir.exists("./figure/03_annotation")) {
    dir.create("./figure/03_annotation", recursive = TRUE)
}

################################################################################
# 1. Load Batch-Corrected Data
################################################################################

spe <- readRDS("./data/spe_batch_corrected.rds")

################################################################################
# 2. Manual Cell Type Annotation
################################################################################

# IMPORTANT: Annotation is based on marker expression patterns from heatmaps
# Adjust cluster-to-celltype mapping based on your specific markers
# This example is for pancreatic cancer IMC data

cluster_celltype <- recode(
    spe$som_clusters_harmonyCorrected,
    "1" = "Epithelial cell1(PanCK+)",
    "2" = "iCAF(SMA-PDGFRA+)",
    "3" = "DC",
    "4" = "myCAF(SMA+PDGFRA-COL1A1+)",
    "5" = "myCAF(SMA+PDGFRA-COL1A1+)",
    "6" = "Proliferating Epithelial cell2(Ki67+EpCAMloVimentinhi)",
    "7" = "CD4 T cell",
    "8" = "myCAF(SMA+PDGFRA-COL1A1+)",
    "9" = "Proliferating Epithelial cell2(Ki67+EpCAMloVimentinhi)",
    "10" = "mixed cell",
    "11" = "CD4 T cell",
    "12" = "Junk",
    "13" = "CTL",
    "14" = "Myeloid cell",
    "15" = "Monocyte",
    "16" = "Monocyte",
    "17" = "CTL",
    "18" = "Myeloid cell",
    "19" = "Neutrophil",
    "20" = "Epithelial cell3(EpCAMloVimentinhi)",
    "21" = "Myeloid cell",
    "22" = "Myeloid cell",
    "23" = "Myeloid cell",
    "24" = "Epithelial cell4(EpCAMhi)",
    "25" = "Myeloid cell",
    "26" = "iCAF(SMA-PDGFRA+)",
    "27" = "Myeloid cell",
    "28" = "Myeloid cell",
    "29" = "Epithelial cell1(PanCK+)",
    "30" = "Apoptotic cell",
    "31" = "Myeloid cell",
    "32" = "Myeloid cell",
    "33" = "Myeloid cell",
    "34" = "Myeloid cell",
    "35" = "Myeloid cell"
)

# Create outline annotation (same as cluster_celltype for segmentation plots)
outline_celltype <- cluster_celltype

# Assign to SPE object
spe$cluster_celltype <- cluster_celltype
spe$outline_celltype <- outline_celltype

################################################################################
# 3. Define Cell Type Color Palette
################################################################################

# IMPORTANT: Define colors for each unique cell type
# This ensures consistent colors across all visualizations
cluster_celltype_colors <- setNames(
    c('#FFD700',  # Epithelial cell1(PanCK+)
      '#DC143C',  # iCAF(SMA-PDGFRA+)
      '#FF00FF',  # DC
      '#8B0000',  # myCAF(SMA+PDGFRA-COL1A1+)
      '#FF6347',  # Proliferating Epithelial cell2(Ki67+EpCAMloVimentinhi)
      '#00BFFF',  # CD4 T cell
      '#228B22',  # mixed cell
      '#808080',  # Junk
      '#0000CD',  # CTL
      '#BA55D3',  # Myeloid cell
      '#483D8B',  # Monocyte
      '#8B008B',  # Neutrophil
      '#FF4500',  # Epithelial cell3(EpCAMloVimentinhi)
      '#FF8C00',  # Epithelial cell4(EpCAMhi)
      '#2F4F4F'), # Apoptotic cell
    c("Epithelial cell1(PanCK+)",
      "iCAF(SMA-PDGFRA+)",
      "DC",
      "myCAF(SMA+PDGFRA-COL1A1+)",
      "Proliferating Epithelial cell2(Ki67+EpCAMloVimentinhi)",
      "CD4 T cell",
      "mixed cell",
      "Junk",
      "CTL",
      "Myeloid cell",
      "Monocyte",
      "Neutrophil",
      "Epithelial cell3(EpCAMloVimentinhi)",
      "Epithelial cell4(EpCAMhi)",
      "Apoptotic cell")
)

outline_celltype_colors <- setNames(
    rep("black", length(unique(spe$cluster_celltype))),
    unique(spe$cluster_celltype)
)

# Store in metadata
metadata(spe)$color_vectors$cluster_celltype <- cluster_celltype_colors
metadata(spe)$color_vectors$outline_celltype <- outline_celltype_colors

################################################################################
# 4. Annotated UMAP Visualization
################################################################################

p_umap_annotated <- dittoDimPlot(
    spe, 
    var = "cluster_celltype",
    reduction.use = "UMAP_harmonyCorrected",
    size = 0.7,
    do.label = TRUE
) +
    scale_color_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme_classic(base_size = 14) +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 10, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_blank(),
        panel.grid = element_blank()
    ) +
    guides(color = guide_legend(
        override.aes = list(size = 4, stroke = 0.5)
    ))

ggsave(
    filename = "./figure/03_annotation/01_UMAP_annotated.tiff",
    plot = p_umap_annotated,
    width = 12,
    height = 8,
    dpi = 400,
    compression = "lzw"
)

################################################################################
# 5. ComplexHeatmap for Annotation Verification
################################################################################

# Calculate mean expression per cell type
celltype_mean <- aggregateAcrossCells(
    as(spe, "SingleCellExperiment"),
    ids = spe$cluster_celltype,
    statistics = "mean",
    use.assay.type = "exprs",
    subset.row = rownames(spe)[rowData(spe)$marker_class %in% c("type", "state")]
)

# Extract expression matrix
mat <- assay(celltype_mean, "exprs")

# Z-score transformation
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# Clip values for better visualization
clip <- 3
mat_z <- pmin(pmax(mat_z, -clip), clip)

# Hierarchical clustering
set.seed(ANNOTATION_PARAMS$seed)
hc <- hclust(dist(t(mat_z)), method = "ward.D2")
k_groups <- 5
grp <- cutree(hc, k = k_groups)
col_ord <- hc$order
col_split <- factor(grp[colnames(mat_z)[col_ord]], levels = sort(unique(grp)))

# Color scale
z_col <- colorRamp2(c(-clip, 0, clip), c("blue", "white", "red"))

# Cell type colors
celltype_cols <- metadata(spe)$color_vectors$cluster_celltype[colnames(mat_z)[col_ord]]

# Get cell counts
ncells_vec <- colData(celltype_mean)$ncells
ncells_vec <- ncells_vec[match(colnames(mat_z)[col_ord], 
                                rownames(colData(celltype_mean)))]

# Top annotation
ha_top <- HeatmapAnnotation(
    `Cell type` = colnames(mat_z)[col_ord],
    `#cells` = anno_barplot(
        ncells_vec,
        border = FALSE,
        height = unit(18, "mm"),
        gp = gpar(fill = "grey50")
    ),
    col = list(`Cell type` = celltype_cols),
    annotation_name_side = "left",
    simple_anno_size = unit(3, "mm")
)

# Block annotation
ha_block <- HeatmapAnnotation(
    block = anno_block(
        gp = gpar(fill = brewer.pal(k_groups, "Set3")),
        labels = paste0("G", levels(col_split)),
        labels_gp = gpar(fontsize = 10, fontface = "bold")
    )
)

# Row split by marker class
if ("marker_class" %in% colnames(rowData(celltype_mean))) {
    row_split <- factor(
        rowData(celltype_mean)$marker_class,
        levels = c("type", "state", "other")
    )
} else {
    row_split <- NULL
}

# Create heatmap
ht <- Heatmap(
    mat_z[, col_ord],
    name = "Z-score",
    col = z_col,
    
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 9),
    show_column_names = FALSE,
    
    rect_gp = gpar(col = "white", lwd = 0.5),
    border = gpar(col = "black", lwd = 1),
    use_raster = FALSE,
    
    column_split = col_split,
    cluster_columns = FALSE,
    column_gap = unit(1.2, "mm"),
    
    top_annotation = c(ha_block, ha_top),
    
    row_split = row_split,
    cluster_rows = TRUE,
    row_title_rot = 0,
    row_gap = unit(1, "mm"),
    
    heatmap_legend_param = list(
        at = c(-clip, 0, clip),
        title = "Z-score",
        legend_width = unit(4, "cm"),
        direction = "horizontal"
    )
)

tiff(
    filename = "./figure/03_annotation/02_annotation_heatmap.tiff",
    width = 14,
    height = 10,
    units = "in",
    res = 600,
    compression = "lzw"
)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()

################################################################################
# 6. Cell Type Proportion Barplots
################################################################################

# By sample (ROI level)
p_bar_sample <- dittoBarPlot(
    spe,
    var = "cluster_celltype",
    group.by = "sample_id",
    legend.show = TRUE
) +
    scale_fill_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.title.x = element_blank(),
        legend.title = element_blank()
    )

ggsave(
    filename = "./figure/03_annotation/03_barplot_by_sample.tiff",
    plot = p_bar_sample,
    width = 12,
    height = 7,
    dpi = 400,
    compression = "lzw"
)

# By group
p_bar_group <- dittoBarPlot(
    spe,
    var = "cluster_celltype",
    group.by = "group",
    legend.show = TRUE
) +
    scale_fill_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
    )

ggsave(
    filename = "./figure/03_annotation/04_barplot_by_group.tiff",
    plot = p_bar_group,
    width = 6,
    height = 7,
    dpi = 400,
    compression = "lzw"
)

# Export raw proportion data
bar_data <- dittoBarPlot(
    spe,
    var = "cluster_celltype",
    group.by = "group",
    data.out = TRUE
)

write.csv(
    bar_data$data,
    file = "./figure/03_annotation/proportion_by_group.csv",
    row.names = FALSE
)

################################################################################
# 7. Line Graphs: Marker Expression per Cell Type
################################################################################

# NOTE: This creates line graphs showing marker expression patterns
# Useful for manuscript supplementary figures

# Get expression data
expr_data <- assay(celltype_mean, "exprs")
expr_data_t <- t(expr_data)

# Convert to long format
expr_df <- data.frame(expr_data_t)
expr_df$cluster_celltype <- rownames(expr_data_t)

expr_long <- melt(
    expr_df,
    id.vars = "cluster_celltype",
    variable.name = "Marker",
    value.name = "Expression"
)

# Define marker order for visualization
desired_markers <- c(
    "Vimentin", "SMA", "CD31", "EpCAM", "Ki67",
    "CD45", "CD68", "CD11b", "MHC-II", "CD206",
    "CD11c", "CD3", "CD4", "CD8", "B220"
)

available_markers <- intersect(desired_markers, unique(expr_long$Marker))

if (length(available_markers) > 0) {
    expr_long_filtered <- subset(expr_long, Marker %in% available_markers)
    expr_long_filtered$Marker <- factor(
        expr_long_filtered$Marker,
        levels = available_markers
    )
    
    # Normalize expression (0-1 scaling per marker)
    expr_long_normalized <- expr_long_filtered %>%
        group_by(Marker) %>%
        mutate(
            Normalized_Expression = (Expression - min(Expression)) / 
                                  (max(Expression) - min(Expression))
        ) %>%
        ungroup()
    
    # Get colors
    color_vector <- metadata(spe)$color_vectors$cluster_celltype
    
    # Create output directory
    dir.create("./figure/03_annotation/line_graphs", 
               showWarnings = FALSE, 
               recursive = TRUE)
    
    # Generate line graph for each cell type
    unique_celltypes <- unique(expr_long_normalized$cluster_celltype)
    
    for (celltype in unique_celltypes) {
        data_subset <- subset(
            expr_long_normalized,
            cluster_celltype == celltype
        )
        
        p <- ggplot(
            data_subset,
            aes(x = Marker, y = Normalized_Expression, group = 1)
        ) +
            geom_line(color = color_vector[celltype], size = 1.5) +
            geom_point(color = color_vector[celltype], size = 4, shape = 15) +
            theme_classic() +
            theme(
                axis.text.x = element_blank(),
                plot.title = element_text(hjust = 0.5),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.line = element_line(size = 1)
            ) +
            scale_y_continuous(
                limits = c(0, 1),
                expand = c(0, 0),
                breaks = c(0, 1)
            ) +
            coord_cartesian(ylim = c(0, 1.02))
        
        # Safe filename
        safe_celltype <- gsub("[^[:alnum:] ]", "", celltype)
        safe_celltype <- gsub(" ", "_", safe_celltype)
        
        ggsave(
            filename = sprintf(
                "./figure/03_annotation/line_graphs/%s_normalized.png",
                safe_celltype
            ),
            plot = p,
            width = 5,
            height = 4,
            dpi = 300
        )
    }
}

################################################################################
# 8. Save Annotated Object
################################################################################

saveRDS(spe, "./data/spe_annotated.rds")

# Clean up
rm(celltype_mean, mat, mat_z, expr_long, ht)
gc()

# NOTE: Next step is 04_quality_control.R
