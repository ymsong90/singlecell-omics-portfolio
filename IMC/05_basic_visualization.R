################################################################################
# IMC Analysis Pipeline
# Step 05: Basic Visualization
# 
# Description:
#   - Publication-quality UMAP/t-SNE plots
#   - Cell type proportion barplots (by sample and group)
#   - ComplexHeatmap with annotations (cell counts, group proportions, spatial features)
#   - Marker expression heatmaps (z-score scaled)
#   - Spatial cell segmentation images
#   - plotPixels for channel visualization
#   - Export raw proportion data
#
# Input:
#   - ./data/spe_filtered.rds (from Step 04)
#   - ./data/images.rds
#   - ./data/masks.rds
#
# Output:
#   - Multiple publication-quality figures in ./figure/05_visualization/
#   - Raw proportion data CSV files
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(SpatialExperiment)
library(SingleCellExperiment)
library(dittoSeq)
library(ComplexHeatmap)
library(circlize)
library(cytomapper)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(scales)
library(scuttle)

# Set parameters
VIZ_PARAMS <- list(
    umap_point_size       = 0.7,
    n_samples_spatial     = 3,
    heatmap_zscore_clip   = 3,
    seed                  = 900223
)

# Create output directory
if (!dir.exists("./figure/05_visualization")) {
    dir.create("./figure/05_visualization", recursive = TRUE)
}

################################################################################
# 1. Load Filtered Data
################################################################################

spe <- readRDS("./data/spe_filtered.rds")

# Load images and masks if spatial visualization is needed
if (file.exists("./data/images.rds") && file.exists("./data/masks.rds")) {
    images <- readRDS("./data/images.rds")
    masks <- readRDS("./data/masks.rds")
    spatial_available <- TRUE
} else {
    spatial_available <- FALSE
    warning("Images/masks not found. Spatial visualization will be skipped.")
}

################################################################################
# 2. Publication-Quality UMAP
################################################################################

# Determine which UMAP to use
umap_name <- if ("UMAP_harmonyCorrected" %in% names(reducedDims(spe))) {
    "UMAP_harmonyCorrected"
} else {
    "UMAP"
}

p_umap <- dittoDimPlot(
    spe,
    var = "cluster_celltype",
    reduction.use = umap_name,
    size = VIZ_PARAMS$umap_point_size,
    do.label = FALSE
) +
    scale_color_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme_classic(base_size = 14) +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.6, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        legend.margin = margin(l = 20),
        
        axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        axis.text = element_text(size = 12, color = "black"),
        
        plot.title = element_blank(),
        plot.margin = margin(20, 20, 20, 20),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank()
    ) +
    guides(color = guide_legend(
        override.aes = list(size = 4, stroke = 0.5)
    ))

ggsave(
    filename = "./figure/05_visualization/01_UMAP_publication.tiff",
    plot = p_umap,
    width = 12,
    height = 8,
    dpi = 400,
    compression = "lzw"
)

# UMAP split by group
p_umap_split <- dittoDimPlot(
    spe,
    var = "cluster_celltype",
    reduction.use = umap_name,
    split.by = "group",
    size = 0.5,
    do.label = FALSE
) +
    scale_color_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme_bw() +
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid = element_blank()
    )

ggsave(
    filename = "./figure/05_visualization/02_UMAP_split_by_group.tiff",
    plot = p_umap_split,
    width = 16,
    height = 6,
    dpi = 400,
    compression = "lzw"
)

################################################################################
# 3. Cell Type Proportion Barplots
################################################################################

# By sample (ROI level)
p_bar_sample <- dittoBarPlot(
    spe,
    var = "cluster_celltype",
    group.by = "sample_id",
    legend.show = TRUE
) +
    scale_fill_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
    )

ggsave(
    filename = "./figure/05_visualization/03_proportion_by_sample.tiff",
    plot = p_bar_sample,
    width = 14,
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
    theme_classic(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.ticks = element_blank(),
        panel.grid = element_blank()
    )

ggsave(
    filename = "./figure/05_visualization/04_proportion_by_group.tiff",
    plot = p_bar_group,
    width = 6,
    height = 7,
    dpi = 400,
    compression = "lzw"
)

# Export raw proportion data
dir.create("./figure/05_visualization/csv_data", 
           showWarnings = FALSE, 
           recursive = TRUE)

# By sample
bar_data_sample <- dittoBarPlot(
    spe,
    var = "cluster_celltype",
    group.by = "sample_id",
    data.out = TRUE
)

write.csv(
    bar_data_sample$data,
    file = "./figure/05_visualization/csv_data/proportion_by_sample.csv",
    row.names = FALSE
)

# By group
bar_data_group <- dittoBarPlot(
    spe,
    var = "cluster_celltype",
    group.by = "group",
    data.out = TRUE
)

write.csv(
    bar_data_group$data,
    file = "./figure/05_visualization/csv_data/proportion_by_group.csv",
    row.names = FALSE
)

################################################################################
# 4. ComplexHeatmap with Multiple Annotations
################################################################################

# Calculate mean expression per cell type
celltype_mean <- aggregateAcrossCells(
    as(spe, "SingleCellExperiment"),
    ids = spe$cluster_celltype,
    statistics = "mean",
    use.assay.type = "exprs",
    subset.row = rownames(spe)[rowData(spe)$marker_class %in% c("type", "state")]
)

# Extract and transform matrix
mat <- assay(celltype_mean, "exprs")
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# Clip extreme values
clip <- VIZ_PARAMS$heatmap_zscore_clip
mat_z <- pmin(pmax(mat_z, -clip), clip)

# Hierarchical clustering
set.seed(VIZ_PARAMS$seed)
hc <- hclust(dist(t(mat_z)), method = "ward.D2")
k_groups <- 5
grp <- cutree(hc, k = k_groups)
col_ord <- hc$order
col_split <- factor(grp[colnames(mat_z)[col_ord]], levels = sort(unique(grp)))

# Color scales
z_col <- colorRamp2(c(-clip, 0, clip), c("blue", "white", "red"))

# Cell counts
ncells_vec <- colData(celltype_mean)$ncells
ncells_vec <- ncells_vec[match(colnames(mat_z)[col_ord],
                                rownames(colData(celltype_mean)))]

# Group proportions per celltype
group_props <- prop.table(table(spe$cluster_celltype, spe$group), margin = 1)
group_props <- group_props[colnames(mat_z)[col_ord], ]

# Top annotation
ha_top <- HeatmapAnnotation(
    `Cell type` = colnames(mat_z)[col_ord],
    `#cells` = anno_barplot(
        ncells_vec,
        border = FALSE,
        height = unit(18, "mm"),
        gp = gpar(fill = "grey50")
    ),
    `Group proportions` = anno_barplot(
        group_props,
        border = FALSE,
        height = unit(18, "mm"),
        gp = gpar(fill = metadata(spe)$color_vectors$group)
    ),
    col = list(
        `Cell type` = metadata(spe)$color_vectors$cluster_celltype[colnames(mat_z)[col_ord]]
    ),
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
row_split <- factor(
    rowData(celltype_mean)$marker_class,
    levels = c("type", "state", "other")
)

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
    filename = "./figure/05_visualization/05_ComplexHeatmap_annotated.tiff",
    width = 16,
    height = 12,
    units = "in",
    res = 600,
    compression = "lzw"
)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()

################################################################################
# 5. Marker Expression on UMAP
################################################################################

# Select key markers for visualization
key_markers <- c("CD3", "CD4", "CD8", "CD68", "Ki67", "EpCAM")
available_key_markers <- intersect(key_markers, rownames(spe))

if (length(available_key_markers) > 0) {
    p_markers <- lapply(available_key_markers, function(marker) {
        dittoDimPlot(
            spe,
            var = marker,
            reduction.use = umap_name,
            assay = "exprs",
            size = 0.3
        ) +
            scale_color_viridis(name = marker, option = "plasma") +
            theme_classic(base_size = 10) +
            theme(
                legend.position = "right",
                plot.title = element_text(face = "bold", size = 11)
            ) +
            ggtitle(marker)
    })
    
    ggsave(
        filename = "./figure/05_visualization/06_marker_expression_UMAP.tiff",
        plot = wrap_plots(p_markers, ncol = 3),
        width = 15,
        height = 10,
        dpi = 400,
        compression = "lzw"
    )
}

################################################################################
# 6. Spatial Visualization
################################################################################

if (spatial_available) {
    
    # Select random samples for visualization
    set.seed(VIZ_PARAMS$seed)
    sample_ids <- unique(spe$sample_id)
    cur_id <- sample(sample_ids, min(VIZ_PARAMS$n_samples_spatial, length(sample_ids)))
    
    images_subset <- images[names(images) %in% cur_id]
    masks_subset <- masks[names(masks) %in% cur_id]
    
    # Ensure order matches
    masks_subset <- masks_subset[names(images_subset)]
    
    # Normalize images
    images_subset <- normalize(images_subset)
    images_subset <- normalize(images_subset, inputRange = c(0, 0.2))
    
    # Plot cell segmentation with cell type colors
    tiff(
        filename = "./figure/05_visualization/07_cell_segmentation.tiff",
        width = 12,
        height = 4 * length(cur_id),
        units = "in",
        res = 400,
        compression = "lzw"
    )
    plotCells(
        masks_subset,
        object = spe,
        cell_id = "ObjectNumber",
        img_id = "sample_id",
        colour_by = "cluster_celltype",
        outline_by = "outline_celltype",
        colour = list(
            cluster_celltype = metadata(spe)$color_vectors$cluster_celltype,
            outline_celltype = metadata(spe)$color_vectors$outline_celltype
        ),
        display = "single",
        legend = list(
            colour_by.legend.cex = 0.7,
            outline_by.legend.cex = 0.5
        )
    )
    dev.off()
    
    # plotPixels for key markers
    plot_markers <- c("CD68", "CD4", "EpCAM")
    available_plot_markers <- intersect(plot_markers, rownames(spe))
    
    if (length(available_plot_markers) >= 2) {
        tiff(
            filename = "./figure/05_visualization/08_spatial_pixels.tiff",
            width = 12,
            height = 4 * length(cur_id),
            units = "in",
            res = 400,
            compression = "lzw"
        )
        plotPixels(
            images_subset,
            colour_by = available_plot_markers[1:min(3, length(available_plot_markers))],
            colour = list(
                CD68 = c("black", "yellow"),
                CD4 = c("black", "green"),
                EpCAM = c("black", "blue")
            ),
            img_id = "sample_id",
            missing_colour = "white",
            legend = list(
                colour_by.title.cex = 0.7,
                colour_by.labels.cex = 0.7
            )
        )
        dev.off()
    }
}

################################################################################
# 7. Cell Size and QC Metrics by Cell Type
################################################################################

# Cell size distribution
p_size <- colData(spe) %>%
    as.data.frame() %>%
    ggplot(aes(x = cluster_celltype, y = area, fill = cluster_celltype)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    scale_fill_manual(values = metadata(spe)$color_vectors$cluster_celltype) +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.x = element_blank(),
        legend.position = "none"
    ) +
    ylab("Cell Area (pixels)")

ggsave(
    filename = "./figure/05_visualization/09_cell_size_by_type.tiff",
    plot = p_size,
    width = 12,
    height = 6,
    dpi = 400,
    compression = "lzw"
)

# Cell counts per sample
cell_counts <- colData(spe) %>%
    as.data.frame() %>%
    group_by(sample_id, group, cluster_celltype) %>%
    summarize(n = n(), .groups = "drop")

write.csv(
    cell_counts,
    file = "./figure/05_visualization/csv_data/cell_counts_per_sample.csv",
    row.names = FALSE
)

################################################################################
# 8. Summary Statistics
################################################################################

# Generate summary table
summary_stats <- colData(spe) %>%
    as.data.frame() %>%
    group_by(cluster_celltype) %>%
    summarize(
        n_cells = n(),
        mean_area = mean(area, na.rm = TRUE),
        median_area = median(area, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(desc(n_cells))

write.csv(
    summary_stats,
    file = "./figure/05_visualization/csv_data/summary_statistics.csv",
    row.names = FALSE
)

################################################################################
# 9. Clean Up
################################################################################

rm(celltype_mean, mat, mat_z, ht, images_subset, masks_subset)
gc()

