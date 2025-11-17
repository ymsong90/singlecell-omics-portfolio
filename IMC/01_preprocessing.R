################################################################################
# IMC Analysis Pipeline
# Step 01: Preprocessing
# 
# Description:
#   - Load Steinbock IMC data
#   - Map sample metadata (group, block)
#   - Transform counts (asinh with cofactor 5)
#   - Define marker classes (type vs state)
#   - Initial dimensionality reduction (UMAP, t-SNE)
#   - Basic QC and visualization
#
# Input:
#   - ./data/steinbock/ directory (Steinbock output)
#   - ./data/steinbock/sample_metadata.csv
#
# Output:
#   - ./data/spe_preprocessed.rds
#   - ./data/images.rds
#   - ./data/masks.rds
#   - QC plots in ./figure/01_preprocessing/
#
# Author: YMS
# Date: 2024-11-17
################################################################################

# Load required packages
library(imcRtools)
library(SingleCellExperiment)
library(SpatialExperiment)
library(S4Vectors)
library(cytomapper)
library(scater)
library(dittoSeq)
library(scuttle)
library(viridis)
library(ggplot2)
library(stringr)
library(dplyr)
library(scales)
library(patchwork)

# Set parameters
PREPROCESSING_PARAMS <- list(
    transform_cofactor  = 5,           # asinh transformation cofactor
    min_cell_area      = 5,           # minimum cell area (pixels)
    n_cells_heatmap    = 2000,        # cells for QC heatmap
    n_samples_spatial  = 3,           # samples for spatial visualization
    seed               = 900223
)

# Create output directories
if (!dir.exists("./figure/01_preprocessing")) {
    dir.create("./figure/01_preprocessing", recursive = TRUE)
}

if (!dir.exists("./data")) {
    dir.create("./data", recursive = TRUE)
}

################################################################################
# 1. Load Steinbock Data
################################################################################

# NOTE: read_steinbock() loads all Steinbock outputs at once
# Including intensities, regionprops, neighbors, images, masks
spe <- read_steinbock("./data/steinbock/")

# Create unique cell identifiers: sample_id + ObjectNumber
colnames(spe) <- paste0(spe$sample_id, "_", spe$ObjectNumber)

################################################################################
# 2. Load and Map Metadata
################################################################################

meta <- read.csv("./data/steinbock/sample_metadata.csv", 
                 stringsAsFactors = FALSE)

# IMPORTANT: Metadata must contain sample_id, group, block columns
stopifnot(all(c("sample_id", "group", "block") %in% colnames(meta)))
stopifnot(!anyDuplicated(meta$sample_id))

# Map metadata to SPE object
idx <- match(spe$sample_id, meta$sample_id)

# WARNING: Stop if any sample_id is missing from metadata
if (anyNA(idx)) {
    missing_samples <- unique(spe$sample_id[is.na(idx)])
    stop("Missing sample_id in metadata: ", 
         paste(head(missing_samples, 10), collapse = ", "))
}

spe$group <- factor(meta$group[idx], levels = unique(meta$group))
spe$block <- as.integer(meta$block[idx])

################################################################################
# 3. Transform Counts
################################################################################

# NOTE: asinh transformation is standard for IMC data
# Formula: asinh(counts / cofactor)
# This stabilizes variance and approximates log transformation
assay(spe, "exprs") <- asinh(counts(spe) / PREPROCESSING_PARAMS$transform_cofactor)

################################################################################
# 4. Define Marker Classes
################################################################################

# Define markers to exclude from analysis
# NOTE: DNA channels are typically excluded from downstream analysis
exclude_markers <- c("DNA1", "DNA2", "CD45", "Vimentin", "MHC-II")

rowData(spe)$use_channel <- !grepl(
    paste(exclude_markers, collapse = "|"), 
    rownames(spe), 
    ignore.case = TRUE
)

################################################################################
# 5. Define Color Palettes
################################################################################

# Group colors (Wt/Ko/Panc)
group_colors <- c(Wt = "#4e79a7", Ko = "#e15759", Panc = "#59a14f")
group_levels <- levels(spe$group)

# Add colors for additional groups if needed
missing_groups <- setdiff(group_levels, names(group_colors))
if (length(missing_groups) > 0) {
    extra_colors <- setNames(
        hue_pal()(length(missing_groups)), 
        missing_groups
    )
    group_colors <- c(group_colors, extra_colors)
}
group_colors <- group_colors[group_levels]

# Block colors
block_levels <- sort(unique(spe$block))
block_colors <- setNames(
    hue_pal()(length(block_levels)), 
    as.character(block_levels)
)

# Store in metadata
metadata(spe)$color_vectors <- list(
    group = group_colors,
    block = block_colors
)

################################################################################
# 6. Load Images and Masks
################################################################################

images <- loadImages("./data/steinbock/img/")
masks <- loadImages("./data/steinbock/masks/", as.is = TRUE)

# Set channel names to match marker names
channelNames(images) <- rownames(spe)

# Add metadata to images
group_img <- meta$group[match(names(images), meta$sample_id)]
mcols(images) <- mcols(masks) <- DataFrame(
    sample_id = names(images),
    group = group_img
)

################################################################################
# 7. Dimensionality Reduction
################################################################################

# NOTE: Only use markers defined in use_channel
# This excludes DNA and other unwanted markers
set.seed(PREPROCESSING_PARAMS$seed)

spe <- runUMAP(
    spe, 
    subset_row = rowData(spe)$use_channel,
    exprs_values = "exprs"
)

spe <- runTSNE(
    spe,
    subset_row = rowData(spe)$use_channel,
    exprs_values = "exprs"
)

################################################################################
# 8. Basic Visualization
################################################################################

# UMAP and t-SNE by group
p1 <- dittoDimPlot(
    spe, 
    var = "group", 
    reduction.use = "UMAP", 
    size = 0.2
) +
    scale_color_manual(values = metadata(spe)$color_vectors$group) +
    ggtitle("UMAP by Group")

p2 <- dittoDimPlot(
    spe, 
    var = "group", 
    reduction.use = "TSNE", 
    size = 0.2
) +
    scale_color_manual(values = metadata(spe)$color_vectors$group) +
    ggtitle("t-SNE by Group")

ggsave(
    filename = "./figure/01_preprocessing/01_UMAP_TSNE_by_group.tiff",
    plot = p1 + p2,
    width = 14,
    height = 6,
    dpi = 300,
    compression = "lzw"
)

# Marker expression visualization
example_markers <- c("CD3", "CD4", "Ki67", "FOXP3")
available_markers <- intersect(example_markers, rownames(spe))

if (length(available_markers) > 0) {
    p_markers <- lapply(available_markers, function(marker) {
        dittoDimPlot(
            spe, 
            var = marker, 
            reduction.use = "UMAP",
            assay = "exprs", 
            size = 0.2
        ) +
            scale_color_viridis(name = marker)
    })
    
    ggsave(
        filename = "./figure/01_preprocessing/02_marker_expression_UMAP.tiff",
        plot = wrap_plots(p_markers, ncol = 2),
        width = 10,
        height = 10,
        dpi = 300,
        compression = "lzw"
    )
}

################################################################################
# 9. QC: Per-Cell Heatmap
################################################################################

# Sample random cells for heatmap (to reduce computation time)
set.seed(PREPROCESSING_PARAMS$seed)
cur_cells <- sample(
    seq_len(ncol(spe)), 
    min(PREPROCESSING_PARAMS$n_cells_heatmap, ncol(spe))
)

tiff(
    filename = "./figure/01_preprocessing/03_per_cell_heatmap.tiff",
    width = 10,
    height = 8,
    units = "in",
    res = 300,
    compression = "lzw"
)
dittoHeatmap(
    spe[, cur_cells],
    genes = rownames(spe)[rowData(spe)$use_channel],
    assay = "exprs",
    cluster_cols = TRUE,
    scale = "none",
    heatmap.colors = viridis(100),
    annot.by = "group",
    annotation_colors = list(
        group = metadata(spe)$color_vectors$group
    )
)
dev.off()

################################################################################
# 10. QC: Per-Image Mean Heatmap
################################################################################

# Calculate mean expression per sample
image_mean <- aggregateAcrossCells(
    spe,
    ids = spe$sample_id,
    statistics = "mean",
    use.assay.type = "counts"
)

# Transform aggregated counts
assay(image_mean, "exprs") <- asinh(
    counts(image_mean) / PREPROCESSING_PARAMS$transform_cofactor
)

tiff(
    filename = "./figure/01_preprocessing/04_per_image_mean_heatmap.tiff",
    width = 12,
    height = 9,
    units = "in",
    res = 300,
    compression = "lzw"
)
dittoHeatmap(
    image_mean,
    genes = rownames(spe)[rowData(spe)$use_channel],
    assay = "exprs",
    cluster_cols = TRUE,
    scale = "none",
    heatmap.colors = viridis(100),
    annot.by = c("group", "sample_id"),
    annotation_colors = list(
        group = metadata(spe)$color_vectors$group
    ),
    show_colnames = TRUE
)
dev.off()

################################################################################
# 11. QC: Image Coverage and Cell Size
################################################################################

# Calculate area covered by cells per image
p_coverage <- colData(spe) %>%
    as.data.frame() %>%
    group_by(sample_id) %>%
    summarize(
        cell_area = sum(area),
        no_pixels = mean(width_px) * mean(height_px)
    ) %>%
    mutate(covered_area = cell_area / no_pixels) %>%
    ggplot() +
    geom_point(aes(reorder(sample_id, covered_area), covered_area)) +
    theme_minimal(base_size = 12) +
    ylim(c(0, 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    ylab("% Covered Area") +
    xlab("")

# Cell size distribution per image
p_cellsize <- colData(spe) %>%
    as.data.frame() %>%
    ggplot() +
    geom_boxplot(aes(sample_id, area)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    ylab("Cell Area (pixels)") +
    xlab("")

ggsave(
    filename = "./figure/01_preprocessing/05_QC_coverage_cellsize.tiff",
    plot = p_coverage / p_cellsize,
    width = 12,
    height = 8,
    dpi = 300,
    compression = "lzw"
)

################################################################################
# 12. Filter Small Cells
################################################################################

# NOTE: Very small cells are likely debris or segmentation artifacts
# Remove cells with area < threshold
n_before <- ncol(spe)
spe <- spe[, spe$area >= PREPROCESSING_PARAMS$min_cell_area]
n_after <- ncol(spe)

# Log filtering results (optional: write to log file)
# cat(sprintf("Filtered %d cells (%.2f%%) with area < %d pixels\n",
#             n_before - n_after,
#             100 * (n_before - n_after) / n_before,
#             PREPROCESSING_PARAMS$min_cell_area))

################################################################################
# 13. Spatial Visualization (Optional)
################################################################################

# Select random samples for visualization
set.seed(PREPROCESSING_PARAMS$seed)
sample_ids <- unique(spe$sample_id)
cur_id <- sample(
    sample_ids, 
    min(PREPROCESSING_PARAMS$n_samples_spatial, length(sample_ids))
)

images_subset <- images[names(images) %in% cur_id]
masks_subset <- masks[names(masks) %in% cur_id]

# Normalize images for visualization
images_subset <- normalize(images_subset)
images_subset <- normalize(images_subset, inputRange = c(0, 0.2))

# Plot pixel intensities for key markers
# NOTE: Only plot if markers exist in dataset
if (all(c("CD68", "CD4", "EpCAM") %in% rownames(spe))) {
    tiff(
        filename = "./figure/01_preprocessing/06_spatial_pixels.tiff",
        width = 12,
        height = 4,
        units = "in",
        res = 300,
        compression = "lzw"
    )
    plotPixels(
        images_subset,
        colour_by = c("CD68", "CD4", "EpCAM"),
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

################################################################################
# 14. Save Objects
################################################################################

saveRDS(spe, "./data/spe_preprocessed.rds")
saveRDS(images, "./data/images.rds")
saveRDS(masks, "./data/masks.rds")

# Clean up
rm(image_mean, images_subset, masks_subset)
gc()

# NOTE: Next step is 02_batch_correction.R
