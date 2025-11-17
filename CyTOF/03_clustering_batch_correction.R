################################################################################
# Step 03: Clustering and Batch Correction
#
# Purpose: FlowSOM clustering and Harmony batch correction
# Dataset: NMIBC bladder cancer CyTOF
#
# Input:
#   - ./data/sce_integrated.RData
#
# Output:
#   - ./data/sce_clustered.RData
#   - ./results/03_clustering/*.png (heatmaps, DR plots)
#
# Author: YMS
# Date: 2025
################################################################################

library(CATALYST)
library(scater)
library(BiocSingular)
library(harmony)
library(ggplot2)

# Create output directory
if (!dir.exists("./results/03_clustering")) {
    dir.create("./results/03_clustering", recursive = TRUE)
}

################################################################################
# 1. Load Integrated Data
################################################################################

load("./data/sce_integrated.RData")

################################################################################
# 2. FlowSOM Clustering
################################################################################

# NOTE: FlowSOM uses self-organizing maps for clustering
# Uses type markers (use_channel == TRUE) for cell identification

set.seed(900223)
sce <- cluster(
    sce,
    features = "type",
    xdim = 10,
    ydim = 10,
    maxK = 50,
    seed = 900223
)

# Print available cluster resolutions
cat("=== Available Cluster Resolutions ===\n")
print(names(cluster_codes(sce)))

################################################################################
# 3. Determine Expression Assay for Downstream Analysis
################################################################################

# NOTE: If batch correction was applied (e.g., cyCombine), use corrected assay
# Otherwise use standard 'exprs' assay

expr_assay <- if ("exprs_bc" %in% assayNames(sce)) {
    "exprs_bc"
} else {
    "exprs"
}

cat("Using assay for analysis:", expr_assay, "\n")

################################################################################
# 4. Define Features for PCA
################################################################################

# NOTE: Use type markers if available, otherwise all use_channel markers
has_type <- "marker_class" %in% colnames(rowData(sce)) &&
    any(rowData(sce)$marker_class == "type", na.rm = TRUE)

features_for_pca <- if (has_type) {
    rownames(sce)[
        rowData(sce)$marker_class == "type" &
        rowData(sce)$use_channel
    ]
} else {
    rownames(sce)[rowData(sce)$use_channel]
}

cat("Features for PCA:", length(features_for_pca), "\n")

################################################################################
# 5. PCA Dimensionality Reduction
################################################################################

set.seed(900223)
sce <- scater::runPCA(
    sce,
    subset_row   = features_for_pca,
    exprs_values = expr_assay,
    ncomponents  = 30,
    BSPARAM      = BiocSingular::ExactParam()
)

################################################################################
# 6. Harmony Batch Correction
################################################################################

# NOTE: Harmony corrects PCA embeddings, not expression values
# This preserves original marker intensities while removing batch effects

if (!("batch" %in% colnames(colData(sce)))) {
    stop("colData(sce)$batch is missing. Check Step 02.")
}

# Extract PCA coordinates
Z <- reducedDim(sce, "PCA")

# Prepare batch metadata
meta_df <- data.frame(batch = as.factor(colData(sce)$batch))

# Apply Harmony correction
Z_harmony <- harmony::HarmonyMatrix(
    data_mat  = Z,
    meta_data = meta_df,
    vars_use  = "batch",
    do_pca    = FALSE
)

# Store Harmony-corrected embeddings
reducedDim(sce, "PCA_harmony") <- Z_harmony

################################################################################
# 7. Generate UMAP and t-SNE on Harmony Embeddings
################################################################################

# Harmony-corrected DR (for visualization)
set.seed(900223)
sce <- scater::runUMAP(
    sce,
    dimred = "PCA_harmony",
    name = "UMAP_harmony"
)

set.seed(900223)
sce <- scater::runTSNE(
    sce,
    dimred = "PCA_harmony",
    name = "TSNE_harmony"
)

# Original DR for comparison (downsampled for speed)
set.seed(900223)
sce <- runDR(sce, "TSNE", cells = 1000, features = "type")

set.seed(900223)
sce <- runDR(sce, "UMAP", cells = 1000, features = "type")

################################################################################
# 8. Initial Visualization
################################################################################

# Define canonical immune markers for visualization
total_marker <- c(
    'CD1c', 'CD3e', 'CD7', 'CD8', 'CD38',
    'HLA-DR', 'CD123', 'CD4', 'CD11c', 'CD68',
    'CD36', 'CD64', 'CD14', 'CD11b', 'CD16',
    'CD66ace', 'CD15', 'CD20'
)

feat_ok <- intersect(total_marker, rownames(sce))

# Expression heatmap by cluster
if (length(feat_ok) > 0) {
    heatmap_plot <- plotExprHeatmap(
        sce,
        features = feat_ok,
        by = "cluster_id",
        k = "meta34",
        col_clust = FALSE,
        bars = TRUE,
        perc = TRUE,
        assay = expr_assay
    )
    
    ggsave(
        filename = "./results/03_clustering/expression_heatmap_meta34.png",
        plot = heatmap_plot,
        width = 12,
        height = 10,
        dpi = 300
    )
}

# UMAP colored by batch (Harmony-corrected)
umap_batch <- plotDR(
    sce,
    "UMAP_harmony",
    color_by = "batch"
) +
    ggtitle("UMAP (Harmony-corrected) by Batch") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave(
    filename = "./results/03_clustering/UMAP_harmony_batch.png",
    plot = umap_batch,
    width = 10,
    height = 7,
    dpi = 300
)

# UMAP colored by cluster
umap_cluster <- plotDR(
    sce,
    "UMAP_harmony",
    color_by = "meta34"
) +
    ggtitle("UMAP (Harmony-corrected) by Cluster") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave(
    filename = "./results/03_clustering/UMAP_harmony_cluster.png",
    plot = umap_cluster,
    width = 10,
    height = 7,
    dpi = 300
)

# t-SNE colored by cluster (Harmony-corrected)
tsne_cluster <- plotDR(
    sce,
    "TSNE_harmony",
    color_by = "meta34"
) +
    ggtitle("t-SNE (Harmony-corrected) by Cluster") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave(
    filename = "./results/03_clustering/TSNE_harmony_cluster.png",
    plot = tsne_cluster,
    width = 10,
    height = 7,
    dpi = 300
)

# t-SNE colored by batch
tsne_batch <- plotDR(
    sce,
    "TSNE_harmony",
    color_by = "batch"
) +
    ggtitle("t-SNE (Harmony-corrected) by Batch") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave(
    filename = "./results/03_clustering/TSNE_harmony_batch.png",
    plot = tsne_batch,
    width = 10,
    height = 7,
    dpi = 300
)

# Multi-heatmap (expression + frequency)
multi_heatmap <- plotMultiHeatmap(
    sce,
    hm1 = "type",
    hm2 = "state",
    k = "meta34",
    row_clust = TRUE,
    col_clust = FALSE,
    row_anno = FALSE,
    bars = TRUE,
    perc = TRUE
)

ggsave(
    filename = "./results/03_clustering/multi_heatmap.png",
    plot = multi_heatmap,
    width = 14,
    height = 10,
    dpi = 300
)

################################################################################
# 9. Save Clustered Object
################################################################################

save(sce, file = "./data/sce_clustered.RData")

# Clean up
gc()

################################################################################
# End of Step 03
################################################################################
