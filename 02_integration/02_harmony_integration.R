################################################################################
# Step 02: Normalization and Harmony Integration
#
# Purpose: Normalize expression data and perform batch correction
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./data/porcn.combined_QC.RData
#
# Output:
#   - ./data/porcn.combined.harmony.RData
#   - ./results/02_integration/*.png
#
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

# Parameters
PARAMS <- list(
    normalization_method = "LogNormalize",
    scale_factor         = 10000,
    n_variable_features  = 2000,
    n_pcs                = 50,      # Compute 50 PCs
    use_pcs              = 1:30,    # Use first 30 for clustering
    cluster_resolution   = 2.5
)

# Create output directory
if (!dir.exists("./results/02_integration")) {
    dir.create("./results/02_integration", recursive = TRUE)
}

################################################################################
# 1. Load QC-Filtered Data
################################################################################

load("./data/porcn.combined_QC.RData")

################################################################################
# 2. Normalization
################################################################################

# NOTE: LogNormalize assumes each cell originally contains
# the same number of molecules (total counts)
porcn.combined <- NormalizeData(
    object               = porcn.combined,
    normalization.method = PARAMS$normalization_method,
    scale.factor         = PARAMS$scale_factor
)

################################################################################
# 3. Identify Highly Variable Features
################################################################################

porcn.combined <- FindVariableFeatures(
    porcn.combined,
    selection.method = "vst",
    nfeatures        = PARAMS$n_variable_features
)

# Plot variable features
top_genes <- head(VariableFeatures(porcn.combined), 20)

plot_var_features <- VariableFeaturePlot(porcn.combined)
plot_var_labeled <- LabelPoints(
    plot   = plot_var_features,
    points = top_genes,
    repel  = TRUE
)

ggsave(
    filename = "./results/02_integration/variable_features.png",
    plot     = plot_var_labeled,
    width    = 10,
    height   = 7,
    dpi      = 300
)

################################################################################
# 4. Scale Data
################################################################################

# NOTE: Scaling is essential for PCA
# - Centers expression to mean 0
# - Scales to unit variance
# - Prevents highly expressed genes from dominating
porcn.combined <- ScaleData(
    object   = porcn.combined,
    features = rownames(porcn.combined)
)

################################################################################
# 5. PCA
################################################################################

porcn.combined <- RunPCA(
    porcn.combined,
    features = VariableFeatures(object = porcn.combined),
    npcs     = PARAMS$n_pcs,
    verbose  = FALSE
)

# Elbow plot to determine optimal PC number
elbow_plot <- ElbowPlot(porcn.combined, ndims = PARAMS$n_pcs) +
    ggtitle("PCA Elbow Plot") +
    geom_vline(xintercept = max(PARAMS$use_pcs), 
               linetype = "dashed", color = "red")

ggsave(
    filename = "./results/02_integration/PCA_elbow_plot.png",
    plot     = elbow_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
)

# PCA visualization before Harmony
pca_before <- DimPlot(
    porcn.combined,
    reduction = "pca",
    group.by  = "ID",
    pt.size   = 0.5
) +
    ggtitle("PCA (Before Harmony)")

################################################################################
# 6. Baseline Clustering (Without Harmony)
################################################################################

# NOTE: This is for comparison purposes
# Helps visualize batch effects before correction
porcn.combined <- FindNeighbors(
    object = porcn.combined,
    dims   = PARAMS$use_pcs
)

porcn.combined <- FindClusters(
    object     = porcn.combined,
    resolution = PARAMS$cluster_resolution
)

porcn.combined <- RunUMAP(
    object = porcn.combined,
    dims   = PARAMS$use_pcs
)

# Store for comparison
porcn.combined[["clusters_no_harmony"]] <- Idents(object = porcn.combined)

umap_no_harmony <- DimPlot(
    object    = porcn.combined,
    reduction = "umap",
    group.by  = "ID",
    pt.size   = 0.5
) +
    ggtitle("UMAP Without Harmony")

################################################################################
# 7. Harmony Batch Correction
################################################################################

# NOTE: Harmony corrects for batch effects while preserving biological variation
# - Iteratively adjusts cell embeddings in PCA space
# - Does not modify raw/normalized counts
porcn.combined.harmony <- porcn.combined %>% 
    RunHarmony(
        "ID",                    # Batch variable (WT vs KO)
        plot_convergence = TRUE,
        max.iter.harmony = 20
    )

# Save Harmony convergence plot
harmony_plot <- porcn.combined.harmony@tools$RunHarmony$plot
ggsave(
    filename = "./results/02_integration/harmony_convergence.png",
    plot     = harmony_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
)

################################################################################
# 8. Clustering With Harmony-Corrected Embeddings
################################################################################

porcn.combined.harmony <- FindNeighbors(
    object    = porcn.combined.harmony,
    dims      = PARAMS$use_pcs,
    reduction = "harmony"
)

porcn.combined.harmony <- FindClusters(
    object     = porcn.combined.harmony,
    resolution = PARAMS$cluster_resolution
)

# UMAP using Harmony embeddings
porcn.combined.harmony <- RunUMAP(
    object    = porcn.combined.harmony,
    dims      = PARAMS$use_pcs,
    reduction = "harmony"
)

# t-SNE using Harmony embeddings
porcn.combined.harmony <- RunTSNE(
    object    = porcn.combined.harmony,
    dims      = PARAMS$use_pcs,
    reduction = "harmony"
)

################################################################################
# 9. Generate Visualizations
################################################################################

# UMAP with Harmony - clusters
umap_harmony_clusters <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    label     = TRUE,
    pt.size   = 0.8
) +
    ggtitle("UMAP With Harmony (Clusters)")

# UMAP with Harmony - by condition
umap_harmony_condition <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    group.by  = "ID",
    pt.size   = 0.5
) +
    ggtitle("UMAP With Harmony (Condition)")

# UMAP with Harmony - split by condition
umap_harmony_split <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "umap",
    split.by  = "ID",
    label     = TRUE,
    pt.size   = 0.5
) +
    ggtitle("UMAP With Harmony (Split)")

# t-SNE visualization
tsne_harmony <- DimPlot(
    object    = porcn.combined.harmony,
    reduction = "tsne",
    label     = TRUE,
    pt.size   = 0.7
) +
    ggtitle("t-SNE With Harmony")

# Comparison: Before vs After Harmony
comparison_plot <- (umap_no_harmony | umap_harmony_condition) /
                   (pca_before | tsne_harmony)

ggsave(
    filename = "./results/02_integration/integration_comparison.png",
    plot     = comparison_plot,
    width    = 16,
    height   = 12,
    dpi      = 300
)

# Save individual plots
ggsave(
    filename = "./results/02_integration/UMAP_harmony_clusters.png",
    plot     = umap_harmony_clusters,
    width    = 10,
    height   = 8,
    dpi      = 300
)

ggsave(
    filename = "./results/02_integration/UMAP_harmony_split.png",
    plot     = umap_harmony_split,
    width    = 14,
    height   = 6,
    dpi      = 300
)

################################################################################
# 10. Post-Integration QC Check
################################################################################

Idents(porcn.combined.harmony) <- 'ID'

vln_qc <- VlnPlot(
    porcn.combined.harmony,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol     = 3,
    pt.size  = 0
)

ggsave(
    filename = "./results/02_integration/QC_metrics_post_integration.png",
    plot     = vln_qc,
    width    = 12,
    height   = 5,
    dpi      = 300
)

################################################################################
# 11. Save Integrated Object
################################################################################

save(porcn.combined.harmony, file = "./data/porcn.combined.harmony.RData")

# Clean up
rm(porcn.combined, umap_no_harmony, umap_harmony_clusters, 
   umap_harmony_condition, umap_harmony_split, tsne_harmony,
   comparison_plot, vln_qc, pca_before, elbow_plot, 
   plot_var_features, plot_var_labeled)
gc()

################################################################################
# End of Step 02
################################################################################
