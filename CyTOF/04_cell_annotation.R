################################################################################
# Step 04: Cell Type Annotation
#
# Purpose: Annotate clusters based on marker expression
# Dataset: NMIBC bladder cancer CyTOF
#
# Input:
#   - ./data/sce_clustered.RData
#   - ./cluster_annotation.xlsx (cluster annotation table)
#
# Output:
#   - ./data/sce_annotated.RData
#   - ./results/04_annotation/*.png (annotated visualizations)
#
# Author: YMS
# Date: 2025
################################################################################

library(CATALYST)
library(dplyr)
library(ggplot2)
library(readxl)
library(scater)

# Create output directory
if (!dir.exists("./results/04_annotation")) {
    dir.create("./results/04_annotation", recursive = TRUE)
}

################################################################################
# 1. Load Clustered Data
################################################################################

load("./data/sce_clustered.RData")

################################################################################
# 2. Apply Cluster Annotations
################################################################################

# NOTE: Cluster annotation based on canonical marker expression
# Load pre-defined annotation table mapping cluster IDs to cell types

merging_table <- read_excel("./cluster_annotation.xlsx")

# Apply cluster annotations
sce <- mergeClusters(
    sce,
    k = "meta18",
    table = merging_table,
    id = "celltype",
    overwrite = TRUE
)

################################################################################
# 3. Filter Cell Types (Optional)
################################################################################

# NOTE: Remove low-quality or unwanted cell populations if needed
# Example: Keep only well-defined immune cell types

keep_celltypes <- c(
    "CD16- Neutrophils",
    "NK cells",
    "CD15+ NK cells",
    "CD8+ T cells",
    "CD16+ Neutrophils",
    "B cells",
    "pDCs",
    "Plasma cells",
    "cDCs",
    "CD36+CD38+CD4+ T cells",
    "CD4+ T cells",
    "DP T cells",
    "Monocytes"
)

sce_fil <- filterSCE(
    sce,
    cluster_id %in% keep_celltypes,
    k = "celltype"
)

################################################################################
# 4. Define Cell Type Colors
################################################################################

# NOTE: Consistent color palette for visualizations
celltype_colors <- c(
    "B cells"                = "#a9a9a9",
    "CD15+ NK cells"         = "#ffe119",
    "CD16- Neutrophils"      = "#f58231",
    "CD16+ Neutrophils"      = "#000075",
    "CD36+CD38+CD4+ T cells" = "#911eb4",
    "CD4+ T cells"           = "#4363d8",
    "CD8+ T cells"           = "#bfef45",
    "DP T cells"             = "#3cb44b",
    "cDCs"                   = "#e6194b",
    "Monocytes"              = "#f032e6",
    "Plasma cells"           = "#fabebe",
    "NK cells"               = "#008080",
    "pDCs"                   = "#9a6324"
)

################################################################################
# 5. Generate Annotated Visualizations
################################################################################

# t-SNE colored by cell type
tsne_celltype <- plotDR(
    sce_fil,
    "TSNE",
    color_by = "celltype"
) +
    scale_color_manual(values = celltype_colors) +
    ggtitle("t-SNE: Cell Types") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank()
    )

ggsave(
    filename = "./results/04_annotation/TSNE_celltype.png",
    plot = tsne_celltype,
    width = 10,
    height = 7,
    dpi = 300
)

# t-SNE faceted by condition
tsne_condition <- plotDR(
    sce_fil,
    "TSNE",
    color_by = "celltype",
    facet_by = "condition"
) +
    scale_color_manual(values = celltype_colors) +
    ggtitle("t-SNE: Cell Types by Condition") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank()
    )

ggsave(
    filename = "./results/04_annotation/TSNE_celltype_by_condition.png",
    plot = tsne_condition,
    width = 16,
    height = 10,
    dpi = 300
)

# Expression heatmap with cell type annotations
expr_heatmap <- plotExprHeatmap(
    sce_fil,
    features = "type",
    by = "cluster_id",
    k = "celltype",
    col_clust = FALSE,
    bars = TRUE,
    perc = TRUE
)

ggsave(
    filename = "./results/04_annotation/expression_heatmap_celltype.png",
    plot = expr_heatmap,
    width = 12,
    height = 10,
    dpi = 300
)

# Frequency heatmap
freq_heatmap <- plotFreqHeatmap(
    sce_fil,
    k = "celltype",
    perc = TRUE,
    row_anno = FALSE,
    col_clust = FALSE
)

ggsave(
    filename = "./results/04_annotation/frequency_heatmap.png",
    plot = freq_heatmap,
    width = 10,
    height = 8,
    dpi = 300
)

################################################################################
# 6. Cell Abundance Analysis
################################################################################

# Cell abundance by sample
abundance_sample <- plotAbundances(
    sce_fil,
    k = "celltype",
    by = "sample_id",
    group_by = NULL
) +
    scale_fill_manual(values = celltype_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
    filename = "./results/04_annotation/abundance_by_sample.png",
    plot = abundance_sample,
    width = 14,
    height = 7,
    dpi = 300
)

# Cell abundance by cluster
abundance_cluster <- plotAbundances(
    sce_fil,
    k = "celltype",
    by = "cluster_id"
)

ggsave(
    filename = "./results/04_annotation/abundance_by_cluster.png",
    plot = abundance_cluster,
    width = 10,
    height = 7,
    dpi = 300
)

################################################################################
# 7. Myeloid Cell Subclustering
################################################################################

# NOTE: Myeloid cells (monocytes, DCs, neutrophils) require finer resolution
# Extract and re-cluster for subtype identification

myeloid_types <- c(
    "Monocytes",
    "cDCs",
    "pDCs",
    "CD16- Neutrophils",
    "CD16+ Neutrophils"
)

sce_myeloid <- filterSCE(
    sce,
    cluster_id %in% myeloid_types,
    k = "celltype"
)

# Re-cluster myeloid cells at higher resolution
set.seed(900223)
sce_myeloid <- cluster(
    sce_myeloid,
    features = "type",
    xdim = 10,
    ydim = 10,
    maxK = 20,
    seed = 900223
)

# Re-run dimensionality reduction
set.seed(900223)
sce_myeloid <- runDR(
    sce_myeloid,
    "TSNE",
    cells = 5000,
    features = "type"
)

# Plot myeloid subclusters
myeloid_heatmap <- plotExprHeatmap(
    sce_myeloid,
    features = "type",
    by = "cluster_id",
    k = "meta13",
    bars = TRUE,
    perc = TRUE
)

ggsave(
    filename = "./results/04_annotation/myeloid_subcluster_heatmap.png",
    plot = myeloid_heatmap,
    width = 12,
    height = 10,
    dpi = 300
)

# t-SNE of myeloid subclusters
myeloid_tsne <- plotDR(
    sce_myeloid,
    "TSNE",
    color_by = "meta13"
) +
    ggtitle("Myeloid Subclustering") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave(
    filename = "./results/04_annotation/myeloid_TSNE.png",
    plot = myeloid_tsne,
    width = 10,
    height = 7,
    dpi = 300
)

################################################################################
# 8. Apply Myeloid Subcluster Annotations
################################################################################

# NOTE: Apply fine-grained annotations to myeloid subclusters
# Example mapping (adjust based on marker expression):

myeloid_merging <- read_excel("./myeloid_annotation.xlsx")

sce_myeloid <- mergeClusters(
    sce_myeloid,
    k = "meta13",
    table = myeloid_merging,
    id = "myeloid_subtype",
    overwrite = TRUE
)

# Visualize myeloid subtypes
myeloid_subtype_tsne <- plotDR(
    sce_myeloid,
    "TSNE",
    color_by = "myeloid_subtype"
) +
    ggtitle("Myeloid Subtypes") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave(
    filename = "./results/04_annotation/myeloid_subtypes_TSNE.png",
    plot = myeloid_subtype_tsne,
    width = 10,
    height = 7,
    dpi = 300
)

################################################################################
# 9. Save Annotated Objects
################################################################################

save(sce_fil, file = "./data/sce_annotated.RData")
save(sce_myeloid, file = "./data/sce_myeloid_annotated.RData")

# Clean up
gc()

################################################################################
# End of Step 04
################################################################################
