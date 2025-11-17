################################################################################
# Step 06: Advanced Visualization and Publication Figures
#
# Purpose: Generate publication-quality figures and comprehensive visualizations
# Dataset: NMIBC bladder cancer CyTOF
#
# Input:
#   - ./data/sce_annotated.RData
#   - ./data/sce_myeloid_annotated.RData
#
# Output:
#   - ./results/06_visualization/*.png (publication figures)
#
# Author: YMS
# Date: 2025
################################################################################

library(CATALYST)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(gtools)
library(ggdendro)

# Create output directory
if (!dir.exists("./results/06_visualization")) {
    dir.create("./results/06_visualization", recursive = TRUE)
}

################################################################################
# 1. Load Data
################################################################################

load("./data/sce_annotated.RData")
load("./data/sce_myeloid_annotated.RData")

################################################################################
# 2. Feature Plots: Marker Expression on t-SNE
################################################################################

# NOTE: Generate individual feature plots for each marker
# Useful for manuscript supplementary figures

feature_markers <- c(
    'CD45', 'CD3e', 'CD4', 'CD8', 'CD11b', 'CD11c',
    'CD14', 'CD15', 'CD16', 'CD20', 'CD32', 'CD36',
    'CD38', 'CD64', 'CD68', 'CD192'
)

# Custom color palette (gray to blue)
feature_colors <- colorRampPalette(c("grey90", "#0033FF"))(100)

for (marker in feature_markers) {
    p <- plotDR(
        sce_fil,
        "TSNE",
        color_by = marker,
        a_pal = feature_colors
    ) +
        ggtitle(marker) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            legend.position = "right"
        )
    
    ggsave(
        filename = paste0(
            "./results/06_visualization/featureplot_",
            marker,
            ".png"
        ),
        plot = p,
        width = 7,
        height = 6,
        dpi = 300
    )
}

################################################################################
# 3. Stacked Barplot: Cell Abundance
################################################################################

# NOTE: Create ordered stacked barplot with condition annotation
# Shows cell type proportions across all samples

# Get abundance data
abundance <- plotAbundances(
    sce_fil,
    k = "celltype",
    by = "sample_id",
    group_by = NULL
)

plot_data <- abundance$data

# Order samples (Control first, then NMIBC sorted naturally)
control_samples <- unique(
    as.character(plot_data$sample_id[grepl("^Control", plot_data$sample_id)])
)
nmibc_samples <- unique(
    as.character(plot_data$sample_id[grepl("^NMIBC", plot_data$sample_id)])
)

control_sorted <- mixedsort(control_samples)
nmibc_sorted <- mixedsort(nmibc_samples)
ordered_samples <- c(control_sorted, nmibc_sorted)

plot_data$sample_id <- factor(plot_data$sample_id, levels = ordered_samples)

# Prepare condition annotation
condition_data <- metadata(sce_fil)$experiment_info %>%
    filter(sample_id %in% ordered_samples) %>%
    select(sample_id, condition) %>%
    distinct() %>%
    mutate(sample_id = factor(sample_id, levels = ordered_samples))

# Define condition colors
condition_colors <- c(
    "ctrl" = "#4DAF4A",
    "UR"   = "#377EB8",
    "RE"   = "#E41A1C",
    "RC"   = "#984EA3",
    "DX"   = "#FF7F00"
)

# Create condition annotation bar
p_anno <- ggplot(
    condition_data,
    aes(x = sample_id, y = "Condition", fill = condition)
) +
    geom_tile(color = NA) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    scale_x_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        plot.margin = margin(0, 0, 3, 0)
    )

# Create stacked barplot
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

p_bar <- ggplot(
    plot_data,
    aes(x = sample_id, y = Freq, fill = cluster_id)
) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_minimal() +
    theme(
        axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(0, 0, 5, 0),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10)
    )

# Combine plots
p_final <- p_anno / p_bar +
    plot_layout(heights = c(0.5, 10))

ggsave(
    filename = "./results/06_visualization/stacked_barplot_ordered.png",
    plot = p_final,
    width = 14,
    height = 8,
    dpi = 300
)

################################################################################
# 4. Hierarchical Clustering of Samples
################################################################################

# NOTE: Cluster samples by cell type composition
# Creates dendrogram based on abundance patterns

# Convert to wide format for clustering
abundance_wide <- plot_data %>%
    pivot_wider(names_from = cluster_id, values_from = Freq, values_fill = 0)

mat <- as.matrix(abundance_wide[, -1])
rownames(mat) <- abundance_wide$sample_id

# Hierarchical clustering
d <- dist(mat, method = "euclidean")
hc <- hclust(d, method = "ward.D2")
sample_order <- hc$order
ordered_samples_clust <- rownames(mat)[sample_order]

# Update plot data with clustered order
plot_data$sample_id <- factor(
    plot_data$sample_id,
    levels = ordered_samples_clust
)

condition_data$sample_id <- factor(
    condition_data$sample_id,
    levels = ordered_samples_clust
)

# Create dendrogram plot
dend_data <- dendro_data(as.dendrogram(hc))

p_dend <- ggplot(segment(dend_data)) +
    geom_segment(
        aes(x = x, y = y, xend = xend, yend = yend),
        linewidth = 0.5
    ) +
    scale_x_continuous(expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(5, 5, 0, 5))

# Recreate annotation and barplot with new order
p_anno_clust <- ggplot(
    condition_data,
    aes(x = sample_id, y = "Condition", fill = condition)
) +
    geom_tile(color = NA) +
    scale_fill_manual(values = condition_colors, name = "Condition") +
    scale_x_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        plot.margin = margin(0, 0, 3, 0)
    )

p_bar_clust <- ggplot(
    plot_data,
    aes(x = sample_id, y = Freq, fill = cluster_id)
) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion (%)") +
    theme_minimal() +
    theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
        panel.grid = element_blank(),
        legend.position = "right"
    )

# Combine with dendrogram
p_final_dend <- p_dend / p_anno_clust / p_bar_clust +
    plot_layout(heights = c(1.5, 0.3, 10), guides = "collect") &
    theme(legend.position = "right")

ggsave(
    filename = "./results/06_visualization/stacked_barplot_clustered.png",
    plot = p_final_dend,
    width = 16,
    height = 9,
    dpi = 300
)

################################################################################
# 5. Multi-Panel Figure: Overview
################################################################################

# NOTE: Create comprehensive overview figure for publication
# Combines UMAP, t-SNE, heatmap, and abundance

# Panel A: UMAP by cell type
p_umap <- plotDR(
    sce_fil,
    "UMAP",
    color_by = "celltype"
) +
    scale_color_manual(values = celltype_colors) +
    ggtitle("A. UMAP: Cell Types") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        legend.position = "none"
    )

# Panel B: t-SNE by condition
p_tsne <- plotDR(
    sce_fil,
    "TSNE",
    color_by = "condition"
) +
    ggtitle("B. t-SNE: Conditions") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        legend.position = "right"
    )

# Combine panels
p_overview <- (p_umap | p_tsne) / p_bar

ggsave(
    filename = "./results/06_visualization/overview_figure.png",
    plot = p_overview,
    width = 18,
    height = 12,
    dpi = 300
)

################################################################################
# 6. Myeloid Cell Visualization
################################################################################

# t-SNE of myeloid subtypes by condition
myeloid_tsne <- plotDR(
    sce_myeloid,
    "TSNE",
    color_by = "myeloid_subtype",
    facet_by = "condition"
) +
    ggtitle("Myeloid Subtypes Across Conditions") +
    theme_bw() +
    theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank()
    )

ggsave(
    filename = "./results/06_visualization/myeloid_TSNE_by_condition.png",
    plot = myeloid_tsne,
    width = 16,
    height = 10,
    dpi = 300
)

################################################################################
# 7. Clean Up
################################################################################

gc()

################################################################################
# End of Step 06
################################################################################
