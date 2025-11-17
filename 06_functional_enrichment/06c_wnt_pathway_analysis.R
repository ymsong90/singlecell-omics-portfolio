################################################################################
# Step 06c: WNT Pathway-Specific Analysis
#
# Purpose: Analyze WNT signaling components in myeloid cells
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./data/Myeloid_subset_filtered.RData
#
# Output:
#   - ./results/06_WNT/*.csv (expression data)
#   - ./results/06_WNT/*.png (visualizations)
#
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Create output directory
dir.create("./results/06_WNT", recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Load Myeloid Subset
################################################################################

load("./data/Myeloid_subset_filtered.RData")

# Focus on WT for baseline expression
myeloid_wt <- subset(Myeloid_subset_fil, subset = ID == "WT")

################################################################################
# 2. Define WNT Pathway Genes
################################################################################

# WNT ligands
wnt_ligands <- c(
    "Wnt1", "Wnt2", "Wnt2b", "Wnt3", "Wnt3a", "Wnt4", "Wnt5a", "Wnt5b",
    "Wnt6", "Wnt7a", "Wnt7b", "Wnt8a", "Wnt8b", "Wnt9a", "Wnt9b",
    "Wnt10a", "Wnt10b", "Wnt11", "Wnt16"
)

# Receptors
wnt_receptors <- c(
    "Lrp1", "Lrp5", "Lrp6",
    "Fzd1", "Fzd2", "Fzd3", "Fzd4", "Fzd5", "Fzd6", "Fzd7", "Fzd8", "Fzd9", "Fzd10",
    "Ror1", "Ror2", "Ryk"
)

# Transcription factors
wnt_tfs <- c("Tcf3", "Tcf4", "Tcf7", "Lef1")

# All WNT pathway genes
wnt_genes <- unique(c(wnt_ligands, wnt_receptors, wnt_tfs))

################################################################################
# 3. Calculate Average Expression per Cluster
################################################################################

# Check which genes are present in the dataset
wnt_genes_present <- intersect(wnt_genes, rownames(myeloid_wt))

# Extract expression data
expr_data <- GetAssayData(myeloid_wt, slot = "data")[wnt_genes_present, ]

# Calculate mean expression per cluster
Idents(myeloid_wt) <- "Myeloid_labels"
cluster_names <- levels(Idents(myeloid_wt))

mean_expr <- data.frame(Gene = wnt_genes_present)

for (cluster in cluster_names) {
    cluster_cells <- WhichCells(myeloid_wt, idents = cluster)
    
    if (length(cluster_cells) > 0) {
        cluster_mean <- rowMeans(expr_data[, cluster_cells, drop = FALSE])
        mean_expr[[cluster]] <- cluster_mean
    }
}

# Reshape for plotting
mean_expr_long <- mean_expr %>%
    pivot_longer(
        cols = -Gene,
        names_to = "Cluster",
        values_to = "Expression"
    ) %>%
    mutate(
        Category = case_when(
            Gene %in% wnt_ligands ~ "Ligand",
            Gene %in% wnt_receptors ~ "Receptor",
            Gene %in% wnt_tfs ~ "Transcription Factor",
            TRUE ~ "Other"
        )
    )

write.csv(
    mean_expr_long,
    "./results/06_WNT/WNT_pathway_expression_by_cluster.csv",
    row.names = FALSE
)

################################################################################
# 4. Visualization - Heatmap
################################################################################

# Prepare matrix for heatmap
expr_matrix <- mean_expr %>%
    column_to_rownames("Gene") %>%
    as.matrix()

# Heatmap
library(pheatmap)

pheatmap(
    expr_matrix,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    filename = "./results/06_WNT/WNT_pathway_heatmap.png",
    width = 8,
    height = 10
)

################################################################################
# 5. Visualization - Dot Plot
################################################################################

# NOTE: Show only highly expressed genes for clarity
top_genes <- mean_expr_long %>%
    group_by(Gene) %>%
    summarise(max_expr = max(Expression)) %>%
    filter(max_expr > 0.5) %>%
    pull(Gene)

plot_data <- mean_expr_long %>%
    filter(Gene %in% top_genes)

p_dot <- ggplot(plot_data, aes(x = Gene, y = Cluster, size = Expression, color = Category)) +
    geom_point() +
    scale_size_continuous(range = c(1, 8)) +
    scale_color_manual(values = c(
        "Ligand" = "#e74c3c",
        "Receptor" = "#3498db",
        "Transcription Factor" = "#2ecc71"
    )) +
    theme_minimal(base_size = 12) +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold")
    ) +
    labs(
        x = "WNT Pathway Gene",
        y = "Myeloid Subtype",
        title = "WNT Pathway Expression in Myeloid Cells"
    )

ggsave(
    "./results/06_WNT/WNT_pathway_dotplot.png",
    plot = p_dot,
    width = 12,
    height = 6,
    dpi = 300
)

################################################################################
# 6. Compare KO vs WT
################################################################################

# Extract KO samples
myeloid_ko <- subset(Myeloid_subset_fil, subset = ID == "KO")

# Calculate fold change for key WNT genes
# ... (implementation for specific comparisons)

# Clean up
rm(myeloid_wt, myeloid_ko, expr_data, mean_expr, mean_expr_long, plot_data)
gc()

################################################################################
# End of Step 06c
################################################################################
