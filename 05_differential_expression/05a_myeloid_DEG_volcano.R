################################################################################
# Step 05a: Myeloid Cell Differential Expression Analysis
#
# Purpose: Identify DEGs between KO and WT in myeloid cell subclusters
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./data/Myeloid_subset_filtered.RData
#
# Output:
#   - ./results/05_DEG/myeloid/*.csv (DEG tables)
#   - ./results/05_DEG/myeloid/volcano/*.png
#
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Parameters
DEG_PARAMS <- list(
    logfc_threshold = 0.25,
    min_pct = 0.10,
    padj_cutoff = 0.05
)

# Create output directories
dir.create("./results/05_DEG/myeloid", recursive = TRUE, showWarnings = FALSE)
dir.create("./results/05_DEG/myeloid/volcano", recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Load Myeloid Subset
################################################################################

load("./data/Myeloid_subset_filtered.RData")

################################################################################
# 2. Prepare for DEG Analysis
################################################################################

# Set identity to combine cell type and condition
# NOTE: This creates groups like "Classical Monocyte_WT", "Classical Monocyte_KO"
Idents(Myeloid_subset_fil) <- "Myeloid_labels"

Myeloid_subset_fil$celltype_condition <- paste(
    Idents(Myeloid_subset_fil),
    Myeloid_subset_fil$ID,
    sep = "_"
)

Idents(Myeloid_subset_fil) <- "celltype_condition"

################################################################################
# 3. Run DEG Analysis for Each Myeloid Subtype
################################################################################

myeloid_subtypes <- c(
    "Classical Monocyte",
    "Activated Monocyte",
    "Trem2+ Macrophage",
    "Selenop+ Macrophage",
    "Hexb+ Macrophage"
)

deg_list <- list()

for (subtype in myeloid_subtypes) {
    
    ident_ko <- paste0(subtype, "_KO")
    ident_wt <- paste0(subtype, "_WT")
    
    # Check if both groups exist
    if (!ident_ko %in% levels(Idents(Myeloid_subset_fil)) | 
        !ident_wt %in% levels(Idents(Myeloid_subset_fil))) {
        next
    }
    
    # Find markers
    deg <- FindMarkers(
        Myeloid_subset_fil,
        ident.1 = ident_ko,
        ident.2 = ident_wt,
        min.pct = DEG_PARAMS$min_pct,
        logfc.threshold = DEG_PARAMS$logfc_threshold,
        test.use = "wilcox"
    )
    
    # Add gene names and cluster info
    deg <- deg %>%
        rownames_to_column("gene") %>%
        mutate(cluster = subtype)
    
    deg_list[[subtype]] <- deg
    
    # Save individual DEG table
    write.csv(
        deg,
        file = paste0("./results/05_DEG/myeloid/DEG_", gsub(" ", "_", subtype), "_KO_vs_WT.csv"),
        row.names = FALSE
    )
}

# Combine all DEGs
deg_all <- bind_rows(deg_list)

write.csv(
    deg_all,
    file = "./results/05_DEG/myeloid/Myeloid_DEG_all_subtypes.csv",
    row.names = FALSE
)

################################################################################
# 4. Generate Volcano Plots
################################################################################

# Volcano plot function
plot_volcano <- function(deg_df, title, 
                         lfc_cut = DEG_PARAMS$logfc_threshold, 
                         padj_cut = DEG_PARAMS$padj_cutoff,
                         label_n = 20) {
    
    # Add significance labels
    deg_df <- deg_df %>%
        mutate(
            log10P = -log10(p_val_adj + 1e-300),
            sig_flag = case_when(
                p_val_adj < padj_cut & avg_log2FC > lfc_cut ~ "Up",
                p_val_adj < padj_cut & avg_log2FC < -lfc_cut ~ "Down",
                TRUE ~ "NS"
            ),
            sig_flag = factor(sig_flag, levels = c("Up", "Down", "NS"))
        )
    
    # Select top genes to label
    top_genes <- deg_df %>%
        filter(sig_flag != "NS") %>%
        arrange(p_val_adj) %>%
        group_by(sig_flag) %>%
        slice_head(n = label_n/2) %>%
        ungroup()
    
    # Plot
    p <- ggplot(deg_df, aes(x = avg_log2FC, y = log10P, color = sig_flag)) +
        geom_point(size = 1.5, alpha = 0.6) +
        scale_color_manual(
            breaks = c("Up", "Down", "NS"),
            values = c("Up" = "#e74c3c", "Down" = "#3498db", "NS" = "grey50"),
            labels = c(
                paste0("Up (", sum(deg_df$sig_flag == "Up"), ")"),
                paste0("Down (", sum(deg_df$sig_flag == "Down"), ")"),
                "NS"
            )
        ) +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), 
                   linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(padj_cut), 
                   linetype = "dashed", alpha = 0.5) +
        geom_text_repel(
            data = top_genes,
            aes(label = gene),
            size = 3,
            max.overlaps = 20,
            show.legend = FALSE
        ) +
        labs(
            title = title,
            x = "log2 Fold Change (KO vs WT)",
            y = "-log10(Adjusted P-value)",
            color = "Regulation"
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right",
            panel.grid.minor = element_blank()
        )
    
    return(p)
}

# Generate volcano plots for each subtype
for (subtype in names(deg_list)) {
    
    deg_df <- deg_list[[subtype]]
    
    p <- plot_volcano(deg_df, title = paste(subtype, "- KO vs WT"))
    
    ggsave(
        filename = paste0("./results/05_DEG/myeloid/volcano/Volcano_", 
                         gsub(" ", "_", subtype), ".png"),
        plot = p,
        width = 10,
        height = 8,
        dpi = 300
    )
}

################################################################################
# 5. Filter Significant DEGs
################################################################################

# Upregulated genes
deg_up <- deg_all %>%
    filter(p_val_adj < DEG_PARAMS$padj_cutoff & 
           avg_log2FC > DEG_PARAMS$logfc_threshold) %>%
    arrange(cluster, desc(avg_log2FC))

write.csv(
    deg_up,
    file = "./results/05_DEG/myeloid/DEG_upregulated.csv",
    row.names = FALSE
)

# Downregulated genes
deg_down <- deg_all %>%
    filter(p_val_adj < DEG_PARAMS$padj_cutoff & 
           avg_log2FC < -DEG_PARAMS$logfc_threshold) %>%
    arrange(cluster, avg_log2FC)

write.csv(
    deg_down,
    file = "./results/05_DEG/myeloid/DEG_downregulated.csv",
    row.names = FALSE
)

# Clean up
rm(deg_list, deg_all, deg_up, deg_down)
gc()

################################################################################
# End of Step 05a
################################################################################
