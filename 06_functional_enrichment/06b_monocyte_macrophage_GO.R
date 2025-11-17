################################################################################
# Step 06b: Monocyte/Macrophage-Specific GO Analysis
#
# Purpose: Detailed functional analysis of monocyte and macrophage subtypes
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./results/05_DEG/myeloid/DEG_*.csv (cluster-specific)
#
# Output:
#   - ./results/06_GO/monocyte_macrophage/*.csv
#   - ./results/06_GO/monocyte_macrophage/*.png
#
# Author: YMS
# Date: 2025
################################################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Create output directory
dir.create("./results/06_GO/monocyte_macrophage", recursive = TRUE, showWarnings = FALSE)

################################################################################
# 1. Load Cluster-Specific DEGs
################################################################################

# NOTE: Focus on monocyte and macrophage populations
clusters_of_interest <- c(
    "Classical_Monocyte",
    "Activated_Monocyte",
    "Trem2+_Macrophage",
    "Selenop+_Macrophage",
    "Hexb+_Macrophage"
)

################################################################################
# 2. Cluster-Specific GO Enrichment
################################################################################

for (cluster in clusters_of_interest) {
    
    deg_file <- paste0("./results/05_DEG/myeloid/DEG_", cluster, "_KO_vs_WT.csv")
    
    if (!file.exists(deg_file)) next
    
    deg <- read.csv(deg_file)
    
    # Filter significant genes
    sig_genes <- deg %>%
        filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
        pull(gene)
    
    if (length(sig_genes) < 10) next
    
    # Convert to Entrez IDs
    entrez_ids <- bitr(
        sig_genes,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Mm.eg.db
    )$ENTREZID
    
    # GO enrichment
    go_result <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        readable = TRUE
    )
    
    # Save results
    if (nrow(as.data.frame(go_result)) > 0) {
        write.csv(
            as.data.frame(go_result),
            paste0("./results/06_GO/monocyte_macrophage/GO_", cluster, ".csv"),
            row.names = FALSE
        )
        
        # Plot
        p <- dotplot(go_result, showCategory = 15) +
            ggtitle(paste("GO Enrichment -", gsub("_", " ", cluster)))
        
        ggsave(
            paste0("./results/06_GO/monocyte_macrophage/GO_", cluster, "_dotplot.png"),
            plot = p, width = 12, height = 8, dpi = 300
        )
    }
}

################################################################################
# 3. Compare Monocyte vs Macrophage GO Terms
################################################################################

# NOTE: This comparison helps identify functional differences
# between monocytes (inflammatory) and macrophages (tissue-resident)

# Load and combine results for comparison
# ... (implementation depends on specific research questions)

# Clean up
gc()

################################################################################
# End of Step 06b
################################################################################
