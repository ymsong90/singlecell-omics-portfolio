################################################################################
# Step 04: Cell Proportion Analysis
#
# Purpose: Analyze cell type composition changes between KO and WT
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./data/porcn.combined.harmony_annotated.RData
#
# Output:
#   - ./results/04_visualization/proportion_*.png
#   - ./results/04_visualization/proportion_*.csv
#
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

# Create output directory
if (!dir.exists("./results/04_visualization")) {
    dir.create("./results/04_visualization", recursive = TRUE)
}

################################################################################
# 1. Load Annotated Data
################################################################################

load("./data/porcn.combined.harmony_annotated.RData")

################################################################################
# 2. Calculate Cell Type Proportions
################################################################################

# Extract metadata
meta_df <- porcn.combined.harmony@meta.data %>%
    as_tibble(rownames = "cell") %>%
    select(cell, ID, NH_labels, Complete_Labels)

# Calculate proportions by condition
prop_nh <- meta_df %>%
    count(ID, NH_labels) %>%
    group_by(ID) %>%
    mutate(prop = n / sum(n) * 100) %>%
    ungroup()

prop_complete <- meta_df %>%
    count(ID, Complete_Labels) %>%
    group_by(ID) %>%
    mutate(prop = n / sum(n) * 100) %>%
    ungroup()

# Save proportion tables
write.csv(prop_nh, 
          "./results/04_visualization/proportions_NH_labels.csv",
          row.names = FALSE)

write.csv(prop_complete,
          "./results/04_visualization/proportions_Complete_labels.csv",
          row.names = FALSE)

################################################################################
# 3. Statistical Testing
################################################################################

# NOTE: Chi-square test for overall distribution differences
# Use Fisher's exact test for specific cell types with low counts

# Overall chi-square test
contingency_table <- table(meta_df$ID, meta_df$Complete_Labels)
chisq_test <- chisq.test(contingency_table)

# Save test results
test_results <- data.frame(
    Test = "Chi-square",
    Statistic = chisq_test$statistic,
    P_value = chisq_test$p.value
)

write.csv(test_results,
          "./results/04_visualization/proportion_test_results.csv",
          row.names = FALSE)

################################################################################
# 4. Stacked Bar Plot
################################################################################

# Define color palette for Complete_Labels
celltype_colors <- c(
    "Basal cancer cell"      = "#0B6A42",
    "Classical cancer cell"  = "#22A1B1",
    "CD8 T cell"             = "#C2E9F9",
    "NK cell"                = "#8F7BB7",
    "CD4 T cell"             = "#0E7DC2",
    "CAF"                    = "#EA6D22",
    "DC"                     = "#F17D97",
    "Mono/Mac"               = "#ED4169",
    "Neutrophil"             = "#778899",
    "B cell"                 = "#808000",
    "Endothelial cell"       = "#F0B817"
)

# Stacked bar plot
p_stacked <- ggplot(prop_complete, aes(x = ID, y = prop, fill = Complete_Labels)) +
    geom_bar(stat = "identity", width = 0.6, color = "grey30") +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    scale_fill_manual(values = celltype_colors, name = "Cell Type") +
    labs(
        x = "Condition",
        y = "Proportion (%)",
        title = "Cell Type Composition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        legend.position = "right"
    )

ggsave(
    filename = "./results/04_visualization/stacked_barplot_proportion.png",
    plot = p_stacked,
    width = 10,
    height = 8,
    dpi = 300
)

################################################################################
# 5. Side-by-Side Comparison
################################################################################

# Calculate fold change (KO/WT)
prop_wide <- prop_complete %>%
    select(ID, Complete_Labels, prop) %>%
    pivot_wider(names_from = ID, values_from = prop, values_fill = 0) %>%
    mutate(
        log2FC = log2((KO + 0.01) / (WT + 0.01)),  # Add pseudocount
        direction = ifelse(log2FC > 0, "Increased", "Decreased")
    )

write.csv(prop_wide,
          "./results/04_visualization/proportion_fold_changes.csv",
          row.names = FALSE)

# Bar plot of fold changes
p_fc <- ggplot(prop_wide, aes(x = reorder(Complete_Labels, log2FC), 
                               y = log2FC, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("Increased" = "#e74c3c", "Decreased" = "#3498db")) +
    coord_flip() +
    labs(
        x = "Cell Type",
        y = "log2(Fold Change) KO/WT",
        title = "Cell Type Proportion Changes"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.grid.major.y = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "top"
    )

ggsave(
    filename = "./results/04_visualization/proportion_fold_changes.png",
    plot = p_fc,
    width = 8,
    height = 6,
    dpi = 300
)

################################################################################
# 6. Cell Count Summary
################################################################################

# Absolute cell counts
cell_counts <- meta_df %>%
    count(ID, Complete_Labels) %>%
    pivot_wider(names_from = ID, values_from = n, values_fill = 0)

write.csv(cell_counts,
          "./results/04_visualization/cell_counts_by_type.csv",
          row.names = FALSE)

# Clean up
rm(meta_df, prop_nh, prop_complete, prop_wide, contingency_table, 
   chisq_test, test_results, p_stacked, p_fc, cell_counts)
gc()

################################################################################
# End of Step 04
################################################################################
