################################################################################
# Single-Cell RNA-seq Analysis Pipeline
# Step 01: Data Loading and Quality Control
#
# Purpose: Load 10X Genomics scRNA-seq data and perform QC filtering
# Dataset: Mouse PORCN KO vs WT
#
# Input:
#   - ./data/wt/filtered_feature_bc_matrix/
#   - ./data/ko/filtered_feature_bc_matrix/
#
# Output:
#   - ./data/porcn.combined_QC.RData
#   - ./results/01_QC/*.png
#
# Author: YMS
# Date: 2025
################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(hdf5r)

# QC Parameters
QC_PARAMS <- list(
    min_features = 200,      # Minimum genes per cell
    max_features = 8000,     # Maximum genes per cell
    max_mt_pct   = 20        # Maximum mitochondrial percentage
)

# Create output directory
if (!dir.exists("./results/01_QC")) dir.create("./results/01_QC", recursive = TRUE)

################################################################################
# 1. Load 10X Data
################################################################################

porcn_wt <- Read10X(data.dir = "./data/wt/filtered_feature_bc_matrix/")
porcn_ko <- Read10X(data.dir = "./data/ko/filtered_feature_bc_matrix/")

################################################################################
# 2. Create Seurat Objects
################################################################################

# NOTE: Initial filtering applied at object creation
# - min.cells: Gene must be detected in ≥3 cells
# - min.features: Cell must have ≥200 genes
porcn_wt <- CreateSeuratObject(
    counts       = porcn_wt, 
    project      = 'porcn_wt',
    min.cells    = 3,
    min.features = 200
)

porcn_ko <- CreateSeuratObject(
    counts       = porcn_ko, 
    project      = 'porcn_ko',
    min.cells    = 3,
    min.features = 200
)

# Add sample ID metadata
porcn_wt$ID <- "WT"
porcn_ko$ID <- "KO"

################################################################################
# 3. Merge Objects
################################################################################

porcn.combined <- merge(porcn_wt, y = porcn_ko)
Idents(porcn.combined) <- 'ID'

################################################################################
# 4. Calculate QC Metrics
################################################################################

# NOTE: Mouse mitochondrial genes start with "mt-" (lowercase)
# Human would use "^MT-" (uppercase)
porcn.combined[["percent.mt"]] <- PercentageFeatureSet(
    porcn.combined, 
    pattern = "^mt-"
)

################################################################################
# 5. Generate Pre-Filtering QC Plots
################################################################################

# Violin plots
vln_plot <- VlnPlot(
    porcn.combined,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
) +
    plot_annotation(title = "QC Metrics Before Filtering")

ggsave(
    filename = "./results/01_QC/QC_violin_plots_before_filtering.png",
    plot     = vln_plot,
    width    = 12,
    height   = 5,
    dpi      = 300
)

# Scatter plots to examine metric relationships
plot1 <- FeatureScatter(
    porcn.combined,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
) + 
    ggtitle("UMI Count vs Mitochondrial %")

plot2 <- FeatureScatter(
    porcn.combined,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
) +
    ggtitle("UMI Count vs Gene Count")

ggsave(
    filename = "./results/01_QC/QC_scatter_plots.png",
    plot     = plot1 + plot2,
    width    = 12,
    height   = 5,
    dpi      = 300
)

################################################################################
# 6. Apply Quality Filters
################################################################################

# NOTE: These thresholds should be adjusted based on:
# - Tissue type
# - Cell type composition
# - Library quality
n_before <- ncol(porcn.combined)

porcn.combined <- subset(
    porcn.combined,
    subset = nFeature_RNA > QC_PARAMS$min_features & 
             nFeature_RNA < QC_PARAMS$max_features & 
             percent.mt < QC_PARAMS$max_mt_pct
)

n_after <- ncol(porcn.combined)

################################################################################
# 7. Generate Post-Filtering QC Plots
################################################################################

vln_plot_filtered <- VlnPlot(
    porcn.combined,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
) +
    plot_annotation(title = "QC Metrics After Filtering")

ggsave(
    filename = "./results/01_QC/QC_violin_plots_after_filtering.png",
    plot     = vln_plot_filtered,
    width    = 12,
    height   = 5,
    dpi      = 300
)

################################################################################
# 8. Save Processed Object
################################################################################

save(porcn.combined, file = "./data/porcn.combined_QC.RData")

# Clean up
rm(porcn_wt, porcn_ko, vln_plot, vln_plot_filtered, plot1, plot2)
gc()

################################################################################
# End of Step 01
################################################################################
